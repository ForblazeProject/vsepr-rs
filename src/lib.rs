pub mod math;
pub mod traits;
pub mod vsepr;

pub use math::Vec3;
pub use traits::{get_covalent_radius, AtomTrait, BondTrait};
pub use vsepr::{calculate_steric_number, Geometry};

/// Main optimizer for VSEPR-based molecular geometry refinement.
pub struct VseprOptimizer {
    /// Number of optimization iterations.
    pub iterations: usize,
    /// Movement scaling factor for each step.
    pub force_constant: f64,
}

impl Default for VseprOptimizer {
    fn default() -> Self {
        Self {
            iterations: 1500, // Increased iterations for torsion convergence
            force_constant: 0.15,
        }
    }
}

impl VseprOptimizer {
    pub fn new() -> Self {
        Self::default()
    }

    /// Refines the molecular geometry based on VSEPR theory and basic force field constraints.
    pub fn optimize<A: AtomTrait, B: BondTrait>(&self, atoms: &mut [A], bonds: &[B]) {
        let n = atoms.len();
        if n == 0 {
            return;
        }

        // 1. Preparation: Create an adjacency list for fast neighbor lookup.
        let mut adj = vec![Vec::new(); n];
        for (idx, bond) in bonds.iter().enumerate() {
            let (i, j) = bond.get_atom_indices();
            if i < n && j < n {
                adj[i].push((j, idx));
                adj[j].push((i, idx));
            }
        }

        // Extract coordinates and initialize with jitter if atoms are overlapping at the origin.
        let mut positions: Vec<Vec3> = atoms.iter().map(|a| Vec3::from(a.get_position())).collect();
        for i in 0..n {
            if positions[i].length_squared() < 0.001 {
                // Spread atoms if they are stacked at (0,0,0) to break symmetry.
                positions[i] = Vec3(
                    (i as f64 * 1.1).cos() * 1.0,
                    (i as f64 * 1.3).sin() * 1.0,
                    (i as f64 * 1.7).cos() * 1.0,
                );
            }
        }

        // Pre-calculate VSEPR geometries (Steric Numbers are topology-dependent).
        let geometries: Vec<Geometry> = (0..n)
            .map(|i| Geometry::from_steric_number(calculate_steric_number(i, atoms, bonds)))
            .collect();

        // 2. Optimization Loop
        for iter in 0..self.iterations {
            let mut forces = vec![Vec3::ZERO; n];
            // Exponential damping to ensure convergence in the final steps.
            let damping = (-3.0 * (iter as f64 / self.iterations as f64)).exp();

            // Apply various force field constraints
            self.apply_bond_constraints(atoms, bonds, &positions, &mut forces);
            self.apply_angle_constraints(n, &adj, &geometries, &positions, &mut forces);
            self.apply_planarity_constraints(n, &adj, &geometries, &positions, &mut forces);
            self.apply_torsion_constraints(n, &adj, bonds, &geometries, &positions, &mut forces); // New!
            self.apply_repulsion_constraints(n, &adj, &positions, &mut forces);

            // 3. Update coordinates
            for i in 0..n {
                positions[i] = positions[i].add(forces[i].mul(self.force_constant * damping));
            }
        }

        // 4. Write back results to user-defined structures.
        for i in 0..n {
            atoms[i].set_position(positions[i].into());
        }
    }

    /// Applies 1-2 interaction forces (Bond Lengths).
    fn apply_bond_constraints<A: AtomTrait, B: BondTrait>(
        &self,
        atoms: &[A],
        bonds: &[B],
        positions: &[Vec3],
        forces: &mut [Vec3],
    ) {
        for bond in bonds {
            let (i, j) = bond.get_atom_indices();
            let p1 = positions[i];
            let p2 = positions[j];
            let dist = p1.dist(p2);
            
            let r1 = get_covalent_radius(atoms[i].atomic_number());
            let r2 = get_covalent_radius(atoms[j].atomic_number());
            let raw_sum = (r1 + r2) / 100.0;

            let order = bond.get_bond_order();
            let scale_factor = if order >= 3.0 {
                0.82
            } else if order >= 2.0 {
                0.88
            } else if order >= 1.5 {
                0.92
            } else {
                1.00
            };

            let target = raw_sum * scale_factor;

            if dist > 0.001 {
                let diff = target - dist;
                let force = p1.sub(p2).normalize().mul(diff * 1.2);
                forces[i] = forces[i].add(force);
                forces[j] = forces[j].sub(force);
            }
        }
    }

    /// Applies 1-3 interaction forces (Bond Angles).
    fn apply_angle_constraints(
        &self,
        n: usize,
        adj: &[Vec<(usize, usize)>],
        geometries: &[Geometry],
        positions: &[Vec3],
        forces: &mut [Vec3],
    ) {
        for i in 0..n {
            let neighbors = &adj[i];
            if neighbors.len() < 2 {
                continue;
            }

            let ideal_angle = geometries[i].ideal_angle();
            let cos_theta = ideal_angle.cos();

            for a_idx in 0..neighbors.len() {
                for b_idx in (a_idx + 1)..neighbors.len() {
                    let ni = neighbors[a_idx].0;
                    let nj = neighbors[b_idx].0;

                    let p_i = positions[ni];
                    let p_j = positions[nj];
                    let p_center = positions[i];

                    let v_i = p_i.sub(p_center);
                    let v_j = p_j.sub(p_center);
                    let r1 = v_i.length();
                    let r2 = v_j.length();

                    let target_d = (r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * cos_theta)
                        .max(0.0)
                        .sqrt();

                    let current_d = p_i.dist(p_j);
                    if current_d > 0.001 {
                        let diff = target_d - current_d;
                        let force = p_i.sub(p_j).normalize().mul(diff * 0.2);
                        forces[ni] = forces[ni].add(force);
                        forces[nj] = forces[nj].sub(force);
                    }
                }
            }
        }
    }

    /// Applies Planarity Constraints (Improper Torsions) for sp2 centers.
    fn apply_planarity_constraints(
        &self,
        n: usize,
        adj: &[Vec<(usize, usize)>],
        geometries: &[Geometry],
        positions: &[Vec3],
        forces: &mut [Vec3],
    ) {
        for i in 0..n {
            if let Geometry::TrigonalPlanar = geometries[i] {
                let neighbors = &adj[i];
                if neighbors.len() == 3 {
                    let p0 = positions[i];
                    let p1 = positions[neighbors[0].0];
                    let p2 = positions[neighbors[1].0];
                    let p3 = positions[neighbors[2].0];

                    let v12 = p2.sub(p1);
                    let v13 = p3.sub(p1);
                    let normal = v12.cross(v13).normalize();

                    if normal.length_squared() < 1e-6 {
                        continue;
                    }

                    let v10 = p0.sub(p1);
                    let dist_from_plane = v10.dot(normal);

                    let correction = normal.mul(-dist_from_plane * 0.5);
                    
                    forces[i] = forces[i].add(correction);
                    
                    let recoil = correction.mul(-0.333);
                    forces[neighbors[0].0] = forces[neighbors[0].0].add(recoil);
                    forces[neighbors[1].0] = forces[neighbors[1].0].add(recoil);
                    forces[neighbors[2].0] = forces[neighbors[2].0].add(recoil);
                }
            }
        }
    }

    /// Applies Dihedral (Torsion) Constraints (1-4 interactions).
    /// Focuses on keeping aromatic and double bonds planar (0 or 180 degrees).
    fn apply_torsion_constraints<B: BondTrait>(
        &self,
        _n: usize,
        adj: &[Vec<(usize, usize)>],
        bonds: &[B],
        geometries: &[Geometry],
        positions: &[Vec3],
        forces: &mut [Vec3],
    ) {
        for bond in bonds {
            let (j, k) = bond.get_atom_indices();
            let order = bond.get_bond_order();

            // Only apply strong torsion constraints for Aromatic (1.5) and Double (2.0) bonds.
            if order < 1.5 { continue; }
            // Both atoms must be sp2 (TrigonalPlanar) for planar constraint.
            if !matches!(geometries[j], Geometry::TrigonalPlanar) || 
               !matches!(geometries[k], Geometry::TrigonalPlanar) {
                continue;
            }

            let neighbors_j = &adj[j];
            let neighbors_k = &adj[k];

            for &(i, _) in neighbors_j {
                if i == k { continue; } // Skip the bond itself
                for &(l, _) in neighbors_k {
                    if l == j { continue; } // Skip the bond itself

                    // Torsion i-j-k-l
                    let p_i = positions[i];
                    let p_j = positions[j];
                    let p_k = positions[k];
                    let p_l = positions[l];

                    // Vectors: b1=i->j, b2=j->k, b3=k->l
                    let b1 = p_j.sub(p_i);
                    let b2 = p_k.sub(p_j);
                    let b3 = p_l.sub(p_k);

                    // Normal vectors to planes defined by (i,j,k) and (j,k,l)
                    let n1 = b1.cross(b2);
                    let n2 = b2.cross(b3);
                    
                    if n1.length_squared() < 1e-6 || n2.length_squared() < 1e-6 { continue; }

                    let n1_norm = n1.normalize();
                    let n2_norm = n2.normalize();

                    // Cosine of dihedral angle
                    let _cos_phi = n1_norm.dot(n2_norm);
                    // Sine of dihedral angle (roughly, using cross product direction)
                    let sin_phi = n1_norm.cross(n2_norm).dot(b2.normalize());

                    // Target: Planar (0 or 180 degrees)
                    // We want to minimize sin_phi (make planes parallel).
                    // Force is proportional to sin(phi).
                    
                    // Simple restoring force towards closest planar configuration (0 or 180).
                    // If sin_phi is positive, torque direction is opposite.
                    
                    let strength = 0.05; // Torsion is weaker than bonds/angles
                    let torque = sin_phi * strength;

                    // Apply rotational forces (simplified perpendicular forces)
                    // Push i and l towards the plane defined by j-k bond and ideal vector.
                    
                    // Direction perpendicular to the planes
                    let f_i = n1_norm.mul(-torque);
                    let f_l = n2_norm.mul(torque);

                    forces[i] = forces[i].add(f_i);
                    forces[l] = forces[l].add(f_l);
                    
                    // Recoil on central atoms (simplified conservation)
                    forces[j] = forces[j].sub(f_i);
                    forces[k] = forces[k].sub(f_l);
                }
            }
        }
    }

    /// Applies Non-bonded Repulsion (Van der Waals-like).
    fn apply_repulsion_constraints(
        &self,
        n: usize,
        adj: &[Vec<(usize, usize)>],
        positions: &[Vec3],
        forces: &mut [Vec3],
    ) {
        let min_repulsion_d = 1.2;
        if n < 100 {
            for i in 0..n {
                for j in (i + 1)..n {
                    if adj[i].iter().any(|&(neighbor, _)| neighbor == j) {
                        continue;
                    }
                    self.apply_repulsion(i, j, positions, forces, min_repulsion_d);
                }
            }
        } else {
            self.apply_spatial_repulsion(n, positions, adj, forces, min_repulsion_d);
        }
    }

    #[inline]
    fn apply_repulsion(
        &self,
        i: usize,
        j: usize,
        positions: &[Vec3],
        forces: &mut [Vec3],
        min_d: f64,
    ) {
        let p1 = positions[i];
        let p2 = positions[j];
        let diff_v = p1.sub(p2);
        let dist_sq = diff_v.length_squared();
        
        if dist_sq < min_d * min_d && dist_sq > 0.0001 {
            let dist = dist_sq.sqrt();
            let force_mag = (min_d - dist).powi(2) * 0.8;
            let force = diff_v.normalize().mul(force_mag);
            forces[i] = forces[i].add(force);
            forces[j] = forces[j].sub(force);
        }
    }

    fn apply_spatial_repulsion(
        &self,
        _n: usize,
        positions: &[Vec3],
        adj: &[Vec<(usize, usize)>],
        forces: &mut [Vec3],
        min_d: f64,
    ) {
        let mut min_p = positions[0];
        let mut max_p = positions[0];
        for p in positions.iter().skip(1) {
            min_p = Vec3(min_p.0.min(p.0), min_p.1.min(p.1), min_p.2.min(p.2));
            max_p = Vec3(max_p.0.max(p.0), max_p.1.max(p.1), max_p.2.max(p.2));
        }

        let cell_size = min_d;
        let dim_x = ((max_p.0 - min_p.0) / cell_size).ceil() as usize + 1;
        let dim_y = ((max_p.1 - min_p.1) / cell_size).ceil() as usize + 1;
        let dim_z = ((max_p.2 - min_p.2) / cell_size).ceil() as usize + 1;

        let get_cell_idx = |p: Vec3| -> (usize, usize, usize) {
            (
                ((p.0 - min_p.0) / cell_size) as usize,
                ((p.1 - min_p.1) / cell_size) as usize,
                ((p.2 - min_p.2) / cell_size) as usize,
            )
        };

        let mut cells = vec![Vec::new(); dim_x * dim_y * dim_z];
        for (i, &p) in positions.iter().enumerate() {
            let (ix, iy, iz) = get_cell_idx(p);
            let idx = ix * (dim_y * dim_z) + iy * dim_z + iz;
            cells[idx].push(i);
        }

        for ix in 0..dim_x {
            for iy in 0..dim_y {
                for iz in 0..dim_z {
                    let cell_idx = ix * (dim_y * dim_z) + iy * dim_z + iz;
                    let cell_atoms = &cells[cell_idx];
                    if cell_atoms.is_empty() {
                        continue;
                    }

                    for dx in -1..=1 {
                        for dy in -1..=1 {
                            for dz in -1..=1 {
                                let nix = ix as i32 + dx;
                                let niy = iy as i32 + dy;
                                let niz = iz as i32 + dz;

                                if nix >= 0
                                    && nix < dim_x as i32
                                    && niy >= 0
                                    && niy < dim_y as i32
                                    && niz >= 0
                                    && niz < dim_z as i32
                                {
                                    let neighbor_idx = (nix as usize) * (dim_y * dim_z)
                                        + (niy as usize) * dim_z
                                        + (niz as usize);
                                    let neighbor_atoms = &cells[neighbor_idx];

                                    for &i in cell_atoms {
                                        for &j in neighbor_atoms {
                                            if i >= j {
                                                continue;
                                            }
                                            if adj[i].iter().any(|&(neighbor, _)| neighbor == j) {
                                                continue;
                                            }
                                            self.apply_repulsion(i, j, positions, forces, min_d);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MyAtom {
        pos: [f64; 3],
        element: usize,
    }

    impl AtomTrait for MyAtom {
        fn get_position(&self) -> [f64; 3] { self.pos }
        fn set_position(&mut self, pos: [f64; 3]) { self.pos = pos; }
        fn atomic_number(&self) -> usize { self.element }
    }

    struct MyBond {
        pair: (usize, usize),
        order: f32,
    }

    impl BondTrait for MyBond {
        fn get_atom_indices(&self) -> (usize, usize) { self.pair }
        fn get_bond_order(&self) -> f32 { self.order }
    }

    #[test]
    fn test_water_geometry() {
        let mut atoms = vec![
            MyAtom { pos: [0.0, 0.0, 0.0], element: 8 },
            MyAtom { pos: [1.0, 0.0, 0.0], element: 1 },
            MyAtom { pos: [0.0, 1.0, 0.0], element: 1 },
        ];
        let bonds = vec![
            MyBond { pair: (0, 1), order: 1.0 },
            MyBond { pair: (0, 2), order: 1.0 },
        ];
        VseprOptimizer::new().optimize(&mut atoms, &bonds);
        let v1 = Vec3::from(atoms[1].pos).sub(atoms[0].pos.into()).normalize();
        let v2 = Vec3::from(atoms[2].pos).sub(atoms[0].pos.into()).normalize();
        let angle = v1.dot(v2).acos().to_degrees();
        assert!(angle > 100.0 && angle < 115.0);
    }

    #[test]
    fn test_benzene_geometry_strict() {
        // Strict Benzene Test: C6H6
        // Initial state is highly distorted (zigzag / puckered)
        let mut atoms = Vec::new();
        for i in 0..6 {
            // Highly distorted initial coordinates
            atoms.push(MyAtom {
                pos: [
                    (i as f64).cos() * 1.4, 
                    (i as f64).sin() * 1.4, 
                    (i % 2) as f64 * 0.8 // Large z-displacement (0.8 A)
                ],
                element: 6
            });
        }
        // Hydrogens attached with random z-offsets
        for i in 0..6 {
             atoms.push(MyAtom { 
                 pos: [
                     (i as f64).cos() * 2.4, 
                     (i as f64).sin() * 2.4, 
                     ((i + 1) % 2) as f64 * -0.5 
                 ], 
                 element: 1 
             });
        }

        let mut bonds = Vec::new();
        // C-C aromatic bonds
        for i in 0..6 {
            bonds.push(MyBond { pair: (i, (i + 1) % 6), order: 1.5 });
        }
        // C-H bonds
        for i in 0..6 {
            bonds.push(MyBond { pair: (i, 6 + i), order: 1.0 });
        }

        // Run optimization with increased iterations for torsion convergence
        let optimizer = VseprOptimizer { iterations: 2000, force_constant: 0.1 };
        optimizer.optimize(&mut atoms, &bonds);

        // 1. Check C-C bond lengths (Target ~1.40 A)
        let p0 = Vec3::from(atoms[0].pos);
        let p1 = Vec3::from(atoms[1].pos);
        let dist = p0.dist(p1);
        println!("Benzene C-C dist: {:.3}", dist);
        assert!((dist - 1.40).abs() < 0.05);

        // 2. Check Planarity (Strict)
        // Fit a plane to all 6 Carbon atoms and check max deviation.
        // Simple check: Normal vector of triangle (0,2,4) vs (1,3,5) should be parallel,
        // and distance of atoms to the average plane should be small.
        
        // Let's check average z-deviation after aligning to principal axes is hard in simple test.
        // Instead, check local torsion angles or point-plane distance.
        // Here we use the cross product method for all consecutive 4 atoms.
        
        let mut max_deviation = 0.0;
        for i in 0..6 {
            let p_i = Vec3::from(atoms[i].pos); // C
            let p_j = Vec3::from(atoms[(i+1)%6].pos); // C
            let p_k = Vec3::from(atoms[(i+2)%6].pos); // C
            let p_l = Vec3::from(atoms[(i+3)%6].pos); // C
            
            // Dihedral i-j-k-l
            let b1 = p_j.sub(p_i);
            let b2 = p_k.sub(p_j);
            let b3 = p_l.sub(p_k);
            let n1 = b1.cross(b2).normalize();
            let n2 = b2.cross(b3).normalize();
            let sin_phi = n1.cross(n2).length(); // If planar, sin_phi ~ 0
            
            if sin_phi > max_deviation {
                max_deviation = sin_phi;
            }
        }
        
        println!("Max dihedral deviation (sin_phi): {:.4}", max_deviation);
        // Expect extremely flat structure (deviation < 0.1 radians ~ 5.7 degrees)
        assert!(max_deviation < 0.1);
    }
}
