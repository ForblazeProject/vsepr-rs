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
            iterations: 1000,
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
            
            // Calculate ideal bond length based on covalent radii.
            let r1 = get_covalent_radius(atoms[i].atomic_number());
            let r2 = get_covalent_radius(atoms[j].atomic_number());
            let raw_sum = (r1 + r2) / 100.0; // Convert pm to Angstroms

            // Scale based on bond order (bond order correction).
            let order = bond.get_bond_order();
            let scale_factor = if order >= 3.0 {
                0.82 // Triple bond
            } else if order >= 2.0 {
                0.88 // Double bond
            } else if order >= 1.5 {
                0.92 // Aromatic bond
            } else {
                1.00 // Single bond
            };

            let target = raw_sum * scale_factor;

            if dist > 0.001 {
                let diff = target - dist;
                // Strong spring constant for bonds
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
            // Pre-calculate cos of ideal angle
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

                    // Law of cosines to find the ideal 1-3 distance based on CURRENT bond lengths.
                    // Note: While user noted this fluctuates, it's a standard approximation in simple relaxation.
                    // Using fixed ideal bond lengths here might conflict with actual bond constraints.
                    let target_d = (r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * cos_theta)
                        .max(0.0)
                        .sqrt();

                    let current_d = p_i.dist(p_j);
                    if current_d > 0.001 {
                        let diff = target_d - current_d;
                        // Angle forces are weaker than bonds
                        let force = p_i.sub(p_j).normalize().mul(diff * 0.2);
                        forces[ni] = forces[ni].add(force);
                        forces[nj] = forces[nj].sub(force);
                    }
                }
            }
        }
    }

    /// Applies Planarity Constraints (Improper Torsions) for sp2 centers.
    /// This is crucial for aromatic rings (Benzene) and double bonds.
    fn apply_planarity_constraints(
        &self,
        n: usize,
        adj: &[Vec<(usize, usize)>],
        geometries: &[Geometry],
        positions: &[Vec3],
        forces: &mut [Vec3],
    ) {
        for i in 0..n {
            // Only apply to Trigonal Planar (SN=3) atoms with exactly 3 neighbors.
            if let Geometry::TrigonalPlanar = geometries[i] {
                let neighbors = &adj[i];
                if neighbors.len() == 3 {
                    let p0 = positions[i];
                    let p1 = positions[neighbors[0].0];
                    let p2 = positions[neighbors[1].0];
                    let p3 = positions[neighbors[2].0];

                    // Define the plane using the three neighbors (p1, p2, p3).
                    // Normal vector N = (p2 - p1) x (p3 - p1)
                    let v12 = p2.sub(p1);
                    let v13 = p3.sub(p1);
                    let normal = v12.cross(v13).normalize();

                    // If normal is zero length (collinear neighbors), skip.
                    if normal.length_squared() < 1e-6 {
                        continue;
                    }

                    // Calculate distance of central atom p0 from the plane.
                    // Distance d = dot(p0 - p1, normal)
                    let v10 = p0.sub(p1);
                    let dist_from_plane = v10.dot(normal);

                    // Force to push p0 back onto the plane.
                    // p0 moves along -normal * dist
                    // Neighbors move along +normal * (dist / 3) to conserve momentum roughly.
                    let correction = normal.mul(-dist_from_plane * 0.5); // 0.5 is stiffness
                    
                    forces[i] = forces[i].add(correction);
                    
                    let recoil = correction.mul(-0.333);
                    forces[neighbors[0].0] = forces[neighbors[0].0].add(recoil);
                    forces[neighbors[1].0] = forces[neighbors[1].0].add(recoil);
                    forces[neighbors[2].0] = forces[neighbors[2].0].add(recoil);
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
            // Small systems: Use simple O(N^2) loop.
            for i in 0..n {
                for j in (i + 1)..n {
                    if adj[i].iter().any(|&(neighbor, _)| neighbor == j) {
                        continue;
                    }
                    self.apply_repulsion(i, j, positions, forces, min_repulsion_d);
                }
            }
        } else {
            // Large systems (Polymers, etc.): Use Spatial Hashing (O(N)).
            self.apply_spatial_repulsion(n, positions, adj, forces, min_repulsion_d);
        }
    }

    /// Calculates and applies repulsion force between two atoms if they are too close.
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
            // Increased repulsion stiffness from 0.3 to 0.8
            let force_mag = (min_d - dist).powi(2) * 0.8;
            let force = diff_v.normalize().mul(force_mag);
            forces[i] = forces[i].add(force);
            forces[j] = forces[j].sub(force);
        }
    }

    /// High-performance repulsion calculation using a grid-based spatial partition.
    fn apply_spatial_repulsion(
        &self,
        _n: usize,
        positions: &[Vec3],
        adj: &[Vec<(usize, usize)>],
        forces: &mut [Vec3],
        min_d: f64,
    ) {
        // Compute bounding box to define grid dimensions.
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

        // Populate the grid with atom indices.
        let mut cells = vec![Vec::new(); dim_x * dim_y * dim_z];
        for (i, &p) in positions.iter().enumerate() {
            let (ix, iy, iz) = get_cell_idx(p);
            let idx = ix * (dim_y * dim_z) + iy * dim_z + iz;
            cells[idx].push(i);
        }

        // Search neighboring 27 cells (self + neighbors) for potential clashes.
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
                                            } // Avoid double counting and self-repulsion.
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
    fn test_benzene_geometry() {
        // Benzene: C6H6. Ring of 6 carbons.
        // Simplified test without Hydrogens for basic ring planarity and distance.
        let mut atoms = Vec::new();
        // Initialize randomly or in a non-planar way
        for i in 0..6 {
            atoms.push(MyAtom {
                pos: [
                    (i as f64).cos() * 1.4 + (i % 2) as f64 * 0.2, 
                    (i as f64).sin() * 1.4, 
                    (i % 2) as f64 * 0.5 // Zig-zag (puckered) initial state
                ],
                element: 6
            });
        }
        let mut bonds = Vec::new();
        for i in 0..6 {
            bonds.push(MyBond { pair: (i, (i + 1) % 6), order: 1.5 }); // Aromatic bond
        }

        // Add dummy Hydrogens to ensure C is sp2 (SN=3).
        // Each C needs 1 H.
        for i in 0..6 {
             atoms.push(MyAtom { pos: [0.0, 0.0, 0.0], element: 1 });
             bonds.push(MyBond { pair: (i, 6 + i), order: 1.0 });
        }

        VseprOptimizer::new().optimize(&mut atoms, &bonds);

        // Check C-C bond lengths
        let p0 = Vec3::from(atoms[0].pos);
        let p1 = Vec3::from(atoms[1].pos);
        let dist = p0.dist(p1);
        println!("Benzene C-C dist: {:.3}", dist);
        // Ideal aromatic C-C is ~1.40 A. Normal single is 1.54.
        // Our logic: (0.76+0.76) * 0.92 = 1.398. Should be close to 1.40.
        assert!((dist - 1.40).abs() < 0.05);

        // Check Planarity
        // Check if z-coordinates are roughly similar (after optimization, it might rotate, 
        // so checking dot product of normals is better, or just visual inspection via debug print).
        // For simplicity in this test, we check if C0, C1, C2, C3 form a flat plane.
        let p2 = Vec3::from(atoms[2].pos);
        let p3 = Vec3::from(atoms[3].pos);
        // Torsion angle or just distance of p3 from plane(p0, p1, p2)
        let v01 = p1.sub(p0);
        let v12 = p2.sub(p1);
        let normal = v01.cross(v12).normalize();
        let v23 = p3.sub(p2);
        let dist_from_plane = v23.dot(normal).abs();
        println!("Benzene flatness error: {:.3}", dist_from_plane);
        assert!(dist_from_plane < 0.1);
    }
}