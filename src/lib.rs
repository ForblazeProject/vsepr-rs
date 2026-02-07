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

            // --- Bond Constraints (1-2 interactions) ---
            // Pulls bonded atoms toward their ideal covalent bond lengths.
            for bond in bonds {
                let (i, j) = bond.get_atom_indices();
                let p1 = positions[i];
                let p2 = positions[j];
                let dist = p1.dist(p2);
                let target = (get_covalent_radius(atoms[i].atomic_number())
                    + get_covalent_radius(atoms[j].atomic_number()))
                    / 100.0; // Convert pm to Angstroms

                if dist > 0.001 {
                    let diff = target - dist;
                    let force = p1.sub(p2).normalize().mul(diff * 1.0);
                    forces[i] = forces[i].add(force);
                    forces[j] = forces[j].sub(force);
                }
            }

            // --- VSEPR Angle Constraints (1-3 interactions) ---
            // Adjusts angles between two bonds sharing a common central atom.
            for i in 0..n {
                let neighbors = &adj[i];
                if neighbors.len() < 2 {
                    continue;
                }

                let ideal_angle = geometries[i].ideal_angle();

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

                        // Law of cosines to find the ideal 1-3 distance.
                        let target_d = (r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * ideal_angle.cos())
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

            // --- Non-bonded Repulsion (Van der Waals-like) ---
            // Prevents steric clashing, especially for Hydrogens.
            let min_repulsion_d = 1.2;
            if n < 100 {
                // Small systems: Use simple O(N^2) loop.
                for i in 0..n {
                    for j in (i + 1)..n {
                        if adj[i].iter().any(|&(neighbor, _)| neighbor == j) {
                            continue;
                        }
                        self.apply_repulsion(i, j, &positions, &mut forces, min_repulsion_d);
                    }
                }
            } else {
                // Large systems (Polymers, etc.): Use Spatial Hashing (O(N)).
                self.apply_spatial_repulsion(n, &positions, &adj, &mut forces, min_repulsion_d);
            }

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
            let force = diff_v.normalize().mul((min_d - dist).powi(2) * 0.3);
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
    fn test_methane_geometry() {
        let mut atoms = vec![
            MyAtom { pos: [0.0, 0.0, 0.0], element: 6 },
            MyAtom { pos: [0.1, 0.1, 0.1], element: 1 },
            MyAtom { pos: [-0.1, 0.1, 0.1], element: 1 },
            MyAtom { pos: [0.1, -0.1, 0.1], element: 1 },
            MyAtom { pos: [0.1, 0.1, -0.1], element: 1 },
        ];
        let bonds = vec![
            MyBond { pair: (0, 1), order: 1.0 },
            MyBond { pair: (0, 2), order: 1.0 },
            MyBond { pair: (0, 3), order: 1.0 },
            MyBond { pair: (0, 4), order: 1.0 },
        ];
        VseprOptimizer::new().optimize(&mut atoms, &bonds);
        let p_c: Vec3 = atoms[0].pos.into();
        for i in 1..=4 {
            for j in (i+1)..=4 {
                let v_i = Vec3::from(atoms[i].pos).sub(p_c).normalize();
                let v_j = Vec3::from(atoms[j].pos).sub(p_c).normalize();
                let angle = v_i.dot(v_j).acos().to_degrees();
                assert!(angle > 105.0 && angle < 115.0);
            }
        }
    }

    #[test]
    fn test_polyethylene_performance_and_geometry() {
        use std::time::Instant;
        let n_c = 200;
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();

        for i in 0..n_c {
            atoms.push(MyAtom { 
                pos: [i as f64 * 0.5, (i % 2) as f64 * 0.1, (i / 2 % 2) as f64 * 0.1], 
                element: 6 
            });
        }
        for i in 0..n_c - 1 {
            bonds.push(MyBond { pair: (i, i+1), order: 1.0 });
        }
        for i in 0..n_c {
            let h_count = if i == 0 || i == n_c - 1 { 3 } else { 2 };
            let c_pos = atoms[i].pos;
            for _ in 0..h_count {
                let h_idx = atoms.len();
                atoms.push(MyAtom { pos: [c_pos[0] + 0.1, c_pos[1], c_pos[2]], element: 1 });
                bonds.push(MyBond { pair: (i, h_idx), order: 1.0 });
            }
        }

        let start_pe = Instant::now();
        VseprOptimizer::new().optimize(&mut atoms, &bonds);
        let dur_pe = start_pe.elapsed();

        println!("Polyethylene ({} atoms) with Spatial Hash: {:?}", atoms.len(), dur_pe);

        let p_start = Vec3::from(atoms[0].pos);
        let p_end = Vec3::from(atoms[n_c - 1].pos);
        let chain_dist = p_start.dist(p_end);
        println!("Chain end-to-end distance: {:.2} A", chain_dist);

        assert!(chain_dist > 50.0);
    }
}
