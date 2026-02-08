//! vsepr-rs: A lightweight, high-performance molecular geometry optimizer.
//!
//! This crate provides a generic engine to refine 3D molecular coordinates 
//! based on VSEPR (Valence Shell Electron Pair Repulsion) theory and a 
//! lightweight force field.

pub mod forcefield;
pub mod math;
pub mod optimizer;
pub mod traits;
pub mod vsepr;

pub use math::Vec3;
pub use optimizer::VseprOptimizer;
pub use traits::{get_covalent_radius, AtomTrait, BondTrait};
pub use vsepr::{calculate_steric_number, Geometry};

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
    fn test_benzene_from_origin() {
        let mut atoms = Vec::new();
        // 1. Initialize all atoms at origin
        for _ in 0..6 {
            atoms.push(MyAtom { pos: [0.0, 0.0, 0.0], element: 6 });
        }
        for _ in 0..6 {
             atoms.push(MyAtom { pos: [0.0, 0.0, 0.0], element: 1 });
        }
        run_benzene_check(atoms, "Origin Initialization");
    }

    #[test]
    fn test_benzene_from_linear() {
        let mut atoms = Vec::new();
        // 2. Initialize all atoms in a straight line
        for i in 0..6 {
            atoms.push(MyAtom { pos: [i as f64 * 1.5, 0.0, 0.0], element: 6 });
        }
        for i in 0..6 {
             atoms.push(MyAtom { pos: [(i as f64 * 1.5) + 0.5, 0.0, 0.0], element: 1 });
        }
        run_benzene_check(atoms, "Linear Initialization");
    }

    fn run_benzene_check(mut atoms: Vec<MyAtom>, test_name: &str) {
        let mut bonds = Vec::new();
        for i in 0..6 {
            bonds.push(MyBond { pair: (i, (i + 1) % 6), order: 1.5 });
        }
        for i in 0..6 {
            bonds.push(MyBond { pair: (i, 6 + i), order: 1.0 });
        }

        let optimizer = VseprOptimizer { iterations: 3000, force_constant: 0.1 };
        optimizer.optimize(&mut atoms, &bonds);

        println!("--- {} ---", test_name);
        let p1 = Vec3::from(atoms[0].pos);
        let p2 = Vec3::from(atoms[1].pos);
        let dist = p1.dist(p2);
        println!("C-C dist: {:.3}", dist);
        assert!((dist - 1.40).abs() < 0.15, "Bond length check failed");
    }
}
