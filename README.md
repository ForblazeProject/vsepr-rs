# vsepr-rs

A lightweight, high-performance, and zero-dependency Rust crate for molecular geometry optimization based on **VSEPR (Valence Shell Electron Pair Repulsion)** theory.

## Features

- **Zero Dependency:** Built purely on Rust's standard library for maximum portability and fast compilation.
- **Generic Interface:** Uses traits (`AtomTrait`, `BondTrait`) so you can optimize your own data structures without conversion.
- **Fast Optimization:**
    - Uses a constraint-based relaxation engine (Force Field style).
    - **Spatial Hashing (O(N)):** Automatically switches to a grid-based spatial partition for large molecules or polymers (e.g., 100+ atoms) to maintain high performance.
- **Chemically Aware:**
    - Calculates Steric Numbers (SN) based on atomic numbers, bond orders, and formal charges.
    - Considers covalent radii for realistic bond lengths.
    - Prevents steric clashing and hydrogen overcrowding.
- **Strained Rings Support:** Correctly handles ring strain (e.g., cyclopropane) by balancing VSEPR ideal angles with distance constraints.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
vsepr-rs = "0.1.0"
```

## Quick Start

Implement `AtomTrait` and `BondTrait` for your structures and run the optimizer.

```rust
use vsepr_rs::{VseprOptimizer, AtomTrait, BondTrait};

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

fn main() {
    let mut atoms = vec![
        MyAtom { pos: [0.0, 0.0, 0.0], element: 8 }, // Oxygen
        MyAtom { pos: [1.0, 0.0, 0.0], element: 1 }, // Hydrogen
        MyAtom { pos: [0.0, 1.0, 0.0], element: 1 }, // Hydrogen
    ];
    let bonds = vec![
        MyBond { pair: (0, 1), order: 1.0 },
        MyBond { pair: (0, 2), order: 1.0 },
    ];

    let optimizer = VseprOptimizer::new();
    optimizer.optimize(&mut atoms, &bonds);

    println!("Optimized O-H1 position: {:?}", atoms[1].pos);
}
```

## Performance

`vsepr-rs` is designed to be fast enough for real-time applications and large-scale polymer research.

| Molecule | Atoms | Time (Release) |
| :--- | :--- | :--- |
| Methane ($CH_4$) | 5 | ~0.3 ms |
| Polyethylene ($C_{200}H_{402}$) | 602 | ~140 ms |

The O(N) spatial partitioning ensures that the calculation time scales linearly with the size of the molecule for large systems.

## License

MIT or Apache-2.0

## Author

**Forblaze Project**  
Website: [https://forblaze-works.com/](https://forblaze-works.com/)
