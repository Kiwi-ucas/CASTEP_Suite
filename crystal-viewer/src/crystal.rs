//! Crystal structure data types and JSON parsing

use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct CrystalData {
    pub lattice: Lattice,
    pub atoms: Vec<AtomData>,
    #[serde(default)]
    pub positions_fractional: bool,
}

#[derive(Debug, Deserialize)]
pub struct Lattice {
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub alpha: f32,
    pub beta: f32,
    pub gamma: f32,
}

impl Lattice {
    /// Compute Cartesian lattice vectors from a, b, c, alpha, beta, gamma
    pub fn to_vectors(&self) -> [Vec3; 3] {
        let a = self.a;
        let b = self.b;
        let c = self.c;
        let alpha = self.alpha.to_radians();
        let beta = self.beta.to_radians();
        let gamma = self.gamma.to_radians();

        let v1 = Vec3::new(a, 0.0, 0.0);
        let v2 = Vec3::new(b * gamma.cos(), b * gamma.sin(), 0.0);
        let v3_x = c * beta.cos();
        let v3_y = c * (alpha.cos() - beta.cos() * gamma.cos()) / gamma.sin();
        let v3_z = (c * c - v3_x * v3_x - v3_y * v3_y).sqrt();
        let v3 = Vec3::new(v3_x, v3_y, v3_z.max(0.0));

        [v1, v2, v3]
    }

    pub fn cell_volume(&self) -> f32 {
        let a = self.a; let b = self.b; let c = self.c;
        let alpha = self.alpha.to_radians();
        let beta = self.beta.to_radians();
        let gamma = self.gamma.to_radians();
        let cos_a = alpha.cos(); let cos_b = beta.cos(); let cos_g = gamma.cos();
        a * b * c * (1.0 - cos_a * cos_a - cos_b * cos_b - cos_g * cos_g
            + 2.0 * cos_a * cos_b * cos_g).sqrt()
    }
}

#[derive(Debug, Deserialize)]
pub struct AtomData {
    pub element: String,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    #[serde(default)]
    pub label: String,
}

impl CrystalData {
    pub fn from_json(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        let data: CrystalData = serde_json::from_str(&content)?;
        Ok(data)
    }

    /// Get Cartesian positions for all atoms
    pub fn cartesian_positions(&self) -> Vec<Vec3> {
        let vecs = self.lattice.to_vectors();
        self.atoms.iter().map(|atom| {
            if self.positions_fractional {
                atom.x * vecs[0] + atom.y * vecs[1] + atom.z * vecs[2]
            } else {
                Vec3::new(atom.x, atom.y, atom.z)
            }
        }).collect()
    }

    /// Get the cell corner vertices (8 corners)
    pub fn cell_corners(&self) -> [Vec3; 8] {
        let [a, b, c] = self.lattice.to_vectors();
        let zero = Vec3::ZERO;
        [
            zero,
            a,
            a + b,
            b,
            c,
            a + c,
            a + b + c,
            b + c,
        ]
    }

    /// Get the 12 edge index pairs
    pub fn cell_edges() -> [(usize, usize); 12] {
        [
            (0, 1), (1, 2), (2, 3), (3, 0), // bottom
            (4, 5), (5, 6), (6, 7), (7, 4), // top
            (0, 4), (1, 5), (2, 6), (3, 7), // vertical
        ]
    }

    /// Compute the center of the structure
    pub fn center(&self) -> Vec3 {
        let corners = self.cell_corners();
        let mut sum = Vec3::ZERO;
        for c in corners.iter() {
            sum += *c;
        }
        sum / 8.0
    }
}

use bevy::prelude::Vec3;
