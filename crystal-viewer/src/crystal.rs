//! Crystal structure data types and JSON parsing

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CrystalData {
    pub lattice: Lattice,
    pub atoms: Vec<AtomData>,
    #[serde(default)]
    pub positions_fractional: bool,
    #[serde(default)]
    pub modified: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
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
        let v3_z = (c * c - v3_x * v3_x - v3_y * v3_y).max(0.0).sqrt();
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

    /// Compute inverse of the 3×3 Cartesian lattice matrix (columns are a, b, c vectors)
    pub fn inverse_vectors(&self) -> [Vec3; 3] {
        let [a, b, c] = self.to_vectors();
        let m = [
            [a.x, b.x, c.x],
            [a.y, b.y, c.y],
            [a.z, b.z, c.z],
        ];
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if det.abs() < 1e-6 {
            return [Vec3::ZERO, Vec3::ZERO, Vec3::ZERO];
        }
        let inv_det = 1.0 / det;
        let col0 = Vec3::new(
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det,
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
        );
        let col1 = Vec3::new(
            -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det,
        );
        let col2 = Vec3::new(
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
            -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        );
        [col0, col1, col2]
    }

    /// Multiply inverse matrix by vector: M⁻¹ * v (rows of M⁻¹ dot v)
    pub fn apply_inverse(inv: &[Vec3; 3], v: Vec3) -> Vec3 {
        Vec3::new(
            inv[0].x * v.x + inv[1].x * v.y + inv[2].x * v.z,
            inv[0].y * v.x + inv[1].y * v.y + inv[2].y * v.z,
            inv[0].z * v.x + inv[1].z * v.y + inv[2].z * v.z,
        )
    }

    /// Convert Cartesian coordinates to fractional: frac = M⁻¹ * cart
    pub fn to_fractional(&self, cart: Vec3) -> Vec3 {
        Self::apply_inverse(&self.inverse_vectors(), cart)
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
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

impl CrystalData {
    /// Expand asymmetric unit to full unit cell by applying ±1 translations.
    /// Returns (Cartesian positions, parent indices, fractional offsets).
    pub fn expand_to_cell(&self) -> (Vec<Vec3>, Vec<usize>, Vec<Vec3>) {
        let vecs = self.lattice.to_vectors();
        let inv = self.lattice.inverse_vectors();
        let mut positions = Vec::new();
        let mut parents = Vec::new();
        let mut offsets = Vec::new();

        for (i, atom) in self.atoms.iter().enumerate() {
            let frac = if self.positions_fractional {
                Vec3::new(atom.x, atom.y, atom.z)
            } else {
                let cart = Vec3::new(atom.x, atom.y, atom.z);
                Lattice::apply_inverse(&inv, cart)
            };

            for di in -1..=1_i32 {
                for dj in -1_i32..=1 {
                    for dk in -1_i32..=1 {
                        let offset = Vec3::new(di as f32, dj as f32, dk as f32);
                        let f = frac + offset;
                        if f.x >= -0.001 && f.x < 1.001
                            && f.y >= -0.001 && f.y < 1.001
                            && f.z >= -0.001 && f.z < 1.001
                        {
                            let cart = f.x * vecs[0] + f.y * vecs[1] + f.z * vecs[2];
                            positions.push(cart);
                            parents.push(i);
                            offsets.push(offset);
                        }
                    }
                }
            }
        }

        (positions, parents, offsets)
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, Box<dyn std::error::Error>> {
        let json = serde_json::to_string_pretty(self)?;
        Ok(json)
    }

    /// Write modified data back to JSON file
    pub fn write_to_file(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let json = self.to_json()?;
        std::fs::write(path, json)?;
        Ok(())
    }
}

use bevy::prelude::Vec3;
