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
    #[serde(default)]
    pub phonon_modes: Option<PhononModesData>,
    /// Slab cross-section provenance (set by the slab builder). Unknown keys
    /// are ignored by the Fortran JSON reader — metadata only.
    #[serde(default)]
    pub slab: Option<SlabMeta>,
    /// Vacuum-layer provenance (set by the vacuum builder).
    #[serde(default)]
    pub vacuum: Option<VacuumMeta>,
    /// Supercell provenance (set by the supercell builder): the last
    /// applied integer expansion (nx, ny, nz) along a/b/c. Metadata only
    /// — the Fortran JSON reader ignores the key.
    #[serde(default)]
    pub supercell: Option<[i32; 3]>,
}

/// In-plane 2D Bravais basis used for a slab cross-section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum InPlaneBasis {
    /// Primitive cell of the true 2D Bravais lattice of the layer
    /// (translation stabilizer of the layer atom set — the finest
    /// in-plane cell, e.g. hexagonal for fcc (111)).
    Primitive,
    /// Right-angle in-plane cell (searched among small-integer
    /// combinations of the primitive vectors; yields a 90°/90°/90°
    /// 3D cell when available).
    Orthogonal,
    /// In-plane cell built from conventional-cell translations only
    /// (legacy behaviour; possibly a supercell of the true primitive).
    Conventional,
}

impl Default for InPlaneBasis {
    fn default() -> Self {
        InPlaneBasis::Primitive
    }
}

impl InPlaneBasis {
    /// UI-facing index (0 = primitive, 1 = orthogonal, 2 = conventional).
    pub const fn as_u8(self) -> u8 {
        match self {
            InPlaneBasis::Primitive => 0,
            InPlaneBasis::Orthogonal => 1,
            InPlaneBasis::Conventional => 2,
        }
    }

    pub const fn from_u8(v: u8) -> InPlaneBasis {
        match v {
            0 => InPlaneBasis::Primitive,
            1 => InPlaneBasis::Orthogonal,
            _ => InPlaneBasis::Conventional,
        }
    }
}

/// Provenance of the slab cross-section applied to this structure.
#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct SlabMeta {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    /// Slab start position along the surface normal, in Å.
    pub start_ang: f32,
    /// Slab thickness (normal component), in Å.
    pub thickness_ang: f32,
    /// In-plane supercell expansion: in-plane cell = U×V of the chosen
    /// 2D basis. `#[serde(default)]` — older JSON reads as 1×1.
    #[serde(default = "default_slab_uv")]
    pub u: u8,
    #[serde(default = "default_slab_uv")]
    pub v: u8,
    /// Which in-plane 2D basis was used.
    #[serde(default)]
    pub basis: InPlaneBasis,
    /// MS-style explicit in-plane vectors `(i j k)` (conventional basis)
    /// — when set, these (rather than `u`×`v` of `basis`) define the
    /// in-plane cell. `#[serde(default)]` — older JSON reads as absent.
    #[serde(default)]
    pub u_vec: Option<[i32; 3]>,
    #[serde(default)]
    pub v_vec: Option<[i32; 3]>,
}

fn default_slab_uv() -> u8 {
    1
}

/// Provenance of the vacuum layer applied to this structure.
#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct VacuumMeta {
    /// 1, 2, 3 → extended along the a / b / c axis.
    pub axis: u8,
    /// Vacuum thickness in Å.
    pub thickness_ang: f32,
    /// false = vacuum on top only, true = split V/2 on each side.
    pub both_sides: bool,
    /// Vertical position of the structure inside the cell, 0..=1:
    /// 0 = structure at the bottom (vacuum on top), 1 = at the top
    /// (vacuum below), 0.5 = centred. Absent in older JSON — derived
    /// from `both_sides` (0.0 / 0.5).
    #[serde(default)]
    pub position: f32,
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

    /// Recover cell parameters (a, b, c, alpha, beta, gamma) from arbitrary
    /// Cartesian lattice vectors (norms + mutual angles). Used after slab /
    /// vacuum construction, where the new basis vectors are arbitrary.
    pub fn from_cartesian_vectors(a: Vec3, b: Vec3, c: Vec3) -> Lattice {
        let a_len = a.length();
        let b_len = b.length();
        let c_len = c.length();
        let cos_a = (b_len * c_len > 0.0).then(|| b.dot(c) / (b_len * c_len)).unwrap_or(0.0).clamp(-1.0, 1.0);
        let cos_b = (a_len * c_len > 0.0).then(|| a.dot(c) / (a_len * c_len)).unwrap_or(0.0).clamp(-1.0, 1.0);
        let cos_g = (a_len * b_len > 0.0).then(|| a.dot(b) / (a_len * b_len)).unwrap_or(0.0).clamp(-1.0, 1.0);
        Lattice {
            a: a_len,
            b: b_len,
            c: c_len,
            alpha: cos_a.acos().to_degrees(),
            beta: cos_b.acos().to_degrees(),
            gamma: cos_g.acos().to_degrees(),
        }
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

#[derive(Debug, Clone, Deserialize, Serialize, Default)]
pub struct PhononModesData {
    pub mode_index: u32,
    pub frequency: f32,
    pub ir_intensity: f32,
    pub raman_activity: f32,
    pub mode_charge_norm: f32,
    #[serde(default)]
    pub atom_displacements: Vec<DisplacementData>,
}

#[derive(Debug, Clone, Deserialize, Serialize, Default)]
pub struct DisplacementData {
    pub dx: f32,
    pub dy: f32,
    pub dz: f32,
    pub contribution: f32,
}

impl CrystalData {
    pub fn from_json(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        let data: CrystalData = serde_json::from_str(&content)?;
        Ok(data)
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

    /// Atoms for display.
    ///
    /// A plain structure is shown as the full periodic cell via
    /// [`expand_to_cell`] (±1 images on every axis, so the equivalent
    /// atoms on the box edges and faces are all drawn). A slab/vacuum
    /// cut region or a confirmed supercell is NOT an infinite tiling:
    /// the display is the STORED atoms plus the equivalent face/edge/
    /// corner copies on every PERIODIC axis ([`display_boundary_complete`]).
    /// An atom sitting ON a periodic cell face (frac ≈ 0 or ≈ 1) is
    /// shared with the neighbouring cell, so a copy is drawn on the
    /// opposite face; interior atoms keep their single in-cell copy.
    /// The out-of-plane / vacuum axis of a cut region is never copied
    /// (no vertical slab stacking) unless a supercell expansion made it
    /// periodic.
    pub fn display_positions(&self) -> (Vec<Vec3>, Vec<usize>, Vec<Vec3>) {
        if self.slab.is_some() || self.vacuum.is_some() || self.supercell.is_some() {
            return self.display_boundary_complete();
        }
        self.expand_to_cell()
    }

    /// Which lattice axes are periodic for display purposes, per axis
    /// index 0=a, 1=b, 2=c.
    ///
    /// A plain structure (or a plain confirmed supercell) tiles on all
    /// three axes. A slab/vacuum cut region only tiles on its in-plane
    /// axes; its out-of-plane axis (the slab normal c, or the vacuum
    /// axis) stays non-periodic unless a supercell expansion covers it
    /// (`supercell[axis] > 1`).
    fn periodic_display_axes(&self) -> [bool; 3] {
        let mut p = [true, true, true];
        if self.slab.is_some() {
            // build_slab re-frames the cut so the surface normal is c.
            let c_periodic = self.supercell.map_or(false, |sc| sc[2] > 1);
            p[2] = c_periodic;
        } else if let Some(v) = &self.vacuum {
            let axis = ((v.axis as u32) - 1).min(2) as usize;
            let periodic = self.supercell.map_or(false, |sc| sc[axis] > 1);
            p[axis] = periodic;
        }
        p
    }

    /// The stored atoms, plus equivalent face/edge/corner copies on
    /// every PERIODIC axis ([`periodic_display_axes`]). An atom on a
    /// periodic face (frac ≈ 0 or ≈ 1) gets copies on BOTH faces of the
    /// drawn box (0 and 1, wrapped so no copy lands outside); interior
    /// atoms keep their single in-cell copy. This is the same ±1
    /// completion that [`expand_to_cell`] does for plain structures,
    /// but axis-restricted so a slab/vacuum cut region is not tiled out
    /// of plane.
    fn display_boundary_complete(&self) -> (Vec<Vec3>, Vec<usize>, Vec<Vec3>) {
        let vecs = self.lattice.to_vectors();
        let inv = self.lattice.inverse_vectors();
        let periodic = self.periodic_display_axes();
        let mut positions = Vec::new();
        let mut parents = Vec::new();
        let mut offsets = Vec::new();
        const FACE_EPS: f32 = 1e-3;
        for (a_i, atom) in self.atoms.iter().enumerate() {
            let frac = if self.positions_fractional {
                Vec3::new(atom.x, atom.y, atom.z)
            } else {
                let cart = Vec3::new(atom.x, atom.y, atom.z);
                Lattice::apply_inverse(&inv, cart)
            };
            // Candidate in-cell coordinate per axis: a face atom on a
            // periodic axis contributes both faces (0 and 1); every other
            // case contributes the single in-cell value.
            let cands: [([f32; 2], usize); 3] = [
                Self::face_candidates(frac.x, periodic[0], FACE_EPS),
                Self::face_candidates(frac.y, periodic[1], FACE_EPS),
                Self::face_candidates(frac.z, periodic[2], FACE_EPS),
            ];
            for xi in 0..cands[0].1 {
                for yi in 0..cands[1].1 {
                    for zi in 0..cands[2].1 {
                        let f = Vec3::new(cands[0].0[xi], cands[1].0[yi], cands[2].0[zi]);
                        let cart = f.x * vecs[0] + f.y * vecs[1] + f.z * vecs[2];
                        positions.push(cart);
                        parents.push(a_i);
                        offsets.push(f - frac);
                    }
                }
            }
        }
        (positions, parents, offsets)
    }

    /// Candidate in-cell coordinate values along one axis: a face atom
    /// (frac within FACE_EPS of 0 or 1) on a PERIODIC axis contributes
    /// both faces [0, 1]; all other cases contribute the single in-cell
    /// value [f, f] (count 1).
    fn face_candidates(f: f32, periodic: bool, eps: f32) -> ([f32; 2], usize) {
        if periodic && (f < eps || f > 1.0 - eps) {
            ([0.0, 1.0], 2)
        } else {
            ([f, f], 1)
        }
    }

    /// The two in-plane lattice-vector indices for a slab/vacuum cut
    /// region, or `None` for plain structures (full ±1 cell tiling) and
    /// for slab/vacuum structures whose in-plane cell was already
    /// supercell-expanded. A slab's surface normal is the c axis, so its
    /// in-plane axes are a,b; a vacuum-only structure is in-plane on the
    /// two axes OTHER than its vacuum axis. (Helper for the preview/
    /// camera logic; the display gate itself is
    /// [`display_positions`]/[`periodic_display_axes`].)
    pub fn display_inplane_axes(&self) -> Option<[usize; 2]> {
        let axes = if self.slab.is_some() {
            [0, 1]
        } else if let Some(v) = &self.vacuum {
            let axis = (v.axis - 1).min(2) as usize;
            let rest: Vec<usize> = (0..3).filter(|&k| k != axis).collect();
            [rest[0], rest[1]]
        } else {
            return None;
        };
        // A confirmed in-plane supercell expansion already spans those
        // in-plane cells — the 2×2 display block would multiply the span
        // again, so disable the block on any in-plane axis the supercell
        // expanded (the cell itself now shows the full 2D pattern).
        if let Some(sc) = &self.supercell {
            if sc[axes[0]] > 1 || sc[axes[1]] > 1 {
                return None;
            }
        }
        Some(axes)
    }

    /// Cell-boundary edges to draw. Plain structures and confirmed
    /// supercells: the full 3D box (12 edges — after a supercell apply
    /// this is the single MERGED box, internal edges gone). A plain slab
    /// cut (no vacuum yet): only the 4 in-plane (bottom-face) edges —
    /// the cut plane is not a 3D periodic cell, so the boundary is the
    /// in-plane rectangle whose shape is set by the U/V/hkl in-plane
    /// cell. Once a VACUUM layer has been built (on the slab or
    /// standalone), the structure is a definite 3D slab+vacuum box of
    /// height c = T + V — draw the complete 12-edge cell frame.
    pub fn display_box_edges(&self) -> Vec<(usize, usize)> {
        if self.vacuum.is_some() {
            Self::cell_edges().to_vec()
        } else if self.slab.is_some() {
            vec![(0, 1), (1, 2), (2, 3), (3, 0)]
        } else {
            Self::cell_edges().to_vec()
        }
    }
    /// Expand a single atom at fractional coords to cell images
    pub fn expand_single_atom(&self, frac: Vec3) -> (Vec<Vec3>, Vec<Vec3>) {
        let vecs = self.lattice.to_vectors();
        let mut positions = Vec::new();
        let mut offsets = Vec::new();
        for di in -1..=1_i32 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let off = Vec3::new(di as f32, dj as f32, dk as f32);
                    let f = frac + off;
                    if f.x >= -0.001 && f.x < 1.001
                        && f.y >= -0.001 && f.y < 1.001
                        && f.z >= -0.001 && f.z < 1.001
                    {
                        positions.push(f.x * vecs[0] + f.y * vecs[1] + f.z * vecs[2]);
                        offsets.push(off);
                    }
                }
            }
        }
        (positions, offsets)
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

#[cfg(test)]
mod tests {
    use super::*;
    fn slab_data() -> CrystalData {
        // (001) Cu slab, T = 2.5 Å, primitive L2 in-plane cell (a′ = a/√2)
        CrystalData {
            lattice: Lattice {
                a: 2.556,
                b: 2.556,
                c: 2.5,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
            atoms: vec![
                AtomData { element: "Cu".into(), x: 0.0, y: 0.0, z: 0.0, label: String::new() },
                AtomData { element: "Cu".into(), x: 1.278, y: 1.278, z: 1.8075, label: String::new() },
            ],
            positions_fractional: false,
            modified: true,
            phonon_modes: None,
            slab: Some(SlabMeta {
                h: 0,
                k: 0,
                l: 1,
                start_ang: 0.0,
                thickness_ang: 2.5,
                u: 1,
                v: 1,
                basis: InPlaneBasis::Primitive,
                u_vec: None,
                v_vec: None,
            }),
            vacuum: None,
            supercell: None,
        }
    }

    #[test]
    fn slab_structure_inplane_display_rule() {
        let d = slab_data();
        let (pos, parents, offsets) = d.display_positions();
        // boundary-sharing rule: the corner atom (0,0) is shared by 4
        // cells → all 4 corners of the 2×2 block; the interior atom
        // (½,½) is drawn once, inside the in-plane cell only
        assert_eq!(pos.len(), 5);
        assert_eq!(parents, vec![0,0,0,0,1]);
        let want: std::collections::BTreeSet<(i32, i32, i32)> =
            [(0,0,0), (1,0,0), (0,1,0), (1,1,0)].iter().cloned().collect();
        let got: std::collections::BTreeSet<(i32, i32, i32)> =
            (0..4).map(|k| (offsets[k].x as i32, offsets[k].y as i32, offsets[k].z as i32)).collect();
        assert_eq!(got, want);
        // the interior atom's single copy: zero offset, in-cell position
        assert_eq!(offsets[4], Vec3::ZERO);
        // no copy extends the c (vacuum/normal) axis
        assert!(offsets.iter().all(|o| o.z == 0.0));
        // the full ±1 tiling of this data keeps 8 + 1 = 9 in-cell images
        let (tiled, _, _) = d.expand_to_cell();
        assert_eq!(tiled.len(), 9);
        assert!(pos.len() < tiled.len());
    }

    #[test]
    fn slab_structure_shows_inplane_rectangle_not_box() {
        let mut d = slab_data();
        let edges = d.display_box_edges();
        // slab cut, no vacuum yet: only the 4 in-plane (bottom-face)
        // edges — no 3D box
        assert_eq!(edges, vec![(0,1),(1,2),(2,3),(3,0)]);
        // once a vacuum layer is built, the full 12-edge slab+vacuum
        // frame comes back
        d.vacuum = Some(VacuumMeta { axis: 3, thickness_ang: 15.0, both_sides: false, position: 0.0 });
        assert_eq!(d.display_box_edges().len(), 12);
        assert_eq!(d.display_positions().0.len(), 5, "display is unaffected by the box rule");
    }

    #[test]
    fn vacuum_only_axes_follow_vacuum_axis() {
        let mut d = slab_data();
        d.slab = None;
        d.vacuum = Some(VacuumMeta { axis: 1, thickness_ang: 15.0, both_sides: false, position: 0.0 });
        // vacuum on the a-axis → in-plane axes are b,c (no copies along a)
        let (pos, _, offsets) = d.display_positions();
        // corner atom → 4 shared (b,c) corners; interior atom → 1; all on b,c
        assert_eq!(pos.len(), 5);
        assert!(offsets.iter().all(|o| o.x == 0.0));
        assert_eq!(d.display_inplane_axes(), Some([1, 2]));
        // a vacuum layer means a definite 3D box → full 12-edge frame
        assert_eq!(d.display_box_edges().len(), 12);
    }

    #[test]
    fn plain_structure_still_tiles_the_cell() {
        let mut d = slab_data();
        d.slab = None; // ordinary asymmetric structure
        d.vacuum = None;
        let (plain, _, _) = d.display_positions();
        let (tiled, _, _) = d.expand_to_cell();
        assert_eq!(plain.len(), tiled.len());
        assert_eq!(d.display_box_edges().len(), 12);
    }

    #[test]
    fn confirmed_plain_supercell_renders_face_edge_equivalents() {
        // A confirmed 2×2×2 supercell of a plain structure: the stored
        // atoms are completed with face/edge/corner equivalents on all
        // three axes (the box is merged, 12 edges), matching the plain
        // ±1-tiling behaviour.  atom0 at (0,0,0) is on all three faces
        // → 2³ = 8 copies; atom1 at (½,½,0.723) is interior → 1 copy.
        let mut d = slab_data();
        d.slab = None;
        d.vacuum = None;
        d.supercell = Some([2, 2, 2]);
        let (pos, _, offsets) = d.display_positions();
        assert_eq!(pos.len(), 9, "8 copies of corner atom + 1 interior");
        // all offsets stay within {0,1}³ — no copy goes outside the box
        assert!(offsets.iter().all(|o| {
            o.x >= 0.0 && o.x <= 1.0 && o.y >= 0.0 && o.y <= 1.0 && o.z >= 0.0 && o.z <= 1.0
        }));
        assert_eq!(d.display_box_edges().len(), 12);
        // a 1×1×1 no-op still renders the same boundary-complete display
        d.supercell = Some([1, 1, 1]);
        assert_eq!(d.display_positions().0.len(), 9);
    }

    #[test]
    fn inplane_supercell_still_completes_inplane_faces() {
        // 2× in-plane supercell on a slab: the in-plane cell is enlarged
        // so the old 2×2 block is replaced by boundary completion on the
        // (now larger) in-plane cell.  atom0 at (0,0,0) → 4 in-plane
        // copies; atom1 interior → 1 copy.  c-axis stays non-periodic.
        let mut d = slab_data();
        d.supercell = Some([2, 1, 1]);
        assert_eq!(d.display_inplane_axes(), None);
        let (pos, _, offsets) = d.display_positions();
        assert_eq!(pos.len(), 5, "4 in-plane corner copies + 1 interior");
        // no offset has a z-component (c is not periodic)
        assert!(offsets.iter().all(|o| o.z == 0.0));
        assert_eq!(d.display_box_edges(), vec![(0, 1), (1, 2), (2, 3), (3, 0)]);

        // c-only expansion makes the c-axis periodic too: atom0 now gets
        // 2×2×2 = 8 copies, atom1 still 1 → 9 total.
        let mut c = slab_data();
        c.supercell = Some([1, 1, 2]);
        assert_eq!(c.display_inplane_axes(), Some([0, 1]));
        assert_eq!(c.display_positions().0.len(), 9);
    }

    #[test]
    fn slab_vacuum_supercell_adds_c_face_copies() {
        // Slab + vacuum + c-supercell: all three axes periodic.
        // atom0 (0,0,0) → 8 copies; atom1 (½,½, 0.723) → 1. Total 9.
        let mut d = slab_data();
        d.vacuum = Some(VacuumMeta { axis: 3, thickness_ang: 15.0, both_sides: false, position: 0.0 });
        d.supercell = Some([1, 1, 2]);
        let periodic = d.periodic_display_axes();
        assert_eq!(periodic, [true, true, true]);
        let (pos, _, _) = d.display_positions();
        assert_eq!(pos.len(), 9);
        assert_eq!(d.display_box_edges().len(), 12);
    }

    #[test]
    fn plain_supercell_display_completes_merged_box_boundary() {
        // Plain structure confirmed as a 2x2x2 supercell: the display is
        // the 16 tiled copies PLUS the face/edge/corner equivalents on the
        // merged-box boundary (the same completion the live preview uses).
        // Corner atom -> 8 box corners, face atoms -> 4 + 4 boundary copies,
        // interior atom -> 1: total 35 display positions.
        let mut d = slab_data();
        d.slab = None;
        d.vacuum = None;
        let previewed = crate::slab::build_supercell(&d, &crate::slab::SupercellParams { x: 2, y: 2, z: 2 }).unwrap();
        let (pos, _, _) = previewed.display_positions();
        assert!(pos.len() > previewed.atoms.len(), "boundary copies must exist beyond the 16 tiled atoms");
        assert_eq!(pos.len(), 35, "16 tiled + 19 face/edge/corner equivalents");
        // every copy lands on or inside the merged box (no wrap-out)
        let inv = previewed.lattice.inverse_vectors();
        for p in &pos {
            let f = Lattice::apply_inverse(&inv, *p);
            assert!((f.x >= -1e-3 && f.x <= 1.001)
                && (f.y >= -1e-3 && f.y <= 1.001)
                && (f.z >= -1e-3 && f.z <= 1.001));
        }
        // the far corner of the merged box carries an equivalent atom
        let [va, vb, vc] = previewed.lattice.to_vectors();
        let far = va + vb + vc; // far corner of the enlarged box
        let hit = pos.iter().any(|p| (p - far).length() < 1e-3);
        assert!(hit, "far corner equivalent missing");
    }

    #[test]
    fn slab_c_supercell_display_adds_top_face_copies() {
        // c-superlaced slab: the c axis becomes periodic, so the display
        // gains copies on the top face of the merged c-box, while the slab
        // keeps its in-plane corner sharing.
        let d = slab_data();
        let previewed = crate::slab::build_supercell(&d, &crate::slab::SupercellParams { x: 1, y: 1, z: 2 }).unwrap();
        assert_eq!(previewed.supercell, Some([1, 1, 2]));
        assert!(previewed.slab.is_some());
        let (pos, _, _) = previewed.display_positions();
        assert_eq!(pos.len(), 14, "4 tiled + 10 face/edge/corner copies");
        let inv = previewed.lattice.inverse_vectors();
        // c copies exist on the top face (z = 1 of the merged box)
        let z_fracs: Vec<f32> = pos.iter()
            .map(|p| Lattice::apply_inverse(&inv, *p).z)
            .collect();
        assert!(z_fracs.iter().any(|z| *z > 0.999), "top-face c copies expected");
        assert!(z_fracs.iter().all(|z| *z >= -1e-3 && *z <= 1.001), "no copy may escape the box");
        // slab-only (no vacuum yet) keeps the 4 in-plane rectangle frame
        assert_eq!(previewed.display_box_edges().len(), 4);
    }
}
