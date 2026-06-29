//! Picking via 3D→2D projection: project atoms to screen, pick nearest to cursor

use bevy::prelude::*;
use bevy::window::PrimaryWindow;
use bevy_egui::EguiContexts;

#[derive(Resource)]
pub struct PickingState {
    pub selected: i32,
    pub hovered: i32,
    pub atom_positions: Vec<Vec3>,
    pub atom_material_handles: Vec<Handle<StandardMaterial>>,
    pub atom_entities: Vec<Entity>,
    pub parent_indices: Vec<usize>,
    pub click_start: Option<Vec2>,
    pub modified: bool,
}

impl PickingState {
    pub fn new(
        positions: Vec<Vec3>,
        handles: Vec<Handle<StandardMaterial>>,
        entities: Vec<Entity>,
        parent_indices: Vec<usize>,
    ) -> Self {
        Self {
            selected: -1, hovered: -1,
            atom_positions: positions,
            atom_material_handles: handles,
            atom_entities: entities,
            parent_indices,
            click_start: None,
            modified: false,
        }
    }

    /// Add new atom images (one or more from cell expansion) for a new parent atom
    pub fn add_images(
        &mut self,
        positions: Vec<Vec3>,
        handles: Vec<Handle<StandardMaterial>>,
        entities: Vec<Entity>,
        parent_idx: usize,
    ) {
        let n = positions.len();
        for entity in entities {
            self.atom_entities.push(entity);
        }
        for pos in positions {
            self.atom_positions.push(pos);
        }
        for handle in handles {
            self.atom_material_handles.push(handle);
        }
        for _ in 0..n {
            self.parent_indices.push(parent_idx);
        }
    }

    /// Remove all images belonging to the given parent atom
    pub fn remove_images(&mut self, parent_idx: usize) {
        let mut indices: Vec<usize> = (0..self.parent_indices.len())
            .filter(|&i| self.parent_indices[i] == parent_idx)
            .collect();
        // Remove in descending order to avoid index shifting
        indices.sort_unstable_by(|a, b| b.cmp(a));
        for i in indices {
            self.atom_entities.remove(i);
            self.atom_positions.remove(i);
            self.atom_material_handles.remove(i);
            self.parent_indices.remove(i);
        }
        // Renumber parent indices for atoms after the removed one
        for p in self.parent_indices.iter_mut() {
            if *p > parent_idx {
                *p -= 1;
            }
        }
    }
}

fn build_mvp(cam_gt: &GlobalTransform, proj: &Projection) -> Mat4 {
    let view_mat = Mat4::from(cam_gt.affine().inverse());
    let proj_mat = match proj {
        Projection::Perspective(p) => perspective_mat(p.fov, p.aspect_ratio, p.near, p.far),
        Projection::Orthographic(o) => orthographic_mat(o),
    };
    proj_mat * view_mat
}

fn perspective_mat(fov: f32, aspect: f32, near: f32, far: f32) -> Mat4 {
    let tan_half = (fov * 0.5).tan();
    Mat4::from_cols(
        Vec4::new(1.0 / (aspect * tan_half), 0.0, 0.0, 0.0),
        Vec4::new(0.0, 1.0 / tan_half, 0.0, 0.0),
        Vec4::new(0.0, 0.0, -(far + near) / (far - near), -1.0),
        Vec4::new(0.0, 0.0, -(2.0 * far * near) / (far - near), 0.0),
    )
}

fn orthographic_mat(o: &OrthographicProjection) -> Mat4 {
    let hw = o.area.width() * 0.5;
    let hh = o.area.height() * 0.5;
    Mat4::from_cols(
        Vec4::new(1.0 / hw, 0.0, 0.0, 0.0),
        Vec4::new(0.0, 1.0 / hh, 0.0, 0.0),
        Vec4::new(0.0, 0.0, -2.0 / (o.far - o.near), 0.0),
        Vec4::new(0.0, 0.0, -(o.far + o.near) / (o.far - o.near), 1.0),
    )
}

fn nearest_atom(mvp: &Mat4, positions: &[Vec3], cursor: Vec2, vp_size: Vec2, max_dist: f32) -> i32 {
    let mut best: i32 = -1;
    let mut best_d2 = max_dist.powi(2);
    for (i, pos) in positions.iter().enumerate() {
        let clip = *mvp * pos.extend(1.0);
        if clip.w <= 0.0 { continue; }
        let sx = (clip.x / clip.w * 0.5 + 0.5) * vp_size.x;
        let sy = (1.0 - (clip.y / clip.w * 0.5 + 0.5)) * vp_size.y;
        let d2 = (sx - cursor.x).powi(2) + (sy - cursor.y).powi(2);
        if d2 < best_d2 { best_d2 = d2; best = i as i32; }
    }
    best
}

pub fn click_pick(
    mut picking: ResMut<PickingState>,
    mouse_btn: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window, With<PrimaryWindow>>,
    camera_q: Query<(&Camera, &GlobalTransform, &Projection), With<crate::MainCamera>>,
    mut contexts: EguiContexts,
) {
    // Ignore 3D picks when cursor is over an egui panel
    if contexts.ctx_mut().is_pointer_over_area() {
        if mouse_btn.just_pressed(MouseButton::Left) || mouse_btn.just_released(MouseButton::Left) {
            picking.click_start = None;
        }
        return;
    }
    let Ok(window) = windows.get_single() else { return };
    let Ok((camera, cam_gt, proj)) = camera_q.get_single() else { return };

    if mouse_btn.just_pressed(MouseButton::Left) {
        picking.click_start = window.cursor_position();
    }
    if mouse_btn.just_released(MouseButton::Left) {
        if let (Some(start), Some(end)) = (picking.click_start, window.cursor_position()) {
            if start.distance(end) < 5.0 {
                let vp = camera.logical_viewport_size().unwrap_or(
                    Vec2::new(window.width() as f32, window.height() as f32)
                );
                let mvp = build_mvp(cam_gt, proj);
                let click_dist = vp.y * 0.06;
                picking.selected = nearest_atom(&mvp, &picking.atom_positions, end, vp, click_dist);
                if picking.selected >= 0 {
                    let pos = picking.atom_positions[picking.selected as usize];
                    println!("Selected atom {}: ({:.4}, {:.4}, {:.4})", picking.selected, pos.x, pos.y, pos.z);
                }
            }
        }
        picking.click_start = None;
    }
}

pub fn hover_pick(
    mut picking: ResMut<PickingState>,
    windows: Query<&Window, With<PrimaryWindow>>,
    camera_q: Query<(&Camera, &GlobalTransform, &Projection), With<crate::MainCamera>>,
    mut contexts: EguiContexts,
) {
    // Clear hover when cursor is over an egui panel
    if contexts.ctx_mut().is_pointer_over_area() {
        picking.hovered = -1;
        return;
    }
    if picking.click_start.is_some() { return; }
    let Ok(window) = windows.get_single() else { return };
    let Ok((camera, cam_gt, proj)) = camera_q.get_single() else { return };
    let cursor = match window.cursor_position() {
        Some(c) => c, None => { picking.hovered = -1; return; }
    };
    let vp = camera.logical_viewport_size().unwrap_or(
        Vec2::new(window.width() as f32, window.height() as f32)
    );
    let mvp = build_mvp(cam_gt, proj);
    let hover_dist = vp.y * 0.06;
    picking.hovered = nearest_atom(&mvp, &picking.atom_positions, cursor, vp, hover_dist);
}

pub fn highlight_atoms(
    picking: Res<PickingState>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let sel_parent = if picking.selected >= 0 && (picking.selected as usize) < picking.parent_indices.len() {
        Some(picking.parent_indices[picking.selected as usize])
    } else { None };
    let hov_parent = if picking.hovered >= 0 && (picking.hovered as usize) < picking.parent_indices.len() {
        Some(picking.parent_indices[picking.hovered as usize])
    } else { None };

    for i in 0..picking.atom_material_handles.len() {
        let handle = &picking.atom_material_handles[i];
        let Some(mat) = materials.get_mut(handle) else { continue };
        let my_parent = if i < picking.parent_indices.len() {
            Some(picking.parent_indices[i])
        } else { None };

        if sel_parent.is_some() && my_parent == sel_parent {
            mat.emissive = LinearRgba::rgb(1.0, 0.5, 0.1);
        } else if hov_parent.is_some() && my_parent == hov_parent {
            mat.emissive = LinearRgba::rgb(0.3, 0.3, 0.3);
        } else {
            mat.emissive = LinearRgba::BLACK;
        }
    }
}
