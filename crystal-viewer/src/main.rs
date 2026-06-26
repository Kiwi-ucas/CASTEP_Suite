//! Crystal Viewer — Interactive 3D crystal structure visualization

mod crystal; mod resources; mod picking; mod ui;

use bevy::prelude::*;
use bevy::input::mouse::{MouseMotion, MouseScrollUnit, MouseWheel};
use bevy::render::mesh::{Mesh, Indices, PrimitiveTopology};
use crystal::{CrystalData, Lattice};
use picking::{PickingState, click_pick, hover_pick, highlight_atoms};
use ui::{AtomInfo, CrystalMeta, ui_system};
use std::f32::consts::PI;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let json_path = args.get(1).map(|s| s.as_str()).unwrap_or("");

    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(bevy_egui::EguiPlugin)
        .insert_resource(CrystalPath(json_path.to_string()))
        .add_systems(Startup, setup)
        .add_systems(Update, orbit_camera)
        .add_systems(Update, (click_pick, hover_pick).chain())
        .add_systems(Update, highlight_atoms)
        .add_systems(Update, move_selected_atom.after(highlight_atoms))
        .add_systems(Update, (add_atom_system, delete_atom_system))
        .add_systems(Update, ui_system.after(highlight_atoms))
        .add_systems(Update, (display_mode_system, sync_atom_radii))
        .add_systems(Update, toggle_projection)
        .run();
}

#[derive(Resource)] struct CrystalPath(String);
#[derive(Resource)] struct CrystalScene { center: Vec3 }
#[derive(Resource, Clone)] struct CrystalStore { data: CrystalData, json_path: String }
#[derive(Resource)] struct LatticeData { vecs: [Vec3; 3], inv: [Vec3; 3] }
#[derive(Resource)] struct ImageOffsets(pub Vec<Vec3>);  // fractional offset for each expanded atom

/// State machine for Add Atom UI flow
#[derive(Resource, Default)]
struct AddAtomState {
    show_table: bool,
    selected_element: Option<String>,
    coord_x: String,
    coord_y: String,
    coord_z: String,
    /// (element, frac_x, frac_y, frac_z) — populated by UI, consumed by add_atom_system
    pending: Option<(String, f32, f32, f32)>,
    next_id: u32,
}

#[derive(Resource)]
struct CachedSphere(Handle<Mesh>);
#[derive(Resource)] struct MoveState { step: f32 }
#[derive(Component)] struct MainCamera;
#[derive(Component)] struct FollowCamera;
#[derive(Component)] struct AtomMarker;
#[derive(Component)] struct BondMarker;
#[derive(Component)] struct CellMarker;

#[derive(Resource)]
struct DisplayMode {
    mode: u8,          // 1=ball-stick, 2=space-filling, 3=wireframe
    show_bonds: bool,
    show_cell: bool,
}

/// Saved initial camera state for R-key reset
#[derive(Resource)]
struct CameraInit(CameraState);

#[derive(Resource, Clone)]
pub struct CameraState { pub focus: Vec3, pub radius: f32, pub yaw: f32, pub pitch: f32 }

#[derive(Default)] struct InputState { rotating: bool }

#[derive(Resource, PartialEq)]
enum ProjMode { Perspective, Orthographic }

fn ortho_projection(scale: f32) -> Projection {
    Projection::Orthographic(OrthographicProjection {
        scaling_mode: bevy::render::camera::ScalingMode::AutoMin {
            min_width: scale, min_height: scale,
        },
        near: -1000.0, far: 1000.0,
        ..OrthographicProjection::default_3d()
    })
}

fn toggle_projection(
    keys: Res<ButtonInput<KeyCode>>,
    mut proj_mode: ResMut<ProjMode>,
    mut ortho_scale: ResMut<OrthoScale>,
    mut cam_state: ResMut<CameraState>,
    mut camera_q: Query<&mut Projection, With<MainCamera>>,
) {
    if keys.just_pressed(KeyCode::KeyP) {
        *proj_mode = match *proj_mode {
            ProjMode::Perspective => {
                // Sync: derive ortho_scale from perspective radius
                ortho_scale.0 = cam_state.radius * 1.09;
                ProjMode::Orthographic
            }
            ProjMode::Orthographic => {
                // Sync: derive perspective radius from ortho_scale
                cam_state.radius = ortho_scale.0 / 1.09;
                ProjMode::Perspective
            }
        };
        let Ok(mut proj) = camera_q.get_single_mut() else { return };
        *proj = match *proj_mode {
            ProjMode::Perspective => Projection::Perspective(PerspectiveProjection {
                fov: 1.0, ..default()
            }),
            ProjMode::Orthographic => ortho_projection(ortho_scale.0),
        };
    }
}

#[derive(Resource)]
struct OrthoScale(f32);  // vertical world units visible in ortho mode

fn orbit_camera(
    mut camera_q: Query<&mut Transform, (With<MainCamera>, Without<FollowCamera>)>,
    mut proj_q: Query<&mut Projection, With<MainCamera>>,
    mut light_q: Query<&mut Transform, With<FollowCamera>>,
    mut cam_state: ResMut<CameraState>,
    cam_init: Res<CameraInit>,
    proj_mode: Res<ProjMode>,
    mut ortho_scale: ResMut<OrthoScale>,
    scene: Option<Res<CrystalScene>>,
    mut input: Local<InputState>,
    mouse_btn: Res<ButtonInput<MouseButton>>,
    mut mouse_motion: EventReader<MouseMotion>,
    mut mouse_wheel: EventReader<MouseWheel>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    let Ok(mut cam) = camera_q.get_single_mut() else { return };

    // Right mouse → rotate (no limit on pitch for infinite rotation)
    if mouse_btn.just_pressed(MouseButton::Right) { input.rotating = true; }
    if mouse_btn.just_released(MouseButton::Right) { input.rotating = false; }

    for motion in mouse_motion.read() {
        let d = motion.delta;
        if input.rotating {
            cam_state.yaw -= d.x * 0.005;
            cam_state.pitch += d.y * 0.005;   // no clamp — infinite rotation
        }
    }

    for ev in mouse_wheel.read() {
        let dy = match ev.unit {
            MouseScrollUnit::Line => ev.y * 0.1,
            MouseScrollUnit::Pixel => ev.y * 0.001,
        };
        if *proj_mode == ProjMode::Orthographic {
            ortho_scale.0 = (ortho_scale.0 - dy * ortho_scale.0 * 0.1).clamp(1.0, 200.0);
            cam_state.radius = ortho_scale.0 / 1.09;  // keep in sync
            // Rebuild ortho projection with new scale
            if let Ok(mut proj) = proj_q.get_single_mut() {
                *proj = ortho_projection(ortho_scale.0);
            }
        } else {
            cam_state.radius = (cam_state.radius - dy * cam_state.radius * 0.1).clamp(1.0, 100.0);
            ortho_scale.0 = cam_state.radius * 1.09;  // keep in sync
        }
    }

    if keys.just_pressed(KeyCode::KeyF) {
        if let Some(s) = scene.as_ref() { cam_state.focus = s.center; cam_state.radius = 10.0; }
    }
    if keys.just_pressed(KeyCode::KeyR) {
        *cam_state = cam_init.0.clone();
        ortho_scale.0 = cam_init.0.radius * 1.09;
    }

    let pos = cam_state.focus + spherical(cam_state.radius, cam_state.yaw, cam_state.pitch);
    cam.translation = pos;
    // Build rotation quaternion directly from yaw/pitch orbit angles.
    // spherical() maps (yaw,pitch) to offset; camera looks from offset toward focus.
    // q_pitch rotates -Z (default look) to the pitched direction;
    // q_yaw then yaws that result around world Y.
    let q_pitch = Quat::from_rotation_x(-cam_state.pitch);
    let q_yaw = Quat::from_rotation_y(cam_state.yaw);
    cam.rotation = q_yaw * q_pitch;

    for mut lt in light_q.iter_mut() {
        let cam_right = cam.rotation * Vec3::X;
        let cam_up = cam.rotation * Vec3::Y;
        let lp = pos + cam_right * cam_state.radius * 0.8 + cam_up * cam_state.radius * 0.6;
        lt.translation = lp;
        lt.look_at(cam_state.focus, Vec3::Y);
    }
}

fn spherical(r: f32, yaw: f32, pitch: f32) -> Vec3 {
    Vec3::new(r * pitch.cos() * yaw.sin(), r * pitch.sin(), r * pitch.cos() * yaw.cos())
}

// ── Atom movement ──

fn move_selected_atom(
    keys: Res<ButtonInput<KeyCode>>,
    mut picking: ResMut<PickingState>,
    mut atoms: Query<&mut Transform, With<AtomMarker>>,
    mut move_state: ResMut<MoveState>,
    mut crystal: ResMut<CrystalStore>,
    lattice: Res<LatticeData>,
    offsets: Res<ImageOffsets>,
) {
    if picking.selected < 0 { return; }
    let i = picking.selected as usize;
    if i >= picking.parent_indices.len() { return; }

    let step = move_state.step;
    let mut dx = 0.0_f32;
    let mut dy = 0.0_f32;
    let mut dz = 0.0_f32;

    if keys.just_pressed(KeyCode::KeyL) { dx = step; }
    if keys.just_pressed(KeyCode::KeyJ) { dx = -step; }
    if keys.just_pressed(KeyCode::KeyU) { dy = step; }
    if keys.just_pressed(KeyCode::KeyO) { dy = -step; }
    if keys.just_pressed(KeyCode::KeyI) { dz = step; }
    if keys.just_pressed(KeyCode::KeyK) { dz = -step; }

    if keys.just_pressed(KeyCode::BracketRight) {
        move_state.step = match move_state.step {
            0.01 => 0.05, 0.05 => 0.1, 0.1 => 0.5, 0.5 => 1.0, _ => 0.01,
        };
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        move_state.step = match move_state.step {
            1.0 => 0.5, 0.5 => 0.1, 0.1 => 0.05, 0.05 => 0.01, _ => 1.0,
        };
    }

    if dx != 0.0 || dy != 0.0 || dz != 0.0 {
        let parent = picking.parent_indices[i];
        let cart_delta = Vec3::new(dx, dy, dz);
        let frac_delta = Lattice::apply_inverse(&lattice.inv, cart_delta);

        if parent < crystal.data.atoms.len() {
            // Get current parent coords in fractional space
            let (mut pf_x, mut pf_y, mut pf_z) = if crystal.data.positions_fractional {
                (crystal.data.atoms[parent].x, crystal.data.atoms[parent].y, crystal.data.atoms[parent].z)
            } else {
                let cart = Vec3::new(crystal.data.atoms[parent].x, crystal.data.atoms[parent].y, crystal.data.atoms[parent].z);
                let frac = Lattice::apply_inverse(&lattice.inv, cart);
                (frac.x, frac.y, frac.z)
            };
            // Apply fractional delta
            pf_x += frac_delta.x;
            pf_y += frac_delta.y;
            pf_z += frac_delta.z;
            // Write back in original storage format
            if crystal.data.positions_fractional {
                crystal.data.atoms[parent].x = pf_x;
                crystal.data.atoms[parent].y = pf_y;
                crystal.data.atoms[parent].z = pf_z;
            } else {
                crystal.data.atoms[parent].x = lattice.vecs[0].x * pf_x + lattice.vecs[1].x * pf_y + lattice.vecs[2].x * pf_z;
                crystal.data.atoms[parent].y = lattice.vecs[0].y * pf_x + lattice.vecs[1].y * pf_y + lattice.vecs[2].y * pf_z;
                crystal.data.atoms[parent].z = lattice.vecs[0].z * pf_x + lattice.vecs[1].z * pf_y + lattice.vecs[2].z * pf_z;
            }
            crystal.data.modified = true;

            // Regenerate all images in fractional space, convert to Cartesian
            let parent_frac = Vec3::new(pf_x, pf_y, pf_z);
            let vecs = &lattice.vecs;
            let siblings: Vec<usize> = picking.parent_indices.iter().enumerate()
                .filter(|(_, &p)| p == parent)
                .map(|(j, _)| j)
                .collect();
            for &j in &siblings {
                let img_frac = parent_frac + offsets.0[j];
                let cart = img_frac.x * vecs[0] + img_frac.y * vecs[1] + img_frac.z * vecs[2];
                picking.atom_positions[j] = cart;
                let entity = picking.atom_entities[j];
                if let Ok(mut transform) = atoms.get_mut(entity) {
                    transform.translation = cart;
                }
            }
        }

        picking.modified = true;

        // Auto-save modified positions to JSON
        if let Err(e) = crystal.data.write_to_file(&crystal.json_path) {
            eprintln!("  Failed to save modified positions: {}", e);
        }
    }
}

fn add_atom_system(
    mut add_state: ResMut<AddAtomState>,
    mut crystal: ResMut<CrystalStore>,
    mut picking: ResMut<PickingState>,
    mut offsets: ResMut<ImageOffsets>,
    mut atom_info: ResMut<ui::AtomInfo>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    sphere: Res<CachedSphere>,
    mut commands: Commands,
) {
    let Some((el, fx, fy, fz)) = add_state.pending.take() else { return };

    let label = format!("new_{}", add_state.next_id);
    add_state.next_id += 1;

    // Append to CrystalData
    crystal.data.atoms.push(crystal::AtomData {
        element: el.clone(), x: fx, y: fy, z: fz, label,
    });
    crystal.data.modified = true;

    let parent = crystal.data.atoms.len() - 1;
    let frac = Vec3::new(fx, fy, fz);
    let (positions, image_offsets) = crystal.data.expand_single_atom(frac);

    // Spawn entities for each image
    let color = resources::element_color(&el);
    let mat_handle = materials.add(StandardMaterial {
        base_color: color, metallic: 0.2, perceptual_roughness: 0.5, ..default()
    });
    let mut handles = Vec::new();
    let mut entities = Vec::new();
    for pos in &positions {
        let entity = commands.spawn((
            Mesh3d(sphere.0.clone()),
            MeshMaterial3d(mat_handle.clone()),
            Transform::from_translation(*pos),
            AtomMarker,
        )).id();
        handles.push(mat_handle.clone());
        entities.push(entity);
    }

    // Update PickingState
    picking.add_images(positions, handles, entities, parent);

    // Update ImageOffsets
    for off in image_offsets {
        offsets.0.push(off);
    }

    // Update AtomInfo
    atom_info.elements.push(el.clone());
    atom_info.radii.push(resources::covalent_radius(&el));

    // Auto-save
    crystal.data.modified = true;
    let _ = crystal.data.write_to_file(&crystal.json_path);
}

fn sync_atom_radii(
    atom_info: Res<ui::AtomInfo>,
    picking: Res<PickingState>,
    mut atoms: Query<&mut Transform, (With<AtomMarker>, Without<BondMarker>)>,
    display: Res<DisplayMode>,
    mut initialized: Local<bool>,
) {
    // Skip first-frame change detection (resource insertion triggers is_changed)
    if !*initialized { *initialized = true; return; }
    if display.mode != 1 || !atom_info.is_changed() { return; }
    for (i, mut t) in atoms.iter_mut().enumerate() {
        if i < picking.parent_indices.len() {
            let p = picking.parent_indices[i];
            if p < atom_info.elements.len() && p < atom_info.radii.len() {
                let default_r = resources::covalent_radius(&atom_info.elements[p]);
                t.scale = Vec3::splat(atom_info.radii[p] / default_r);
            }
        }
    }
}

fn delete_atom_system(
    keys: Res<ButtonInput<KeyCode>>,
    mut picking: ResMut<PickingState>,
    mut crystal: ResMut<CrystalStore>,
    mut offsets: ResMut<ImageOffsets>,
    mut atom_info: ResMut<ui::AtomInfo>,
    mut commands: Commands,
) {
    if !keys.just_pressed(KeyCode::KeyD) || picking.selected < 0 { return; }
    let i = picking.selected as usize;
    if i >= picking.parent_indices.len() { return; }
    let parent = picking.parent_indices[i];
    if parent >= crystal.data.atoms.len() { return; }

    // Collect entities and image-offset indices before mutation
    let to_despawn: Vec<Entity> = (0..picking.parent_indices.len())
        .filter(|&j| picking.parent_indices[j] == parent)
        .map(|j| picking.atom_entities[j])
        .collect();
    let off_indices: Vec<usize> = (0..picking.parent_indices.len())
        .filter(|&j| picking.parent_indices[j] == parent)
        .collect();

    // Despawn entities
    for entity in &to_despawn {
        commands.entity(*entity).despawn();
    }

    // Remove from ImageOffsets (descending)
    let mut sorted = off_indices.clone();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    for j in &sorted {
        offsets.0.remove(*j);
    }

    // Remove from PickingState (handles renumbering)
    picking.remove_images(parent);

    // Remove from canonical data
    crystal.data.atoms.remove(parent);
    if parent < atom_info.elements.len() {
        atom_info.elements.remove(parent);
        atom_info.labels.remove(parent);
        atom_info.radii.remove(parent);
    }

    picking.selected = -1;
    picking.modified = true;
    crystal.data.modified = true;
    let _ = crystal.data.write_to_file(&crystal.json_path);
}

// ── Geometry ──

fn uv_sphere(radius: f32, sec: u32, stk: u32) -> Mesh {
    let mut pos = Vec::new(); let mut nrm = Vec::new(); let mut idx = Vec::new();
    for i in 0..=stk {
        let phi = PI * i as f32 / stk as f32;
        let y = -radius * phi.cos(); let rr = radius * phi.sin();
        for j in 0..=sec {
            let t = 2.0 * PI * j as f32 / sec as f32;
            let x = rr * t.cos(); let z = rr * t.sin();
            pos.push([x, y, z]);
            let len = (x*x + y*y + z*z).sqrt().max(0.0001);
            nrm.push([x/len, y/len, z/len]);
        }
    }
    for i in 0..stk { for j in 0..sec {
        let f = i * (sec + 1) + j; let s = f + sec + 1;
        idx.extend_from_slice(&[f, s, f+1, s, s+1, f+1]);
    }}
    let mut m = Mesh::new(PrimitiveTopology::TriangleList, Default::default());
    m.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    m.insert_attribute(Mesh::ATTRIBUTE_NORMAL, nrm);
    m.insert_indices(Indices::U32(idx));
    m
}

// ── Scene ──

fn setup(
    mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>, crystal_path: Res<CrystalPath>,
) {
    let data = if crystal_path.0.is_empty() { default_cu_fcc() } else {
        CrystalData::from_json(&crystal_path.0).unwrap_or_else(|e| {
            eprintln!("Failed to load {}: {}. Using Cu FCC.", crystal_path.0, e);
            default_cu_fcc()
        })
    };

    let (positions, parent_indices, image_offsets) = data.expand_to_cell();
    let center = data.center();
    let corners = data.cell_corners();
    let edges = CrystalData::cell_edges();
    let n = positions.len();

    commands.insert_resource(LatticeData {
        vecs: data.lattice.to_vectors(),
        inv: data.lattice.inverse_vectors(),
    });

    let sphere = meshes.add(uv_sphere(0.5, 32, 32));
    commands.insert_resource(CachedSphere(sphere.clone()));
    let mut handles = Vec::with_capacity(n);
    let mut entities = Vec::with_capacity(n);

    for (i, pos) in positions.iter().enumerate() {
        let parent = parent_indices[i];
        let el = &data.atoms[parent].element;
        let color = resources::element_color(el);
        let mat_handle = materials.add(StandardMaterial {
            base_color: color, metallic: 0.2, perceptual_roughness: 0.5, ..default()
        });
        let entity = commands.spawn((
            Mesh3d(sphere.clone()), MeshMaterial3d(mat_handle.clone()),
            Transform::from_translation(*pos), AtomMarker,
        )).id();
        handles.push(mat_handle);
        entities.push(entity);
    }

    // Bonds
    let bond_mat = materials.add(StandardMaterial {
        base_color: Color::srgba(0.65, 0.65, 0.65, 0.5),
        alpha_mode: AlphaMode::Blend, ..default()
    });
    for i in 0..n { for j in (i+1)..n {
        // Skip bonds between images of the same parent atom
        if parent_indices[i] == parent_indices[j] { continue; }
        let pi = parent_indices[i];
        let pj = parent_indices[j];
        if resources::has_bond(&data.atoms[pi].element, &data.atoms[pj].element,
                                positions[i].distance(positions[j]), 1.2) {
            spawn_bond(&mut commands, &mut meshes, &bond_mat, positions[i], positions[j], BondMarker);
        }
    }}

    // Cell edges
    let edge_mat = materials.add(StandardMaterial {
        base_color: Color::srgba(0.4, 0.4, 0.4, 0.35),
        alpha_mode: AlphaMode::Blend, ..default()
    });
    for (a, b) in edges {
        let s = corners[a]; let e = corners[b];
        let mid = (s + e) * 0.5; let dir = (e - s).normalize_or_zero(); let len = s.distance(e);
        commands.spawn((
            Mesh3d(meshes.add(Cuboid::new(0.04, 0.04, len))),
            MeshMaterial3d(edge_mat.clone()),
            Transform::from_translation(mid).with_rotation(Quat::from_rotation_arc(Vec3::Z, dir)),
            CellMarker,
        ));
    }

    commands.insert_resource(PickingState::new(positions, handles, entities, parent_indices));
    commands.insert_resource(ImageOffsets(image_offsets));

    // Store crystal data + path for auto-save on modification
    commands.insert_resource(CrystalStore {
        data: data.clone(),
        json_path: crystal_path.0.clone(),
    });
    commands.insert_resource(MoveState { step: 0.1 });
    commands.insert_resource(AddAtomState::default());

    // Atom metadata for UI
    let elements: Vec<String> = data.atoms.iter().map(|a| a.element.clone()).collect();
    let labels: Vec<String> = data.atoms.iter().enumerate()
        .map(|(_i, a)| if a.label.is_empty() { a.element.clone() } else { a.label.clone() })
        .collect();
    let radii: Vec<f32> = elements.iter().map(|el| resources::covalent_radius(el)).collect();
    commands.insert_resource(AtomInfo::new(elements, labels, radii));
    commands.insert_resource(DisplayMode { mode: 1, show_bonds: false, show_cell: true });

    // Crystal metadata for UI (strip .json extension from filename)
    let fname = if crystal_path.0.is_empty() {
        "Cu FCC (demo)".to_string()
    } else {
        let raw = std::path::Path::new(&crystal_path.0)
            .file_name().map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| crystal_path.0.clone());
        raw.strip_suffix(".json").map(|s| s.to_string()).unwrap_or(raw)
    };
    commands.insert_resource(CrystalMeta {
        filename: fname,
        a: data.lattice.a, b: data.lattice.b, c: data.lattice.c,
        alpha: data.lattice.alpha, beta: data.lattice.beta, gamma: data.lattice.gamma,
    });

    // Compute default ortho scale from cell bounding box
    let cell_corners = data.cell_corners();
    let mut cmin = Vec3::splat(f32::MAX);
    let mut cmax = Vec3::splat(f32::MIN);
    for c in &cell_corners {
        cmin = cmin.min(*c);
        cmax = cmax.max(*c);
    }
    let cell_diag = (cmax - cmin).length();
    let ortho_scale = (cell_diag * 1.3).max(5.0);

    // Init camera (orthographic by default, looking down Z at XY plane)
    let radius = ortho_scale; let yaw = 0.0; let pitch = 0.0;
    commands.insert_resource(CameraState { focus: center, radius, yaw, pitch });
    commands.insert_resource(CameraInit(CameraState { focus: center, radius, yaw, pitch }));
    commands.insert_resource(CrystalScene { center });
    commands.insert_resource(ProjMode::Orthographic);
    commands.insert_resource(OrthoScale(ortho_scale));

    commands.spawn((
        DirectionalLight { illuminance: 8000.0, shadows_enabled: false, ..default() },
        Transform::default(), FollowCamera,
    ));
    commands.spawn((
        DirectionalLight { illuminance: 2000.0, ..default() },
        Transform::from_xyz(0.0, -5.0, 0.0).looking_at(center, Vec3::Y),
    ));

    commands.spawn((
        Camera3d::default(),
        ortho_projection(ortho_scale),
        Transform::default(),
        MainCamera,
    ));

    println!("Loaded {} atoms. Right-drag to rotate, scroll to zoom, click to select.", n);
}

fn spawn_bond(
    commands: &mut Commands, meshes: &mut ResMut<Assets<Mesh>>,
    material: &Handle<StandardMaterial>, a: Vec3, b: Vec3, marker: BondMarker,
) {
    let mid = (a + b) * 0.5; let dir = (b - a).normalize_or_zero(); let len = a.distance(b);
    commands.spawn((
        Mesh3d(meshes.add(Cylinder::new(0.08, len))),
        MeshMaterial3d(material.clone()),
        Transform::from_translation(mid).with_rotation(Quat::from_rotation_arc(Vec3::Y, dir)),
        marker,
        Visibility::Hidden,  // hidden by default, press B to show
    ));
}

fn display_mode_system(
    keys: Res<ButtonInput<KeyCode>>,
    mut display: ResMut<DisplayMode>,
    mut atoms: Query<&mut Transform, (With<AtomMarker>, Without<BondMarker>)>,
    picking: Res<PickingState>,
    atom_info: Res<ui::AtomInfo>,
    bonds: Query<Entity, With<BondMarker>>,
    cells: Query<Entity, With<CellMarker>>,
    mut commands: Commands,
) {
    let update_scales = |mode: &mut u8| -> bool {
        let changed = matches!(
            (keys.just_pressed(KeyCode::Digit1), keys.just_pressed(KeyCode::Digit2), keys.just_pressed(KeyCode::Digit3)),
            (true, _, _) | (_, true, _) | (_, _, true)
        );
        if !changed { return false; }
        if keys.just_pressed(KeyCode::Digit1) { *mode = 1; }
        if keys.just_pressed(KeyCode::Digit2) { *mode = 2; }
        if keys.just_pressed(KeyCode::Digit3) { *mode = 3; }
        true
    };

    if update_scales(&mut display.mode) {
        // Collect parent atom scales
        let n = atom_info.elements.len();
        let mut scales: Vec<f32> = vec![1.0; n];
        for i in 0..n {
            scales[i] = match display.mode {
                1 => 1.0,   // ball-stick: uniform viewing size
                2 => resources::ionic_radius(&atom_info.elements[i]) / 0.5,
                3 => resources::atomic_radius(&atom_info.elements[i]) / 0.5,
                _ => 1.0,
            };
        }
        for (i, mut t) in atoms.iter_mut().enumerate() {
            if i < picking.parent_indices.len() {
                let p = picking.parent_indices[i];
                if p < scales.len() {
                    t.scale = Vec3::splat(scales[p]);
                }
            }
        }
        // Show/hide bonds per mode
        display.show_bonds = display.mode != 2;
        let vis = if display.show_bonds { Visibility::Visible } else { Visibility::Hidden };
        for e in bonds.iter() { commands.entity(e).insert(vis); }
    }
    if keys.just_pressed(KeyCode::KeyB) {
        display.show_bonds = !display.show_bonds;
        let vis = if display.show_bonds { Visibility::Visible } else { Visibility::Hidden };
        for e in bonds.iter() { commands.entity(e).insert(vis); }
    }
    if keys.just_pressed(KeyCode::KeyC) {
        display.show_cell = !display.show_cell;
        let vis = if display.show_cell { Visibility::Visible } else { Visibility::Hidden };
        for e in cells.iter() { commands.entity(e).insert(vis); }
    }
}

fn default_cu_fcc() -> CrystalData {
    crystal::CrystalData {
        lattice: crystal::Lattice { a: 3.615, b: 3.615, c: 3.615, alpha: 90.0, beta: 90.0, gamma: 90.0 },
        atoms: vec![
            crystal::AtomData { element: "Cu".into(), x: 0.0, y: 0.0, z: 0.0, label: "1".into() },
            crystal::AtomData { element: "Cu".into(), x: 1.8075, y: 1.8075, z: 0.0, label: "2".into() },
            crystal::AtomData { element: "Cu".into(), x: 1.8075, y: 0.0, z: 1.8075, label: "3".into() },
            crystal::AtomData { element: "Cu".into(), x: 0.0, y: 1.8075, z: 1.8075, label: "4".into() },
        ],
        positions_fractional: false,
        modified: false,
    }
}
