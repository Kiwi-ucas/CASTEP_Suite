//! egui UI panels: atom list, info panel, toolbar

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use crate::picking::PickingState;
use crate::MoveState;
use crate::AddAtomState;
use crate::CrystalStore;
use crate::Lattice;
use crate::PhononState;
use crate::RotateState;
use crate::PanelRects;
use crate::PesState;
use crate::Pes3dState;
use crate::DisplayMode;
use crate::VisMode;
use crate::IsoMaterial;
use crate::CubeResource;
use crate::RenderSettings;
use crate::resources;
use crate::SlabState;

#[derive(Resource)]
pub struct AtomInfo {
    pub elements: Vec<String>,
    pub labels: Vec<String>,
    pub radii: Vec<f32>,
}

#[derive(Resource)]
pub struct CrystalMeta {
    pub filename: String,
    pub a: f32, pub b: f32, pub c: f32,
    pub alpha: f32, pub beta: f32, pub gamma: f32,
}

impl AtomInfo {
    pub fn new(elements: Vec<String>, labels: Vec<String>, radii: Vec<f32>) -> Self {
        Self { elements, labels, radii }
    }
}

pub fn ui_system(
    mut contexts: EguiContexts,
    mut picking: ResMut<PickingState>,
    mut atom_info: ResMut<AtomInfo>,
    meta: Option<Res<CrystalMeta>>,
    mut move_state: ResMut<MoveState>,
    mut add_state: ResMut<AddAtomState>,
    crystal: Option<Res<CrystalStore>>,
    mut phonon_state: Option<ResMut<PhononState>>,
    mut rotate_state: ResMut<RotateState>,
    mut panel_rects: ResMut<PanelRects>,
    pes_state: Option<Res<PesState>>,
    mut pes3d_state: Option<ResMut<Pes3dState>>,
    cube: Option<Res<CubeResource>>,
    mut display: ResMut<DisplayMode>,
    mut render_settings: ResMut<RenderSettings>,
    mut slab_state: ResMut<SlabState>,
) {
    let Some(ctx) = contexts.try_ctx_mut() else { return; };

    // Pre-compute: parent indices for selected/hovered images
    let sel_parent = if picking.selected >= 0 && (picking.selected as usize) < picking.parent_indices.len() {
        Some(picking.parent_indices[picking.selected as usize])
    } else { None };
    let hov_parent = if picking.hovered >= 0 && (picking.hovered as usize) < picking.parent_indices.len() {
        Some(picking.parent_indices[picking.hovered as usize])
    } else { None };

    // Pre-compute first expanded-atom index for each asymmetric-unit atom
    let count = atom_info.elements.len();
    let first_images: Vec<Option<usize>> = (0..count)
        .map(|i| picking.parent_indices.iter().position(|&p| p == i))
        .collect();

    // ── Bottom toolbar ──
    let bottom_resp = egui::TopBottomPanel::bottom("toolbar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            if let Some(p) = sel_parent {
                let el = if p < atom_info.elements.len() {
                    atom_info.elements[p].as_str()
                } else { "?" };
                let img_count = picking.parent_indices.iter().filter(|&&x| x == p).count();
                ui.label(egui::RichText::new(
                    format!("Atom {} ({}) selected, {} images   H/K:\u{b1}X  U/M:\u{b1}Y  I/N:\u{b1}Z   D: delete",
                        p, el, img_count)
                ).strong());
            } else {
                let mut base = "\u{1f5b0} Right-drag: rotate | Scroll: zoom | Click: select | 1/2/3: mode | P: proj | A: axes | B: bonds | C: cell | H: atoms | R: reset".to_string();
                if let Some(ref pes) = pes_state {
                    if pes.has_energies {
                        base.push_str(" | +/-: color range");
                    }
                }
                ui.label(base);
            }
        });
    });
    panel_rects.bottom = Some(bottom_resp.response.rect);

    // ── Left panel: atom list (asymmetric unit only) ──
    let left_resp = egui::SidePanel::left("atom_list")
        .resizable(false).default_width(205.0)
        .show(ctx, |ui| {
            ui.heading("Atoms");
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("\u{2795} Add Atom").clicked() {
                    add_state.show_table = true;
                }
                let hide_label = if display.show_atoms { "Hide Atoms" } else { "Show Atoms" };
                if ui.button(hide_label).clicked() {
                    display.show_atoms = !display.show_atoms;
                }
            });
            ui.separator();
            egui::ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                for i in 0..count {
                    let label = format!("{:2}. {:3}", i, atom_info.elements[i]);
                    let selected = sel_parent == Some(i);
                    let hovered = hov_parent == Some(i);

                    let mut text = egui::RichText::new(&label);
                    if selected {
                        text = text.color(egui::Color32::from_rgb(255, 170, 30)).strong();
                    } else if hovered {
                        text = text.color(egui::Color32::from_rgb(180, 180, 180));
                    }
                    ui.horizontal(|ui| {
                        let response = ui.selectable_label(selected, text);
                        if response.clicked() {
                            if let Some(img) = first_images[i] {
                                picking.selected = img as i32;
                            }
                        }
                        if response.hovered() {
                            if let Some(img) = first_images[i] {
                                picking.hovered = img as i32;
                            }
                        }
                        // Radius edit
                        if i < atom_info.radii.len() {
                            let mut r = atom_info.radii[i];
                            let resp = ui.add(
                                egui::DragValue::new(&mut r)
                                    .speed(0.01)
                                    .range(0.1..=5.0)
                                    .suffix(" \u{c5}")
                            );
                            if resp.changed() && i < atom_info.radii.len() {
                                atom_info.radii[i] = r;
                                picking.modified = true;
                            }
                        }
                    });
                }
            });
        });
    panel_rects.left = Some(left_resp.response.rect);

    // ── Right panel ──
    let right_resp = egui::SidePanel::right("info_panel")
        .resizable(false).default_width(220.0)
        .show(ctx, |ui| {
            // Crystal info header
            if let Some(m) = meta.as_ref() {
                ui.heading(&m.filename);
                ui.separator();
                ui.label("Cell Parameters:");
                ui.label(format!("a = {:.3} \u{c5}", m.a));
                ui.label(format!("b = {:.3} \u{c5}", m.b));
                ui.label(format!("c = {:.3} \u{c5}", m.c));
                ui.label(format!("\u{3b1} = {:.2}\u{b0}", m.alpha));
                ui.label(format!("\u{3b2} = {:.2}\u{b0}", m.beta));
                ui.label(format!("\u{3b3} = {:.2}\u{b0}", m.gamma));
                // Cell volume from lattice
                if let Some(c) = crystal.as_ref() {
                    let vol = c.data.lattice.cell_volume();
                    ui.label(format!("Volume: {:.2} \u{c5}\u{b3}", vol));
                }
                ui.label(format!("Asym. atoms: {}", atom_info.elements.len()));
                ui.label(format!("Displayed: {}", picking.atom_positions.len()));
                ui.separator();

                // ── Slab / vacuum / supercell: settings live in popup windows ──
                ui.horizontal(|ui| {
                    if ui.button("Slab cut").clicked() {
                        slab_state.slab_open = true;
                    }
                    if ui.button("Vacuum").clicked() {
                        slab_state.vacuum_open = true;
                    }
                    if ui.button("Supercell").clicked() {
                        slab_state.supercell_open = true;
                    }
                });
                let hkl_ok: Option<(i32, i32, i32)> = match (
                    slab_state.h_str.trim().parse::<i32>(),
                    slab_state.k_str.trim().parse::<i32>(),
                    slab_state.l_str.trim().parse::<i32>(),
                ) {
                    (Ok(h), Ok(k), Ok(l)) if !(h == 0 && k == 0 && l == 0) => Some((h, k, l)),
                    _ => None,
                };

                // ── Dual-input synchronisation (last-edited field wins) ──
                if let Some(c) = crystal.as_ref() {
                    let hkl = hkl_ok;
                    if let Some((h, k, l)) = hkl {
                        if slab_state.last_hkl != (h, k, l) {
                            slab_state.last_hkl = (h, k, l);
                            // Miller indices changed: re-derive the fractional
                            // displays from the current absolute values.
                            if let Some(per) = crate::slab::slab_period(&c.data, h, k, l) {
                                let s0 = slab_state.start_str.trim().parse::<f32>().unwrap_or(0.0);
                                slab_state.s_frac_str = if s0.abs() < 1e-6 {
                                    String::new()
                                } else {
                                    format!("{:.4}", s0 / per)
                                };
                                slab_state.t_frac_str = if slab_state.thickness <= 0.0 {
                                    String::new()
                                } else {
                                    format!("{:.4}", slab_state.thickness / per)
                                };
                                slab_state.s_frac_dirty = false;
                                slab_state.start_dirty = false;
                                slab_state.t_frac_dirty = false;
                                slab_state.thick_dirty = false;
                            }
                        }
                        if let Some(period) = crate::slab::slab_period(&c.data, h, k, l) {
                            if slab_state.s_frac_dirty {
                                if let Ok(f) = slab_state.s_frac_str.trim().parse::<f32>() {
                                    slab_state.start_str = format!("{:.4}", f * period);
                                }
                                slab_state.s_frac_dirty = false;
                            } else if slab_state.start_dirty {
                                let a = slab_state.start_str.trim().parse::<f32>().unwrap_or(0.0);
                                slab_state.s_frac_str = if a.abs() < 1e-6 {
                                    String::new()
                                } else {
                                    format!("{:.4}", a / period)
                                };
                                slab_state.start_dirty = false;
                            }
                            if slab_state.t_frac_dirty {
                                let t = slab_state.t_frac_str.trim();
                                slab_state.thickness = if t.is_empty() {
                                    0.0
                                } else {
                                    t.parse::<f32>().unwrap_or(0.0) * period
                                };
                                slab_state.t_frac_dirty = false;
                            } else if slab_state.thick_dirty {
                                slab_state.t_frac_str = if slab_state.thickness <= 0.0 {
                                    String::new()
                                } else {
                                    format!("{:.4}", slab_state.thickness / period)
                                };
                                slab_state.thick_dirty = false;
                            }
                        }
                        // V ↔ c coupling: c = T_eff + V (last edited wins)
                        let [va, vb, vc] = c.data.lattice.to_vectors();
                        let axis_len = [va.length(), vb.length(), vc.length()]
                            .get((slab_state.vac_axis - 1) as usize)
                            .copied()
                            .unwrap_or(0.0);
                        let old_vac = if c.data.vacuum.as_ref().is_some_and(|vm| vm.axis == slab_state.vac_axis) {
                            c.data.vacuum.as_ref().unwrap().thickness_ang
                        } else {
                            0.0
                        };
                        let base = axis_len - old_vac;
                        // V ↔ c coupling: c = base + V (last edited wins).
                        // Auto mode (vac_cell <= 0): c follows V implicitly,
                        // the c field keeps showing "Auto"; only a manually
                        // typed c value pins it (and V then follows c).
                        if slab_state.vac_v_dirty {
                            if slab_state.vac_cell > 0.0 {
                                slab_state.vac_cell = base + slab_state.vac_thickness;
                                slab_state.vac_c_str = format!("{:.1}", slab_state.vac_cell);
                            } else {
                                slab_state.vac_c_str.clear(); // stay on "Auto"
                            }
                            slab_state.vac_v_dirty = false;
                        } else if slab_state.vac_c_dirty {
                            if slab_state.vac_cell > 0.0 {
                                if slab_state.vac_cell < base {
                                    slab_state.vac_cell = base; // V ≥ 0
                                    slab_state.vac_c_str = format!("{:.1}", base);
                                }
                                slab_state.vac_thickness = (slab_state.vac_cell - base).max(0.0);
                            } else {
                                slab_state.vac_c_str.clear(); // back to auto
                            }
                            slab_state.vac_c_dirty = false;
                        }
                    }
                }

                if !slab_state.error.is_empty() {
                    ui.label(egui::RichText::new(&slab_state.error)
                        .color(egui::Color32::from_rgb(255, 90, 90)));
                }
                // Applied slab / vacuum provenance
                if let Some(c) = crystal.as_ref() {
                    if let Some(ref sm) = c.data.slab {
                        let b = match sm.basis {
                            crate::InPlaneBasis::Primitive => "prim",
                            crate::InPlaneBasis::Orthogonal => "90\u{b0}",
                            crate::InPlaneBasis::Conventional => "conv",
                        };
                        let uv_extra = match (&sm.u_vec, &sm.v_vec) {
                            (Some(uv), Some(vv)) => Some(format!(
                                " U=({} {} {}) V=({} {} {})",
                                uv[0], uv[1], uv[2], vv[0], vv[1], vv[2]
                            )),
                            _ => None,
                        };
                        ui.label(format!(
                            "Slab: ({},{},{}) s={:.2} T={:.2} \u{c5} {} \u{00d7} {} [{}]{}",
                            sm.h, sm.k, sm.l, sm.start_ang, sm.thickness_ang, sm.u, sm.v, b,
                            uv_extra.as_deref().unwrap_or("")
                        ));
                    }
                    if let Some(ref vm) = c.data.vacuum {
                        let ax = ["?", "a", "b", "c"].get((vm.axis as usize).clamp(0, 3)).copied().unwrap_or("?");
                        let pos = if vm.position > 0.0 {
                            format!("{:.2}", vm.position)
                        } else if vm.both_sides {
                            "center".into()
                        } else {
                            "top".into()
                        };
                        ui.label(format!(
                            "Vacuum: {} +{:.1} \u{c5} (pos {})",
                            ax,
                            vm.thickness_ang,
                            pos
                        ));
                    }
                    if let Some(ref sc) = c.data.supercell {
                        ui.label(format!(
                            "Supercell: {}\u{00d7}{}\u{00d7}{}",
                            sc[0], sc[1], sc[2]
                        ));
                    }
                }
                ui.separator();
            } else {
                ui.heading("Crystal Viewer");
            }

            // Selected atom details
            if let Some(p) = sel_parent {
                if p < atom_info.elements.len() {
                    ui.label(egui::RichText::new("Selected Atom").strong());
                    ui.label(format!("Element: {}", atom_info.elements[p]));
                    ui.label(format!("Asym. index: {}", p));
                    if picking.selected >= 0 && (picking.selected as usize) < picking.atom_positions.len() {
                        let cart = picking.atom_positions[picking.selected as usize];
                        ui.label("Cartesian (\u{c5}):");
                        ui.label(format!("  X: {:.4}", cart.x));
                        ui.label(format!("  Y: {:.4}", cart.y));
                        ui.label(format!("  Z: {:.4}", cart.z));
                        // Show fractional coords if lattice data available
                        if let Some(c) = crystal.as_ref() {
                            let inv = c.data.lattice.inverse_vectors();
                            let frac = Lattice::apply_inverse(&inv, cart);
                            ui.label("Fractional:");
                            ui.label(format!("  x: {:.4}", frac.x));
                            ui.label(format!("  y: {:.4}", frac.y));
                            ui.label(format!("  z: {:.4}", frac.z));
                        }
                    }
                }
            } else {
                ui.label("Click an atom\nto see details");
            }

            // ── PES scan info (2D) ──
            if let Some(ref pes) = pes_state {
                ui.separator();
                ui.label(egui::RichText::new("PES Scan").strong());
                ui.label(format!("Plane: {}", pes.plane.to_uppercase()));
                ui.label(format!("Grid: {}×{}", pes.nx, pes.ny));
                ui.label(format!("Mode: {}", pes.scan_mode));
                if pes.has_energies {
                    ui.label(format!("E range: {:.4} – {:.4} eV", pes.e_min, pes.e_max));
                    let pct = (crate::step_to_clip(pes.color_step) * 100.0 + 0.5) as u32;
                    ui.label(format!("Color range: {}% (+/-)", pct));
                    ui.label(format!("Surface: {}", if pes.show_surface { "ON (S)" } else { "OFF (S)" }));
                } else {
                    ui.label("No energies (run CASTEP first)");
                }
            }

            // ── PES 3D scan info ──
            if let Some(ref ps) = pes3d_state {
                ui.separator();
                ui.label(egui::RichText::new("PES 3D Scan").strong());
                ui.label(format!("Grid: {}×{}×{}", ps.nx, ps.ny, ps.nz));
                if ps.has_energies {
                    ui.label(format!("E range: {:.2} – {:.2} eV", ps.e_min, ps.e_max));
                    let mode_str = match ps.vis_mode {
                        VisMode::Isosurface => format!("Isosurface (4) iso={:.2}", ps.iso_value),
                        VisMode::Volume => format!("Volume (5) layer={:.1}", ps.vol_iso_ref),
                        VisMode::Slice => format!("Slice (6) axis={} pos={:.2}", ps.slice_axis, ps.slice_pos),
                        VisMode::Sphere => format!("Sphere (7) R={:.2} Å", ps.sphere_radius),
                        VisMode::Migration => format!("Migration (8) cap={:.2} eV", ps.mig_e_cap),
                    };
                    ui.label(format!("Mode: {}", mode_str));
                } else {
                    ui.label("No energies (run CASTEP first)");
                }
            }

            // ── Isosurface controls (only when in isosurface mode) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && ps.vis_mode == VisMode::Isosurface {
                    ui.horizontal(|ui| {
                        ui.label("Isosurface step");
                        ui.add(egui::DragValue::new(&mut ps.iso_step)
                            .speed(0.1)
                            .range(0.01..=1000.0)
                            .suffix(" eV"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Material");
                        egui::ComboBox::from_id_salt("iso_material")
                            .selected_text(ps.iso_material.name())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Opaque, "Opaque");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::SemiTransparent, "Semi-transparent");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Transparent, "Transparent");
                            });
                    });
                }
            }

            // ── Sphere section controls (mode 7) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && ps.vis_mode == VisMode::Sphere {
                    ui.separator();
                    ui.label(egui::RichText::new("Sphere Section").strong());
                    if let Some(ref cube) = cube {
                        let atoms = &cube.0.atoms;
                        let sel = if ps.sphere_center_idx == usize::MAX {
                            "Custom…".to_string()
                        } else {
                            let i = ps.sphere_center_idx.min(atoms.len().saturating_sub(1));
                            format!("{}: {}", i, crate::atom_z_to_symbol(atoms[i].z))
                        };
                        ui.horizontal(|ui| {
                            ui.label("Center");
                            egui::ComboBox::from_id_salt("sphere_center")
                                .selected_text(sel)
                                .show_ui(ui, |ui| {
                                    for (i, a) in atoms.iter().enumerate() {
                                        ui.selectable_value(
                                            &mut ps.sphere_center_idx, i,
                                            format!("{}: {} ({:.2},{:.2},{:.2})",
                                                i, crate::atom_z_to_symbol(a.z), a.x, a.y, a.z_coord));
                                    }
                                    ui.selectable_value(&mut ps.sphere_center_idx, usize::MAX, "Custom…");
                                });
                        });
                        if ps.sphere_center_idx == usize::MAX {
                            ui.horizontal(|ui| {
                                ui.label("Frac");
                                ui.add(egui::DragValue::new(&mut ps.sphere_center_custom[0])
                                    .speed(0.005).range(0.0..=1.0));
                                ui.add(egui::DragValue::new(&mut ps.sphere_center_custom[1])
                                    .speed(0.005).range(0.0..=1.0));
                                ui.add(egui::DragValue::new(&mut ps.sphere_center_custom[2])
                                    .speed(0.005).range(0.0..=1.0));
                            });
                        }
                    }
                    ui.horizontal(|ui| {
                        ui.label("Radius");
                        ui.add(egui::DragValue::new(&mut ps.sphere_radius)
                            .speed(0.05)
                            .range(0.1..=5.0)
                            .suffix(" Å"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Material");
                        egui::ComboBox::from_id_salt("sphere_material")
                            .selected_text(ps.iso_material.name())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Opaque, "Opaque");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::SemiTransparent, "Semi-transparent");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Transparent, "Transparent");
                            });
                    });
                }
            }

            // ── Migration surface controls (mode 8) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && ps.vis_mode == VisMode::Migration {
                    ui.separator();
                    ui.label(egui::RichText::new("Migration Surface").strong());
                    ui.label("Cage shells = 1st radial min, welded across window bulges");
                    let e_max = ps.e_max.max(1.0);
                    ui.horizontal(|ui| {
                        ui.label("E cap");
                        ui.add(egui::Slider::new(&mut ps.mig_e_cap, 0.05..=e_max)
                            .suffix(" eV"));
                    });
                    ui.checkbox(&mut ps.mig_show_shell, "Cage shells");
                    ui.horizontal(|ui| {
                        ui.label("Material");
                        egui::ComboBox::from_id_salt("mig_material")
                            .selected_text(ps.iso_material.name())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Opaque, "Opaque");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::SemiTransparent, "Semi-transparent");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Transparent, "Transparent");
                            });
                    });
                    if let Some(ref cube) = cube {
                        let n = crate::sphere_section::detect_cage_centers(&cube.0).len();
                        ui.label(format!("Cage centers (auto): {}", n));
                    }
                }
            }

            // ── Volume render controls (only when in volume mode) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && ps.vis_mode == VisMode::Volume {
                    ui.separator();
                    ui.label(egui::RichText::new("Volume").strong());
                    ui.horizontal(|ui| {
                        ui.label("Opacity");
                        ui.add(egui::Slider::new(&mut ps.alpha_scale, 0.1..=2.0)
                            .logarithmic(true));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Falloff");
                        ui.add(egui::Slider::new(&mut ps.alpha_falloff, 0.05..=1.0));
                    });
                    let e_min = ps.e_min;
                    let e_max = ps.e_max;
                    ui.horizontal(|ui| {
                        ui.label("Layer (eV)");
                        ui.add(egui::DragValue::new(&mut ps.vol_iso_ref)
                            .speed(1.0)
                            .range(e_min..=e_max)
                            .suffix(" eV"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Quality");
                        egui::ComboBox::from_id_salt("vol_quality")
                            .selected_text(format!("{} steps", ps.vol_steps))
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut ps.vol_steps, 64u32, "64 steps");
                                ui.selectable_value(&mut ps.vol_steps, 128u32, "128 steps");
                                ui.selectable_value(&mut ps.vol_steps, 192u32, "192 steps");
                            });
                    });
                }
            }

            // ── Color mapping range (for all 3D PES modes) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies {
                    ui.separator();
                    ui.label(egui::RichText::new("Color Mapping").strong());
                    let e_min = ps.e_min;
                    let e_max = ps.e_max;
                    let color_min = ps.color_min;
                    let color_max = ps.color_max;
                    ui.horizontal(|ui| {
                        ui.label("Min");
                        ui.add(egui::DragValue::new(&mut ps.color_min)
                            .speed(1.0)
                            .range(e_min..=color_max)
                            .suffix(" eV"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Max");
                        ui.add(egui::DragValue::new(&mut ps.color_max)
                            .speed(1.0)
                            .range(color_min..=e_max)
                            .suffix(" eV"));
                    });
                    if ui.button("Reset to full range").clicked() {
                        ps.color_min = ps.e_min;
                        ps.color_max = ps.e_max;
                    }
                }
            }

            // ── Spatial clipping (XYZ ranges) — compact: one row per axis ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && matches!(ps.vis_mode, VisMode::Isosurface | VisMode::Volume) {
                    ui.separator();
                    ui.label(egui::RichText::new("Spatial Clipping").strong());

                    let clip_x = ps.clip_x;
                    let clip_y = ps.clip_y;
                    let clip_z = ps.clip_z;

                    // Each row: axis label + "min"/"max" labels OUTSIDE the
                    // DragValue (labels are not part of the editable value)
                    ui.horizontal(|ui| {
                        ui.label("X");
                        ui.label("min");
                        ui.add(egui::DragValue::new(&mut ps.clip_x[0])
                            .speed(0.01).range(0.0..=clip_x[1]));
                        ui.label("max");
                        ui.add(egui::DragValue::new(&mut ps.clip_x[1])
                            .speed(0.01).range(clip_x[0]..=1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Y");
                        ui.label("min");
                        ui.add(egui::DragValue::new(&mut ps.clip_y[0])
                            .speed(0.01).range(0.0..=clip_y[1]));
                        ui.label("max");
                        ui.add(egui::DragValue::new(&mut ps.clip_y[1])
                            .speed(0.01).range(clip_y[0]..=1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Z");
                        ui.label("min");
                        ui.add(egui::DragValue::new(&mut ps.clip_z[0])
                            .speed(0.01).range(0.0..=clip_z[1]));
                        ui.label("max");
                        ui.add(egui::DragValue::new(&mut ps.clip_z[1])
                            .speed(0.01).range(clip_z[0]..=1.0));
                    });

                    if ui.button("Reset clipping").clicked() {
                        ps.clip_x = [0.0, 1.0];
                        ps.clip_y = [0.0, 1.0];
                        ps.clip_z = [0.0, 1.0];
                    }
                }
            }

            ui.with_layout(egui::Layout::bottom_up(egui::Align::Min), |ui| {
                ui.horizontal(|ui| {
                    ui.label("Move Step");
                    ui.add(egui::DragValue::new(&mut move_state.step)
                        .speed(0.05)
                        .suffix(" \u{c5}"));
                });
                ui.horizontal(|ui| {
                    ui.label("Rotation Angle");
                    ui.add(egui::DragValue::new(&mut rotate_state.angle_deg)
                        .suffix("\u{b0}"));
                });
                ui.horizontal(|ui| {
                    if ui.button("Render").clicked() {
                        render_settings.show_dialog = true;
                    }
                    // Export PLY lives next to Render (mode 8 with a cube only)
                    let mig_ok = pes3d_state.as_ref()
                        .map(|p| p.has_energies && p.vis_mode == VisMode::Migration)
                        .unwrap_or(false);
                    if mig_ok && cube.is_some() {
                        if ui.button("Export PLY").clicked() || render_settings.request_ply {
                            if let Some(ps) = pes3d_state.as_ref() {
                                export_migration_ply(&cube, ps);
                            }
                            render_settings.request_ply = false;
                        }
                    }
                    if render_settings.rendering.load(std::sync::atomic::Ordering::Relaxed) {
                        ui.label(egui::RichText::new("Rendering…")
                            .color(egui::Color32::LIGHT_BLUE));
                    }
                });
                ui.separator();
            });
        });
    panel_rects.right = Some(right_resp.response.rect);

    // ── Slab / vacuum settings popups (opened from the right panel buttons) ──
    // Same pattern as the Render dialog: egui::Area + window Frame — the popup
    // is draggable from anywhere on its surface and remembers its position.
    if slab_state.slab_open || slab_state.vacuum_open || slab_state.supercell_open {
        if std::env::var("CRYSTAL_VIEWER_DEBUG_UI").is_ok() {
            eprintln!(
                "[debug-ui] popup branch: slab_open={} vacuum_open={} supercell_open={} screen={:?}",
                slab_state.slab_open,
                slab_state.vacuum_open,
                slab_state.supercell_open,
                ctx.screen_rect()
            );
        }
        let sr = ctx.screen_rect();
        let pop_x = (sr.max.x - 240.0 - 330.0).max(sr.min.x + 10.0);
        if slab_state.slab_open {
            let slab_resp = egui::Area::new(egui::Id::new("slab_settings"))
                .default_pos(egui::pos2(pop_x, 60.0))
                .movable(true)
                .show(ctx, |ui| {
                    egui::Frame::window(ui.style()).show(ui, |ui| {
                        // The window hugs its content (no forced min width,
                        // no inner scroll area) — everything shows at once.
                        ui.horizontal(|ui| {
                            ui.label(egui::RichText::new("Slab cut (hkl)").strong());
                            // Plain left-aligned close button: a right_to_left
                            // sub-layout claims the whole available width on
                            // the sizing pass, which would stretch the
                            // content-hugging frame to the screen edge.
                            if ui.small_button("\u{00d7}").clicked() {
                                slab_state.slab_open = false;
                            }
                        });
                        // Bounded separator: an unbounded one would stretch
                        // the frame to the whole available (screen) width,
                        // defeating the content-hugging window size.
                        ui.separator();
                        slab_section_ui(ui, &mut slab_state, crystal.as_deref());
                    });
                });
            let slab_rect = slab_resp.response.rect;
            if std::env::var("CRYSTAL_VIEWER_DEBUG_UI").is_ok() {
                eprintln!("[debug-ui] slab popup rect={:?}", slab_rect);
            }
        }
        if slab_state.vacuum_open {
            let vac_resp = egui::Area::new(egui::Id::new("vacuum_settings"))
                .default_pos(egui::pos2(pop_x, 420.0))
                .movable(true)
                .show(ctx, |ui| {
                    egui::Frame::window(ui.style()).show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(egui::RichText::new("Vacuum layer").strong());
                            if ui.small_button("\u{00d7}").clicked() {
                                slab_state.vacuum_open = false;
                            }
                        });
                        ui.separator();
                        vacuum_section_ui(ui, &mut slab_state);
                    });
                });
            let vac_rect = vac_resp.response.rect;
            if std::env::var("CRYSTAL_VIEWER_DEBUG_UI").is_ok() {
                eprintln!("[debug-ui] vacuum popup rect={:?}", vac_rect);
            }
        }
        if slab_state.supercell_open {
            let sc_resp = egui::Area::new(egui::Id::new("supercell_settings"))
                .default_pos(egui::pos2(pop_x, 240.0))
                .movable(true)
                .show(ctx, |ui| {
                    egui::Frame::window(ui.style()).show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(egui::RichText::new("Supercell (nx ny nz)").strong());
                            // Plain left-aligned close button (a right_to_left
                            // sub-layout would stretch the frame to the
                            // screen edge — see the slab popup comment).
                            if ui.small_button("\u{00d7}").clicked() {
                                slab_state.supercell_open = false;
                            }
                        });
                        ui.separator();
                        supercell_section_ui(ui, &mut slab_state, crystal.as_deref());
                    });
                });
            let sc_rect = sc_resp.response.rect;
            if std::env::var("CRYSTAL_VIEWER_DEBUG_UI").is_ok() {
                eprintln!("[debug-ui] supercell popup rect={:?}", sc_rect);
            }
        }
    }

    // ── Render dialog (Bevy-native offscreen export — WYSIWYG) ──
    // Built on egui::Area instead of egui::Window: the Area is movable from
    // ANYWHERE on its surface (not just a title bar) and remembers its
    // position across frames.
    if render_settings.show_dialog {
        let sr = ctx.screen_rect();
        egui::Area::new(egui::Id::new("render_dialog"))
            .default_pos(egui::pos2(sr.center().x - 180.0, sr.center().y - 240.0))
            .movable(true)
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).show(ui, |ui| {
                ui.set_min_width(340.0);
                ui.label(egui::RichText::new("Render").strong());
                ui.separator();
                if render_settings.rendering.load(std::sync::atomic::Ordering::Relaxed) {
                    ui.colored_label(egui::Color32::LIGHT_BLUE, "Rendering…");
                }
                let rs = &mut *render_settings;
                ui.horizontal(|ui| {
                    ui.label("Resolution");
                    ui.add(egui::DragValue::new(&mut rs.width).range(320..=8192).suffix(" x"));
                    ui.add(egui::DragValue::new(&mut rs.height).range(240..=8192));
                });
                // NOTE: 8x is deliberately NOT offered — WebGPU only
                // guarantees [1,4] samples for Rgba8UnormSrgb, and Apple
                // silicon (Metal) supports at most 4x. Requesting 8x aborts
                // the whole viewer (wgpu validation error).
                ui.horizontal(|ui| {
                    ui.label("MSAA");
                    egui::ComboBox::from_id_salt("render_msaa")
                        .selected_text(format!("{}x", rs.msaa_samples))
                        .show_ui(ui, |ui| {
                            for s in [1u32, 2, 4] {
                                ui.selectable_value(&mut rs.msaa_samples, s, if s == 1 { "Off (1x)".to_string() } else { format!("{s}x") });
                            }
                        });
                });
                ui.horizontal(|ui| {
                    ui.label("Format");
                    egui::ComboBox::from_id_salt("render_format")
                        .selected_text(rs.format.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut rs.format, crate::render_export::ImgFormat::Png, "PNG");
                            ui.selectable_value(&mut rs.format, crate::render_export::ImgFormat::Tiff, "TIFF");
                        });
                });
                ui.separator();
                ui.label(egui::RichText::new("Scene parameters (live — edits apply to the viewer immediately)").strong());
                ui.horizontal(|ui| {
                    ui.label("Key light");
                    ui.add(egui::Slider::new(&mut rs.key_lux, 0.0..=20000.0).suffix(" lx"));
                });
                ui.horizontal(|ui| {
                    ui.label("Fill light");
                    ui.add(egui::Slider::new(&mut rs.fill_lux, 0.0..=20000.0).suffix(" lx"));
                });
                ui.horizontal(|ui| {
                    ui.label("Ambient");
                    ui.add(egui::Slider::new(&mut rs.ambient_lux, 0.0..=2000.0).suffix(" lx"));
                });
                ui.checkbox(&mut rs.shadows_enabled, "Shadows (key light)");
                ui.horizontal(|ui| {
                    ui.label("Atom roughness");
                    ui.add(egui::Slider::new(&mut rs.atom_roughness, 0.0..=1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("Atom metallic");
                    ui.add(egui::Slider::new(&mut rs.atom_metallic, 0.0..=1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("Tonemapping");
                    egui::ComboBox::from_id_salt("render_tonemap")
                        .selected_text(rs.tonemap.label())
                        .show_ui(ui, |ui| {
                            for m in [crate::render_export::TonemapChoice::TonyMcMapface,
                                      crate::render_export::TonemapChoice::Aces,
                                      crate::render_export::TonemapChoice::AgX,
                                      crate::render_export::TonemapChoice::Reinhard] {
                                ui.selectable_value(&mut rs.tonemap, m, m.label());
                            }
                        });
                });
                ui.separator();
                ui.label(egui::RichText::new("Background color").strong());
                let mut br = (rs.bg_r * 255.0).round() as u8;
                let mut bg = (rs.bg_g * 255.0).round() as u8;
                let mut bb = (rs.bg_b * 255.0).round() as u8;
                ui.horizontal(|ui| { ui.label("R"); ui.add(egui::Slider::new(&mut br, 0..=255)); });
                ui.horizontal(|ui| { ui.label("G"); ui.add(egui::Slider::new(&mut bg, 0..=255)); });
                ui.horizontal(|ui| { ui.label("B"); ui.add(egui::Slider::new(&mut bb, 0..=255)); });
                rs.bg_r = br as f32 / 255.0;
                rs.bg_g = bg as f32 / 255.0;
                rs.bg_b = bb as f32 / 255.0;
                ui.horizontal(|ui| {
                    let presets: [(&str, [f32; 3]); 5] = [
                        ("Black", [0.0, 0.0, 0.0]),
                        ("White", [1.0, 1.0, 1.0]),
                        ("Light gray", [170.0 / 255.0, 170.0 / 255.0, 170.0 / 255.0]),
                        ("Dark gray", [0.25, 0.25, 0.25]),
                        ("Viewer default", [43.0 / 255.0, 44.0 / 255.0, 47.0 / 255.0]),
                    ];
                    for (label, rgb) in presets {
                        if ui.button(label).clicked() {
                            rs.bg_r = rgb[0]; rs.bg_g = rgb[1]; rs.bg_b = rgb[2];
                        }
                    }
                });
                let status = rs.last_status.lock().map(|s| s.clone()).unwrap_or_default();
                if !status.is_empty() {
                    ui.label(egui::RichText::new(&status).color(egui::Color32::from_gray(200)));
                }
                ui.separator();
                let busy = rs.rendering.load(std::sync::atomic::Ordering::Relaxed);
                ui.horizontal(|ui| {
                    if ui.add_enabled(!busy, egui::Button::new(
                            format!("Render -> render.{}", rs.format.ext()))).clicked() {
                        rs.request = true;
                    }
                    if ui.button("Reset all").clicked() {
                        rs.reset_params();
                    }
                    if ui.button("Close").clicked() {
                        rs.show_dialog = false;
                    }
                });
                });
            });
    }

    // ── Add Atom popup: periodic table ──
    if add_state.show_table {
        egui::Window::new("Periodic Table — Select Element")
            .collapsible(false).resizable(false)
            .default_width(800.0).default_height(400.0)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                let mut pt_map = std::collections::HashMap::new();
                for &(sym, r, c) in resources::PERIODIC_TABLE {
                    pt_map.insert((r, c), sym);
                }
                let rows = 9;
                let cols = 18;
                let cell_size = 36.0;
                ui.vertical_centered(|ui| {
                    ui.heading("Periodic Table of Elements");
                });
                ui.add_space(6.0);
                for r in 0..rows {
                    ui.horizontal(|ui| {
                        for c in 0..cols {
                            if let Some(sym) = pt_map.get(&(r, c)) {
                                let bg = resources::category_color(sym);
                                let rich = egui::RichText::new(*sym)
                                    .size(11.0)
                                    .color(egui::Color32::WHITE);
                                let btn = egui::Button::new(rich)
                                    .fill(bg)
                                    .min_size(egui::vec2(cell_size, cell_size));
                                if ui.add(btn).clicked() {
                                    add_state.selected_element = Some(sym.to_string());
                                    add_state.show_table = false;
                                    add_state.coord_x = "0.0".into();
                                    add_state.coord_y = "0.0".into();
                                    add_state.coord_z = "0.0".into();
                                }
                            } else {
                                ui.add_sized([cell_size, cell_size], egui::Label::new(""));
                            }
                        }
                    });
                    if r == 6 || r == 7 {
                        ui.add_space(20.0);
                    }
                }
                ui.add_space(8.0);
                // Legend
                ui.horizontal(|ui| {
                    for (label, color) in [
                        ("Alkali", egui::Color32::from_rgb(0x8B, 0x3A, 0x3A)),
                        ("Alk.earth", egui::Color32::from_rgb(0x8B, 0x69, 0x3A)),
                        ("Transition", egui::Color32::from_rgb(0x3A, 0x5A, 0x7A)),
                        ("Halogen", egui::Color32::from_rgb(0x3A, 0x6B, 0x3A)),
                        ("Noble gas", egui::Color32::from_rgb(0x5A, 0x3A, 0x7A)),
                        ("Non-metal", egui::Color32::from_rgb(0x3A, 0x6B, 0x6B)),
                        ("Lanthanide", egui::Color32::from_rgb(0x7A, 0x4A, 0x5A)),
                        ("Actinide", egui::Color32::from_rgb(0x6A, 0x3A, 0x4A)),
                    ] {
                        ui.add(egui::Button::new(
                            egui::RichText::new(label).size(9.0).color(egui::Color32::WHITE)
                        ).fill(color).min_size(egui::vec2(60.0, 18.0)));
                    }
                });
                ui.add_space(6.0);
                if ui.button("Cancel").clicked() {
                    add_state.show_table = false;
                    add_state.selected_element = None;
                }
            });
    }

    // ── Coordinate input popup ──
    if !add_state.show_table && add_state.selected_element.is_some() {
        let el = add_state.selected_element.as_ref().unwrap().clone();
        egui::Window::new(format!("Add {} atom — fractional coordinates", el))
            .collapsible(false).resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.label("Fractional coordinates:");
                ui.horizontal(|ui| {
                    ui.label("x:"); ui.text_edit_singleline(&mut add_state.coord_x);
                });
                ui.horizontal(|ui| {
                    ui.label("y:"); ui.text_edit_singleline(&mut add_state.coord_y);
                });
                ui.horizontal(|ui| {
                    ui.label("z:"); ui.text_edit_singleline(&mut add_state.coord_z);
                });
                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Confirm").clicked() {
                        if let (Ok(x), Ok(y), Ok(z)) = (
                            add_state.coord_x.parse::<f32>(),
                            add_state.coord_y.parse::<f32>(),
                            add_state.coord_z.parse::<f32>(),
                        ) {
                            add_state.pending = Some((el, x, y, z));
                            add_state.selected_element = None;
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        add_state.selected_element = None;
                    }
                });
            });

        // ── Phonon mode info panel ──
        if let Some(ref mut ps) = phonon_state {
            egui::SidePanel::right("phonon_panel")
                .resizable(false).default_width(200.0)
                .show(ctx, |ui| {
                    ui.heading("Phonon Mode");
                    ui.separator();
                    ui.label(format!("Mode: {}", ps.mode_index));
                    ui.label(format!("Freq: {:.2} cm⁻¹", ps.frequency));
                    ui.label(format!("IR: {:.4} (D/A)²/amu", ps.ir_intensity));
                    ui.label(format!("|p_m|: {:.4}", ps.mode_charge_norm));
                    ui.separator();
                    ui.label("Arrow scale:");
                    let mut sf = ps.scale_factor;
                    if ui.add(egui::Slider::new(&mut sf, 0.1..=10.0)).changed() {
                        ps.scale_factor = sf;
                        // Note: scale change requires re-spawning arrows;
                        // For now, the scale applies on next viewer relaunch
                    }
                    let mut show = ps.show_arrows;
                    if ui.checkbox(&mut show, "Show Arrows").changed() {
                        ps.show_arrows = show;
                    }
                });
        }
    }
}

/// Export the mode-8 migration surface as `migration_surface.ply` (same
/// rebuild the viewer renders, so the export matches the screen).
fn export_migration_ply(cube: &Option<Res<crate::CubeResource>>, ps: &crate::Pes3dState) {
    let result = match cube.as_ref() {
        Some(cube) => {
            let lattice = cube.0.to_lattice();
            let centers = crate::sphere_section::detect_cage_centers(&cube.0);
            match crate::sphere_section::migration_surface_mesh(
                &cube.0.field, cube.0.nx, &lattice, &centers,
                ps.mig_e_cap, ps.mig_show_shell,
                ps.color_min, ps.color_max, ps.iso_material)
            {
                Some(mesh) => crate::write_mesh_ply(&mesh, "migration_surface.ply"),
                None => Err("migration surface build returned None".to_string()),
            }
        }
        None => Err("no cube loaded".to_string()),
    };
    match result {
        Ok((nv, nt)) => {
            let abs = std::env::current_dir()
                .map(|d| d.join("migration_surface.ply"))
                .unwrap_or_else(|_| std::path::PathBuf::from("migration_surface.ply"));
            eprintln!("[mode8] exported {} ({} verts, {} tris)", abs.display(), nv, nt);
        }
        Err(e) => eprintln!("[mode8] export failed: {}", e),
    }
}

/// Slab settings content — shown in the "Slab cut" popup window
/// (moved out of the right panel, which now only shows the button).
fn slab_section_ui(ui: &mut egui::Ui, s: &mut SlabState, c: Option<&crate::CrystalStore>) {
    let hkl_ok: Option<(i32, i32, i32)> = match (
        s.h_str.trim().parse::<i32>(),
        s.k_str.trim().parse::<i32>(),
        s.l_str.trim().parse::<i32>(),
    ) {
        (Ok(h), Ok(k), Ok(l)) if !(h == 0 && k == 0 && l == 0) => Some((h, k, l)),
        _ => None,
    };

    ui.horizontal(|ui| {
        ui.label("h");
        ui.add(egui::TextEdit::singleline(&mut s.h_str).desired_width(40.0));
        ui.label("k");
        ui.add(egui::TextEdit::singleline(&mut s.k_str).desired_width(40.0));
        ui.label("l");
        ui.add(egui::TextEdit::singleline(&mut s.l_str).desired_width(40.0));
    });
    // s: fraction of the layer period + absolute Å, last-edited field wins
    ui.horizontal(|ui| {
        ui.label("s frac");
        let r = ui.add(
            egui::TextEdit::singleline(&mut s.s_frac_str)
                .hint_text("0-1")
                .desired_width(56.0),
        );
        if r.changed() {
            s.s_frac_dirty = true;
        }
        ui.label("s \u{c5}");
        let r = ui.add(
            egui::TextEdit::singleline(&mut s.start_str)
                .hint_text("0")
                .desired_width(64.0),
        );
        if r.changed() {
            s.start_dirty = true;
        }
        if let (Some(c), Some((h, k, l))) = (c, hkl_ok) {
            let s_cur = s.start_str.trim().parse::<f32>().unwrap_or(0.0);
            if ui.add_enabled(
                crate::slab::slab_layer_snap(&c.data, h, k, l, s_cur, true).is_some(),
                egui::Button::new("up"),
            ).on_hover_text("Snap s to the nearest atomic plane above").clicked() {
                if let Some(snap) = crate::slab::slab_layer_snap(&c.data, h, k, l, s_cur, true) {
                    s.start_str = format!("{:.6}", snap);
                    if let Some(per) = crate::slab::slab_period(&c.data, h, k, l) {
                        s.s_frac_str = format!("{:.4}", snap / per);
                    }
                    s.start_dirty = false;
                    s.s_frac_dirty = false;
                }
            }
            if ui.add_enabled(
                crate::slab::slab_layer_snap(&c.data, h, k, l, s_cur, false).is_some(),
                egui::Button::new("down"),
            ).on_hover_text("Snap s to the nearest atomic plane below").clicked() {
                if let Some(snap) = crate::slab::slab_layer_snap(&c.data, h, k, l, s_cur, false) {
                    s.start_str = format!("{:.6}", snap);
                    if let Some(per) = crate::slab::slab_period(&c.data, h, k, l) {
                        s.s_frac_str = format!("{:.4}", snap / per);
                    }
                    s.start_dirty = false;
                    s.s_frac_dirty = false;
                }
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("T frac");
        let r = ui.add(
            egui::TextEdit::singleline(&mut s.t_frac_str)
                .hint_text("3")
                .desired_width(56.0),
        );
        if r.changed() {
            s.t_frac_dirty = true;
        }
        ui.label("T \u{c5}");
        let r = ui.add(
            egui::DragValue::new(&mut s.thickness)
                .speed(0.5)
                .range(0.0..=100.0),
        );
        if r.changed() {
            s.thick_dirty = true;
        }
        if s.thickness == 0.0 {
            ui.label("(0 = 3\u{00b7}period)");
        }
    });
    // Live d_hkl / layer period
    if let (Some(c), Some((h, k, l))) = (c, hkl_ok) {
        if let (Some(d), Some(per)) = (
            crate::slab::slab_d_hkl(&c.data, h, k, l),
            crate::slab::slab_period(&c.data, h, k, l),
        ) {
            ui.label(format!(
                "d = {:.3} \u{c5}   period = {:.3} \u{c5}",
                d, per
            ));
        }
    }
    // U/V definition mode: 0 = integer (in-plane basis × U/V expansion),
    // 1 = explicit (i j k) vectors. Only the selected form's fields show.
    ui.horizontal(|ui| {
        ui.label("U/V def");
        egui::ComboBox::from_id_salt("slab_uv_mode")
            .selected_text(match s.uv_mode {
                0 => "integer \u{00d7} basis",
                _ => "vector (i j k)",
            })
            .show_ui(ui, |cui| {
                let mut cur = s.uv_mode;
                if cui.radio_value(&mut cur, 0, "integer \u{00d7} basis").clicked() {
                    s.uv_mode = 0;
                }
                if cui.radio_value(&mut cur, 1, "vector (i j k)").clicked() {
                    s.uv_mode = 1;
                }
            });
    });
    if s.uv_mode == 0 {
        // ── Integer form: in-plane 2D basis + U/V expansion factors ──
        ui.horizontal(|ui| {
            ui.label("In-plane");
            egui::ComboBox::from_id_salt("slab_basis")
                .selected_text(match s.basis {
                    crate::InPlaneBasis::Primitive => "primitive L2",
                    crate::InPlaneBasis::Orthogonal => "90\u{b0} orthogonal",
                    crate::InPlaneBasis::Conventional => "conventional",
                })
                .show_ui(ui, |cui| {
                    let mut cur = s.basis.as_u8();
                    for i in 0..3u8 {
                        let name = ["primitive", "90\u{b0} orthogonal", "conventional"][i as usize];
                        if cui.radio_value(&mut cur, i, name).clicked() {
                            s.basis = crate::InPlaneBasis::from_u8(i);
                        }
                    }
                });
        });
        ui.horizontal(|ui| {
            ui.label("U");
            ui.add(egui::DragValue::new(&mut s.u).range(1..=8).speed(1));
            ui.label("V");
            ui.add(egui::DragValue::new(&mut s.v).range(1..=8).speed(1));
            if s.basis == crate::InPlaneBasis::Orthogonal {
                ui.label("90\u{b0} #");
                ui.add(egui::DragValue::new(&mut s.orth_idx).range(0..=16).speed(1));
            }
        });
    } else {
        // ── Vector form: MS-style explicit (i j k) in-plane vectors ──
        ui.horizontal(|ui| {
            ui.label("U vec (i j k)");
            ui.add(
                egui::TextEdit::singleline(&mut s.u_vec_str)
                    .hint_text("0 1 0")
                    .desired_width(60.0),
            );
            ui.label("V vec (i j k)");
            ui.add(
                egui::TextEdit::singleline(&mut s.v_vec_str)
                    .hint_text("0 0 1")
                    .desired_width(60.0),
            );
        });
        if s.u_vec_str.trim().is_empty() && s.v_vec_str.trim().is_empty() {
            ui.label(egui::RichText::new("(empty = fall back to integer \u{00d7} basis)")
                .weak()
                .size(11.0));
        }
    }
    let vec_active = s.uv_mode == 1
        && (!s.u_vec_str.trim().is_empty() || !s.v_vec_str.trim().is_empty());
    // Live in-plane cell info
    if let (Some(c), Some((h, k, l))) = (c, hkl_ok) {
        let _ = (h, k, l);
        if let Some(params) = crate::slab_params_from_state(s) {
            if let Some(info) = crate::slab::slab_layer_info(&c.data, &params) {
                let explicit = vec_active;
                let bname = if explicit {
                    "U/V vec"
                } else {
                    match s.basis {
                        crate::InPlaneBasis::Primitive => "prim",
                        crate::InPlaneBasis::Orthogonal => "90\u{b0}",
                        crate::InPlaneBasis::Conventional => "conv",
                    }
                };
                let (a_show, b_show) = if explicit {
                    // slab_layer_info already reports the explicit cell
                    (info.a_prim, info.b_prim)
                } else {
                    (info.a_prim * s.u as f32, info.b_prim * s.v as f32)
                };
                ui.label(format!(
                    "in-plane [{}]: {:.3} \u{00d7} {:.3} \u{c5}, \u{b3} = {:.1}\u{b0}",
                    bname, a_show, b_show, info.gamma_prim
                ));
                let cell_note = if explicit {
                    String::from("explicit cell")
                } else {
                    format!("prim cell \u{00d7} {}x{} basis", s.u, s.v)
                };
                ui.label(format!(
                    "{} atom(s)/layer ({})",
                    info.atoms_per_layer, cell_note
                ));
                if s.uv_mode == 0 && s.basis == crate::InPlaneBasis::Orthogonal {
                    let n = info.ortho.len();
                    ui.label(format!("90\u{b0} candidates: {}", n));
                    if n > 0 {
                        let (oa, ob, _) = info.ortho[(s.orth_idx as usize).min(n - 1)];
                        ui.label(format!("  [cur {:.2}\u{00d7}{:.2}]", oa, ob));
                    }
                }
            }
        }
    }
    ui.horizontal(|ui| {
        if ui.button("Preview").clicked() {
            s.preview_slab = !s.preview_slab;
            s.preview_vacuum = false;
        }
        if ui.button("Apply slab").clicked() {
            s.req_slab = true;
        }
    });
    if !s.error.is_empty() {
        ui.label(egui::RichText::new(&s.error)
            .color(egui::Color32::from_rgb(255, 90, 90)));
    }
}

/// Vacuum settings content — shown in the "Vacuum" popup window.
fn vacuum_section_ui(ui: &mut egui::Ui, s: &mut SlabState) {
    ui.horizontal(|ui| {
        ui.label("Thick V (\u{c5})");
        let r = ui.add(
            egui::DragValue::new(&mut s.vac_thickness)
                .speed(0.5)
                .range(0.0..=300.0),
        );
        if r.changed() {
            s.vac_v_dirty = true;
        }
    });
    ui.horizontal(|ui| {
        ui.label("Cell c (\u{c5})");
        // Text input: empty = auto (c = base + V, shown via the "Auto" hint);
        // type a number to pin c (V then follows); type "Auto" or clear to
        // go back to auto mode.
        let r = ui.add(
            egui::TextEdit::singleline(&mut s.vac_c_str)
                .hint_text("Auto")
                .desired_width(72.0),
        );
        if r.changed() {
            let t = s.vac_c_str.trim().to_lowercase();
            if t.is_empty() || t == "auto" || t == "0" {
                s.vac_cell = 0.0; // auto: c follows V
            } else if let Ok(v) = t.parse::<f32>() {
                s.vac_cell = v.max(0.0);
            }
            s.vac_c_dirty = true;
        }
        if s.vac_cell <= 0.0 {
            ui.label("(Auto = base + V)");
        }
    });
    ui.horizontal(|ui| {
        ui.label("Dir");
        let mut ax = (s.vac_axis - 1).clamp(0, 2) as usize;
        for i in 0..3 {
            let name = ["a", "b", "c"][i];
            if ui.radio_value(&mut ax, i, name).clicked() {
                s.vac_axis = (i + 1) as u8;
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Pos 0\u{2013}1");
        ui.add(egui::Slider::new(&mut s.vac_pos, 0.0..=1.0));
    });
    ui.horizontal(|ui| {
        if ui.button("bottom").clicked() {
            s.vac_pos = 0.0;
            s.vac_both = false;
        }
        if ui.button("center").clicked() {
            s.vac_pos = 0.5;
            s.vac_both = true;
        }
        if ui.button("top").clicked() {
            s.vac_pos = 1.0;
            s.vac_both = false;
        }
        if s.vac_pos > 0.001 {
            ui.label(format!("pos = {:.2}", s.vac_pos));
        }
    });
    ui.horizontal(|ui| {
        if ui.button("Preview").clicked() {
            s.preview_vacuum = !s.preview_vacuum;
            s.preview_slab = false;
        }
        if ui.button("Apply vacuum").clicked() {
            s.req_vacuum = true;
        }
        if ui.add_enabled(
            s.snapshot.is_some(),
            egui::Button::new("Reset"),
        )
        .clicked()
        {
            s.req_reset = true;
        }
    });
    if !s.error.is_empty() {
        ui.label(egui::RichText::new(&s.error)
            .color(egui::Color32::from_rgb(255, 90, 90)));
    }
}

/// Supercell settings content — shown in the "Supercell" popup window.
/// x/y/z are integer multipliers along a/b/c (1,1,1 = the current cell,
/// 2,1,1 = one extra cell along a). Preview keeps every ORIGINAL cell's
/// own edges (boxes side by side) plus all atom copies so the tiling is
/// visible; only "Apply" merges the boxes into one cube.
fn supercell_section_ui(ui: &mut egui::Ui, s: &mut SlabState, c: Option<&crate::CrystalStore>) {
    ui.horizontal(|ui| {
        ui.label("x");
        ui.add(egui::DragValue::new(&mut s.sc_x).range(1..=8));
        ui.label("y");
        ui.add(egui::DragValue::new(&mut s.sc_y).range(1..=8));
        ui.label("z");
        ui.add(egui::DragValue::new(&mut s.sc_z).range(1..=8));
    });
    let n = s.sc_x.max(1) * s.sc_y.max(1) * s.sc_z.max(1);
    if let Some(c) = c {
        let atoms = c.data.atoms.len();
        ui.label(format!(
            "{}\u{00d7}{}\u{00d7}{} = {} cells \u{00d7} {} atoms = {} atoms",
            s.sc_x.max(1), s.sc_y.max(1), s.sc_z.max(1), n, atoms, n as usize * atoms
        ));
        if let Some(sc) = c.data.supercell {
            ui.label(format!("Applied: {}\u{00d7}{}\u{00d7}{}", sc[0], sc[1], sc[2]));
        }
        if n > 512 {
            ui.label("preview is capped at 512 cells (Apply still works up to the atom limit)");
        }
    }
    ui.horizontal(|ui| {
        if s.preview_supercell {
            if ui.button("Stop preview").clicked() {
                s.preview_supercell = false;
            }
        } else {
            if ui.button("Preview").clicked() {
                s.preview_supercell = true;
                s.preview_slab = false;
                s.preview_vacuum = false;
            }
        }
        if ui.button("Apply supercell").clicked() {
            s.req_supercell = true;
            s.preview_supercell = false;
        }
    });
    if s.preview_supercell {
        ui.label("Preview: each original cell keeps its own edges until you apply.");
    }
    if !s.error.is_empty() {
        ui.label(egui::RichText::new(&s.error)
            .color(egui::Color32::from_rgb(255, 90, 90)));
    }
}
