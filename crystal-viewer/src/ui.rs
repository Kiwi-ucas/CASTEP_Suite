//! egui UI panels: atom list, info panel, toolbar

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use crate::picking::PickingState;

#[derive(Resource)]
pub struct AtomInfo {
    pub elements: Vec<String>,
    pub labels: Vec<String>,
}

#[derive(Resource)]
pub struct CrystalMeta {
    pub filename: String,
    pub a: f32, pub b: f32, pub c: f32,
    pub alpha: f32, pub beta: f32, pub gamma: f32,
}

impl AtomInfo {
    pub fn new(elements: Vec<String>, labels: Vec<String>) -> Self {
        Self { elements, labels }
    }
}

pub fn ui_system(
    mut contexts: EguiContexts,
    mut picking: ResMut<PickingState>,
    atom_info: Res<AtomInfo>,
    meta: Option<Res<CrystalMeta>>,
) {
    let ctx = contexts.ctx_mut();

    // ── Bottom toolbar ──
    egui::TopBottomPanel::bottom("toolbar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            ui.label("\u{1f5b0} Right-drag: rotate | Scroll: zoom | Click: select | 1/2/3: mode | P: proj | B: bonds | C: cell | R: reset");
        });
    });

    // ── Left panel: atom list ──
    egui::SidePanel::left("atom_list")
        .resizable(true).default_width(180.0)
        .show(ctx, |ui| {
            ui.heading("Atoms");
            ui.separator();
            let count = atom_info.elements.len();
            egui::ScrollArea::vertical().show(ui, |ui| {
                for i in 0..count {
                    let label = format!("{:2}. {:3}", i, atom_info.elements[i]);
                    let selected = picking.selected == i as i32;
                    let hovered = picking.hovered == i as i32;

                    let mut text = egui::RichText::new(&label);
                    if selected {
                        text = text.color(egui::Color32::from_rgb(255, 170, 30)).strong();
                    } else if hovered {
                        text = text.color(egui::Color32::from_rgb(180, 180, 180));
                    }
                    let response = ui.selectable_label(selected, text);
                    if response.clicked() { picking.selected = i as i32; }
                    if response.hovered() { picking.hovered = i as i32; }
                }
            });
        });

    // ── Right panel ──
    egui::SidePanel::right("info_panel")
        .resizable(true).default_width(220.0)
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
                ui.label(format!("Atoms: {}", atom_info.elements.len()));
                ui.separator();
            } else {
                ui.heading("Crystal Viewer");
            }

            // Selected atom details
            if picking.selected >= 0 && picking.selected < atom_info.elements.len() as i32 {
                let i = picking.selected as usize;
                ui.label(egui::RichText::new("Selected Atom").strong());
                ui.label(format!("Element: {}", atom_info.elements[i]));
                ui.label(format!("Index: {}", i));
                if i < picking.atom_positions.len() {
                    let p = picking.atom_positions[i];
                    ui.label(format!("X: {:.4}", p.x));
                    ui.label(format!("Y: {:.4}", p.y));
                    ui.label(format!("Z: {:.4}", p.z));
                }
                if ui.button("Deselect").clicked() { picking.selected = -1; }
            } else {
                ui.label("Click an atom\nto see details");
            }
        });
}
