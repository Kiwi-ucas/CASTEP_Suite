//! Lookup tables for element properties (radius, color)

use bevy::prelude::Color;
use bevy_egui::egui;

/// Covalent radii in Angstrom (from Cordero et al., 2008)
pub fn covalent_radius(element: &str) -> f32 {
    match element {
        "H" => 0.31, "He" => 0.28, "Li" => 1.28, "Be" => 0.96, "B" => 0.84,
        "C" => 0.76, "N" => 0.71, "O" => 0.66, "F" => 0.57, "Ne" => 0.58,
        "Na" => 1.66, "Mg" => 1.41, "Al" => 1.21, "Si" => 1.11, "P" => 1.07,
        "S" => 1.05, "Cl" => 0.99, "Ar" => 1.06, "K" => 2.03, "Ca" => 1.76,
        "Sc" => 1.70, "Ti" => 1.60, "V" => 1.53, "Cr" => 1.39, "Mn" => 1.39,
        "Fe" => 1.32, "Co" => 1.26, "Ni" => 1.24, "Cu" => 1.28, "Zn" => 1.22,
        "Ga" => 1.22, "Ge" => 1.20, "As" => 1.19, "Se" => 1.20, "Br" => 1.20,
        "Kr" => 1.16, "Rb" => 2.20, "Sr" => 1.95, "Y" => 1.90, "Zr" => 1.75,
        "Nb" => 1.64, "Mo" => 1.54, "Tc" => 1.47, "Ru" => 1.46, "Rh" => 1.42,
        "Pd" => 1.39, "Ag" => 1.45, "Cd" => 1.44, "In" => 1.42, "Sn" => 1.39,
        "Sb" => 1.39, "Te" => 1.38, "I" => 1.39, "Xe" => 1.40,
        "Cs" => 2.44, "Ba" => 2.15, "La" => 2.07, "Ce" => 2.04, "Pr" => 2.03,
        "Nd" => 2.01, "Pm" => 1.99, "Sm" => 1.98, "Eu" => 1.98, "Gd" => 1.96,
        "Tb" => 1.94, "Dy" => 1.92, "Ho" => 1.92, "Er" => 1.89, "Tm" => 1.90,
        "Yb" => 1.87, "Lu" => 1.87, "Hf" => 1.75, "Ta" => 1.70, "W" => 1.62,
        "Re" => 1.51, "Os" => 1.44, "Ir" => 1.41, "Pt" => 1.36, "Au" => 1.36,
        "Hg" => 1.32, "Tl" => 1.45, "Pb" => 1.46, "Bi" => 1.48, "Po" => 1.40,
        "At" => 1.50, "Rn" => 1.50, "Fr" => 2.60, "Ra" => 2.21, "Ac" => 2.15,
        "Th" => 2.06, "Pa" => 2.00, "U" => 1.96, "Np" => 1.90, "Pu" => 1.87,
        "Am" => 1.80, "Cm" => 1.69,
        _ => 1.50,
    }
}

/// CPK / Jmol colors
pub fn element_color(element: &str) -> Color {
    match element {
        "H"  => Color::srgb(1.00, 1.00, 1.00),
        "He" => Color::srgb(0.85, 1.00, 1.00),
        "Li" => Color::srgb(0.80, 0.50, 1.00),
        "Be" => Color::srgb(0.76, 1.00, 0.00),
        "B"  => Color::srgb(1.00, 0.71, 0.71),
        "C"  => Color::srgb(0.20, 0.20, 0.20),
        "N"  => Color::srgb(0.14, 0.14, 1.00),
        "O"  => Color::srgb(1.00, 0.05, 0.05),
        "F"  => Color::srgb(0.56, 0.88, 0.31),
        "Ne" => Color::srgb(0.70, 1.00, 1.00),
        "Na" => Color::srgb(0.67, 0.36, 0.95),
        "Mg" => Color::srgb(0.54, 1.00, 0.00),
        "Al" => Color::srgb(0.75, 0.65, 0.65),
        "Si" => Color::srgb(0.94, 0.78, 0.63),
        "P"  => Color::srgb(1.00, 0.50, 0.00),
        "S"  => Color::srgb(1.00, 1.00, 0.00),
        "Cl" => Color::srgb(0.12, 0.94, 0.12),
        "Ar" => Color::srgb(0.50, 1.00, 1.00),
        "K"  => Color::srgb(0.56, 0.25, 0.83),
        "Ca" => Color::srgb(0.24, 1.00, 0.00),
        "Sc" => Color::srgb(0.90, 0.90, 0.90),
        "Ti" => Color::srgb(0.75, 0.76, 0.78),
        "V"  => Color::srgb(0.65, 0.65, 0.67),
        "Cr" => Color::srgb(0.54, 0.60, 0.78),
        "Mn" => Color::srgb(0.61, 0.48, 0.78),
        "Fe" => Color::srgb(0.88, 0.40, 0.20),
        "Co" => Color::srgb(0.94, 0.56, 0.63),
        "Ni" => Color::srgb(0.31, 0.82, 0.31),
        "Cu" => Color::srgb(0.72, 0.45, 0.20),
        "Zn" => Color::srgb(0.49, 0.50, 0.69),
        "Ga" => Color::srgb(0.76, 0.56, 0.56),
        "Ge" => Color::srgb(0.40, 0.56, 0.56),
        "As" => Color::srgb(0.74, 0.50, 0.89),
        "Se" => Color::srgb(1.00, 0.63, 0.00),
        "Br" => Color::srgb(0.65, 0.16, 0.16),
        "Kr" => Color::srgb(0.36, 0.72, 0.82),
        "Rb" => Color::srgb(0.44, 0.18, 0.69),
        "Sr" => Color::srgb(0.00, 1.00, 0.00),
        "Y"  => Color::srgb(0.58, 1.00, 1.00),
        "Zr" => Color::srgb(0.58, 0.88, 0.88),
        "Nb" => Color::srgb(0.45, 0.76, 0.79),
        "Mo" => Color::srgb(0.33, 0.71, 0.71),
        "Tc" => Color::srgb(0.23, 0.62, 0.62),
        "Ru" => Color::srgb(0.14, 0.56, 0.56),
        "Rh" => Color::srgb(0.04, 0.49, 0.55),
        "Pd" => Color::srgb(0.00, 0.41, 0.52),
        "Ag" => Color::srgb(0.88, 0.88, 1.00),
        "Cd" => Color::srgb(1.00, 0.85, 0.56),
        "In" => Color::srgb(0.65, 0.46, 0.45),
        "Sn" => Color::srgb(0.40, 0.50, 0.50),
        "Sb" => Color::srgb(0.62, 0.39, 0.71),
        "Te" => Color::srgb(0.83, 0.48, 0.00),
        "I"  => Color::srgb(0.58, 0.00, 0.58),
        "Xe" => Color::srgb(0.26, 0.62, 0.69),
        "Cs" => Color::srgb(0.34, 0.09, 0.56),
        "Ba" => Color::srgb(0.00, 0.79, 0.00),
        "La" => Color::srgb(0.44, 0.83, 1.00),
        "Ce" => Color::srgb(1.00, 1.00, 0.78),
        "Pr" => Color::srgb(0.85, 1.00, 0.78),
        "Nd" => Color::srgb(0.78, 1.00, 0.78),
        "Pm" => Color::srgb(0.64, 1.00, 0.78),
        "Sm" => Color::srgb(0.56, 1.00, 0.78),
        "Eu" => Color::srgb(0.38, 1.00, 0.78),
        "Gd" => Color::srgb(0.27, 1.00, 0.78),
        "Tb" => Color::srgb(0.19, 1.00, 0.78),
        "Dy" => Color::srgb(0.12, 1.00, 0.78),
        "Ho" => Color::srgb(0.00, 1.00, 0.61),
        "Er" => Color::srgb(0.00, 0.90, 0.46),
        "Tm" => Color::srgb(0.00, 0.83, 0.32),
        "Yb" => Color::srgb(0.00, 0.75, 0.22),
        "Lu" => Color::srgb(0.00, 0.67, 0.14),
        "Hf" => Color::srgb(0.30, 0.76, 1.00),
        "Ta" => Color::srgb(0.30, 0.65, 1.00),
        "W"  => Color::srgb(0.13, 0.58, 0.84),
        "Re" => Color::srgb(0.15, 0.49, 0.67),
        "Os" => Color::srgb(0.15, 0.40, 0.59),
        "Ir" => Color::srgb(0.09, 0.33, 0.53),
        "Pt" => Color::srgb(0.82, 0.82, 0.88),
        "Au" => Color::srgb(1.00, 0.84, 0.14),
        "Hg" => Color::srgb(0.72, 0.72, 0.82),
        "Tl" => Color::srgb(0.65, 0.33, 0.30),
        "Pb" => Color::srgb(0.34, 0.35, 0.38),
        "Bi" => Color::srgb(0.62, 0.31, 0.71),
        "Po" => Color::srgb(0.67, 0.36, 0.00),
        "At" => Color::srgb(0.46, 0.31, 0.27),
        "Rn" => Color::srgb(0.26, 0.51, 0.59),
        "Fr" => Color::srgb(0.26, 0.00, 0.40),
        "Ra" => Color::srgb(0.00, 0.49, 0.00),
        "Ac" => Color::srgb(0.44, 0.67, 0.98),
        "Th" => Color::srgb(0.00, 0.73, 1.00),
        "Pa" => Color::srgb(0.00, 0.63, 1.00),
        "U"  => Color::srgb(0.00, 0.56, 1.00),
        "Np" => Color::srgb(0.00, 0.50, 1.00),
        "Pu" => Color::srgb(0.00, 0.42, 1.00),
        "Am" => Color::srgb(0.33, 0.36, 0.95),
        "Cm" => Color::srgb(0.47, 0.36, 0.89),
        _ => Color::srgb(1.00, 0.08, 0.58), // pink for unknown
    }
}

/// Atom display radius (scaled for visibility, not physical)
/// Category background color for periodic table (muted, white text readable)
pub fn category_color(element: &str) -> egui::Color32 {
    match element {
        // Alkali metals — muted red
        "H" | "Li" | "Na" | "K" | "Rb" | "Cs" | "Fr" =>
            egui::Color32::from_rgb(0x8B, 0x3A, 0x3A),
        // Alkaline earth — muted orange
        "Be" | "Mg" | "Ca" | "Sr" | "Ba" | "Ra" =>
            egui::Color32::from_rgb(0x8B, 0x69, 0x3A),
        // Transition metals — steel blue
        "Sc" | "Ti" | "V" | "Cr" | "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Zn" |
        "Y" | "Zr" | "Nb" | "Mo" | "Tc" | "Ru" | "Rh" | "Pd" | "Ag" | "Cd" |
        "Hf" | "Ta" | "W" | "Re" | "Os" | "Ir" | "Pt" | "Au" | "Hg" |
        "Rf" | "Db" | "Sg" | "Bh" | "Hs" | "Mt" | "Ds" | "Rg" | "Cn" =>
            egui::Color32::from_rgb(0x3A, 0x5A, 0x7A),
        // Halogens — muted green
        "F" | "Cl" | "Br" | "I" | "At" | "Ts" =>
            egui::Color32::from_rgb(0x3A, 0x6B, 0x3A),
        // Noble gases — muted purple
        "He" | "Ne" | "Ar" | "Kr" | "Xe" | "Rn" | "Og" =>
            egui::Color32::from_rgb(0x5A, 0x3A, 0x7A),
        // Non-metals — muted teal
        "C" | "N" | "O" | "P" | "S" | "Se" =>
            egui::Color32::from_rgb(0x3A, 0x6B, 0x6B),
        // Metalloids — gray-blue
        "B" | "Si" | "Ge" | "As" | "Sb" | "Te" | "Po" =>
            egui::Color32::from_rgb(0x4A, 0x5A, 0x6A),
        // Post-transition metals — gray
        "Al" | "Ga" | "In" | "Sn" | "Tl" | "Pb" | "Bi" | "Nh" | "Fl" | "Mc" | "Lv" =>
            egui::Color32::from_rgb(0x5A, 0x5A, 0x5A),
        // Lanthanides — muted pink
        "La" | "Ce" | "Pr" | "Nd" | "Pm" | "Sm" | "Eu" |
        "Gd" | "Tb" | "Dy" | "Ho" | "Er" | "Tm" | "Yb" | "Lu" =>
            egui::Color32::from_rgb(0x7A, 0x4A, 0x5A),
        // Actinides — darker pink
        "Ac" | "Th" | "Pa" | "U" | "Np" | "Pu" | "Am" |
        "Cm" | "Bk" | "Cf" | "Es" | "Fm" | "Md" | "No" | "Lr" =>
            egui::Color32::from_rgb(0x6A, 0x3A, 0x4A),
        _ => egui::Color32::from_rgb(0x50, 0x50, 0x50),
    }
}

/// Periodic table layout: (symbol, row, col) — 9 rows × 18 columns
pub const PERIODIC_TABLE: &[(&str, usize, usize)] = &[
    ("H", 0, 0),  ("He", 0, 17),
    ("Li", 1, 0), ("Be", 1, 1), ("B", 1, 12), ("C", 1, 13), ("N", 1, 14), ("O", 1, 15), ("F", 1, 16), ("Ne", 1, 17),
    ("Na", 2, 0), ("Mg", 2, 1), ("Al", 2, 12), ("Si", 2, 13), ("P", 2, 14), ("S", 2, 15), ("Cl", 2, 16), ("Ar", 2, 17),
    ("K", 3, 0),  ("Ca", 3, 1), ("Sc", 3, 2), ("Ti", 3, 3), ("V", 3, 4), ("Cr", 3, 5), ("Mn", 3, 6),
    ("Fe", 3, 7), ("Co", 3, 8), ("Ni", 3, 9), ("Cu", 3, 10), ("Zn", 3, 11), ("Ga", 3, 12), ("Ge", 3, 13),
    ("As", 3, 14), ("Se", 3, 15), ("Br", 3, 16), ("Kr", 3, 17),
    ("Rb", 4, 0),  ("Sr", 4, 1), ("Y", 4, 2), ("Zr", 4, 3), ("Nb", 4, 4), ("Mo", 4, 5), ("Tc", 4, 6),
    ("Ru", 4, 7), ("Rh", 4, 8), ("Pd", 4, 9), ("Ag", 4, 10), ("Cd", 4, 11), ("In", 4, 12), ("Sn", 4, 13),
    ("Sb", 4, 14), ("Te", 4, 15), ("I", 4, 16), ("Xe", 4, 17),
    ("Cs", 5, 0),  ("Ba", 5, 1),
    ("Hf", 5, 3), ("Ta", 5, 4), ("W", 5, 5), ("Re", 5, 6), ("Os", 5, 7), ("Ir", 5, 8), ("Pt", 5, 9),
    ("Au", 5, 10), ("Hg", 5, 11), ("Tl", 5, 12), ("Pb", 5, 13), ("Bi", 5, 14), ("Po", 5, 15), ("At", 5, 16), ("Rn", 5, 17),
    ("Fr", 6, 0),  ("Ra", 6, 1),
    ("Rf", 6, 3), ("Db", 6, 4), ("Sg", 6, 5), ("Bh", 6, 6), ("Hs", 6, 7), ("Mt", 6, 8), ("Ds", 6, 9),
    ("Rg", 6, 10), ("Cn", 6, 11), ("Nh", 6, 12), ("Fl", 6, 13), ("Mc", 6, 14), ("Lv", 6, 15), ("Ts", 6, 16), ("Og", 6, 17),
    // Lanthanides (row 5, col 2+)
    ("La", 7, 2), ("Ce", 7, 3), ("Pr", 7, 4), ("Nd", 7, 5), ("Pm", 7, 6), ("Sm", 7, 7), ("Eu", 7, 8),
    ("Gd", 7, 9), ("Tb", 7, 10), ("Dy", 7, 11), ("Ho", 7, 12), ("Er", 7, 13), ("Tm", 7, 14), ("Yb", 7, 15), ("Lu", 7, 16),
    // Actinides (row 8)
    ("Ac", 8, 2), ("Th", 8, 3), ("Pa", 8, 4), ("U", 8, 5), ("Np", 8, 6), ("Pu", 8, 7), ("Am", 8, 8),
    ("Cm", 8, 9), ("Bk", 8, 10), ("Cf", 8, 11), ("Es", 8, 12), ("Fm", 8, 13), ("Md", 8, 14), ("No", 8, 15), ("Lr", 8, 16),
];

/// Slater empirical atomic radii (Å). Slater, JCP 1964, 41, 3199.
/// Noble gases from Clementi. Missing superheavies fall back to 1.50.
pub fn atomic_radius(element: &str) -> f32 {
    match element {
        "H" => 0.25, "He" => 0.31, "Li" => 1.45, "Be" => 1.05, "B" => 0.85,
        "C" => 0.70, "N" => 0.65, "O" => 0.60, "F" => 0.50, "Ne" => 0.38,
        "Na" => 1.80, "Mg" => 1.50, "Al" => 1.25, "Si" => 1.10, "P" => 1.00,
        "S" => 1.00, "Cl" => 1.00, "Ar" => 0.71, "K" => 2.20, "Ca" => 1.80,
        "Sc" => 1.60, "Ti" => 1.40, "V" => 1.35, "Cr" => 1.40, "Mn" => 1.40,
        "Fe" => 1.40, "Co" => 1.35, "Ni" => 1.35, "Cu" => 1.35, "Zn" => 1.35,
        "Ga" => 1.30, "Ge" => 1.25, "As" => 1.15, "Se" => 1.15, "Br" => 1.15,
        "Kr" => 0.88, "Rb" => 2.35, "Sr" => 2.00, "Y" => 1.80, "Zr" => 1.55,
        "Nb" => 1.45, "Mo" => 1.45, "Tc" => 1.35, "Ru" => 1.30, "Rh" => 1.35,
        "Pd" => 1.40, "Ag" => 1.60, "Cd" => 1.55, "In" => 1.55, "Sn" => 1.45,
        "Sb" => 1.45, "Te" => 1.40, "I" => 1.40, "Xe" => 1.08, "Cs" => 2.60,
        "Ba" => 2.15, "La" => 1.95, "Ce" => 1.85, "Pr" => 1.85, "Nd" => 1.85,
        "Pm" => 1.85, "Sm" => 1.85, "Eu" => 1.85, "Gd" => 1.80, "Tb" => 1.75,
        "Dy" => 1.75, "Ho" => 1.75, "Er" => 1.75, "Tm" => 1.75, "Yb" => 1.75,
        "Lu" => 1.75, "Hf" => 1.55, "Ta" => 1.45, "W" => 1.35, "Re" => 1.35,
        "Os" => 1.30, "Ir" => 1.35, "Pt" => 1.35, "Au" => 1.35, "Hg" => 1.50,
        "Tl" => 1.90, "Pb" => 1.80, "Bi" => 1.60, "Po" => 1.90, "At" => 1.50,
        "Rn" => 1.50, "Fr" => 2.60, "Ra" => 2.15, "Ac" => 1.95, "Th" => 1.80,
        "Pa" => 1.80, "U" => 1.75, "Np" => 1.75, "Pu" => 1.75, "Am" => 1.75,
        _ => 1.50,
    }
}

/// Shannon ionic radii (Å) for CN=6, most common oxidation state.
/// Shannon, Acta Cryst. 1976, A32, 751.
/// Fallback: atomic_radius
pub fn ionic_radius(element: &str) -> f32 {
    match element {
        "Li" => 0.76, "Na" => 1.02, "K" => 1.38, "Rb" => 1.52, "Cs" => 1.67,
        "Be" => 0.45, "Mg" => 0.72, "Ca" => 1.00, "Sr" => 1.18, "Ba" => 1.35,
        "Sc" => 0.745, "Ti" => 0.67, "V" => 0.64, "Cr" => 0.615, "Mn" => 0.645,
        "Fe" => 0.645, "Co" => 0.61, "Ni" => 0.69, "Cu" => 0.73, "Zn" => 0.74,
        "Ga" => 0.62, "Ge" => 0.53, "Al" => 0.535, "Si" => 0.40,
        "Y" => 0.90, "Zr" => 0.72, "Nb" => 0.68, "Mo" => 0.65, "Tc" => 0.645,
        "Ru" => 0.62, "Rh" => 0.665, "Pd" => 0.76, "Ag" => 1.15, "Cd" => 0.95,
        "In" => 0.80, "Sn" => 0.69, "Sb" => 0.60, "Hf" => 0.71, "Ta" => 0.64,
        "W" => 0.60, "Re" => 0.63, "Os" => 0.63, "Ir" => 0.625, "Pt" => 0.625,
        "Au" => 0.85, "Hg" => 1.02, "Tl" => 0.885, "Pb" => 0.775, "Bi" => 0.76,
        // Lanthanides (3+)
        "La" => 1.032, "Ce" => 1.01, "Pr" => 0.99, "Nd" => 0.983, "Pm" => 0.97,
        "Sm" => 0.958, "Eu" => 0.947, "Gd" => 0.938, "Tb" => 0.923,
        "Dy" => 0.912, "Ho" => 0.901, "Er" => 0.890, "Tm" => 0.880,
        "Yb" => 0.868, "Lu" => 0.861,
        // Actinides
        "Ac" => 1.12, "Th" => 0.94, "Pa" => 0.90, "U" => 0.89, "Np" => 0.87,
        "Pu" => 0.86, "Am" => 0.975, "Cm" => 0.97,
        // Common p-block
        "B" => 0.27, "Te" => 0.97, "I" => 2.20, "Xe" => 0.48,
        "As" => 0.58,
        _ => atomic_radius(element),
    }
}

/// Detect bonds between atoms: true if distance < sum of covalent radii * factor
pub fn has_bond(el1: &str, el2: &str, distance: f32, factor: f32) -> bool {
    let r1 = covalent_radius(el1);
    let r2 = covalent_radius(el2);
    distance < (r1 + r2) * factor
}
