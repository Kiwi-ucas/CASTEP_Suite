# CASTEP_Suite Makefile
FC      = gfortran
FCFLAGS = -std=f2018 -g -O2 -ffree-form -fimplicit-none -Wall -Wextra -Wno-unused-dummy-argument
LDFLAGS =

SRCDIR  = src
OBJDIR  = obj
TARGET  = CASTEP_Suite

VIEWER_DIR    = crystal-viewer
VIEWER_TARGET = $(VIEWER_DIR)/target/release/crystal-viewer

.PHONY: all clean run debug fortran viewer

all: $(TARGET)
	@if command -v cargo > /dev/null 2>&1; then \
		$(MAKE) viewer; \
	else \
		echo "  [skip] Rust not found — crystal-viewer not built"; \
	fi

viewer:
	@echo "  Building crystal-viewer..."
	@cargo build --release --manifest-path $(VIEWER_DIR)/Cargo.toml
	@echo "  Cleaning intermediate files..."
	@rm -rf $(VIEWER_DIR)/target/release/deps
	@rm -rf $(VIEWER_DIR)/target/release/build
	@rm -rf $(VIEWER_DIR)/target/release/incremental
	@rm -f $(VIEWER_DIR)/target/release/*.d
	@du -sh $(VIEWER_DIR)/target/release/crystal-viewer

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(TARGET): $(SRCDIR)/config.f90 $(SRCDIR)/term_utils.f90 $(SRCDIR)/symmetry.f90 \
           $(SRCDIR)/parser.f90 \
           $(SRCDIR)/cell_writer.f90 $(SRCDIR)/param_writer.f90 \
           $(SRCDIR)/bands_parser.f90 $(SRCDIR)/bands_plotter.f90 \
           $(SRCDIR)/pdos_parser.f90 $(SRCDIR)/dos_compute.f90 \
           $(SRCDIR)/dos_plotter.f90 $(SRCDIR)/cli_menu.f90 \
           $(SRCDIR)/phonon_dos.f90 $(SRCDIR)/polarizability.f90 \
           $(SRCDIR)/phonon_modes.f90 $(SRCDIR)/crystal_json.f90 \
           $(SRCDIR)/thermodynamics.f90 $(SRCDIR)/castep_vib.f90 \
           $(SRCDIR)/poscastep_menu.f90 $(SRCDIR)/main.f90 \
           | $(OBJDIR)
	@$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $(OBJDIR)/config.o $(SRCDIR)/config.f90
	@$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $(OBJDIR)/term_utils.o $(SRCDIR)/term_utils.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/symmetry.o $(SRCDIR)/symmetry.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/parser.o $(SRCDIR)/parser.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/cell_writer.o $(SRCDIR)/cell_writer.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/param_writer.o $(SRCDIR)/param_writer.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/bands_parser.o $(SRCDIR)/bands_parser.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/bands_plotter.o $(SRCDIR)/bands_plotter.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/pdos_parser.o $(SRCDIR)/pdos_parser.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/dos_compute.o $(SRCDIR)/dos_compute.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/dos_plotter.o $(SRCDIR)/dos_plotter.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/cli_menu.o $(SRCDIR)/cli_menu.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/phonon_dos.o $(SRCDIR)/phonon_dos.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/polarizability.o $(SRCDIR)/polarizability.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/phonon_modes.o $(SRCDIR)/phonon_modes.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/crystal_json.o $(SRCDIR)/crystal_json.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/thermodynamics.o $(SRCDIR)/thermodynamics.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/castep_vib.o $(SRCDIR)/castep_vib.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/poscastep_menu.o $(SRCDIR)/poscastep_menu.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/main.o $(SRCDIR)/main.f90
	@$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $(OBJDIR)/config.o $(OBJDIR)/term_utils.o \
	    $(OBJDIR)/symmetry.o $(OBJDIR)/parser.o $(OBJDIR)/cell_writer.o \
	    $(OBJDIR)/param_writer.o \
	    $(OBJDIR)/bands_parser.o $(OBJDIR)/bands_plotter.o $(OBJDIR)/pdos_parser.o \
	    $(OBJDIR)/dos_compute.o $(OBJDIR)/dos_plotter.o $(OBJDIR)/cli_menu.o \
	    $(OBJDIR)/phonon_dos.o $(OBJDIR)/polarizability.o $(OBJDIR)/phonon_modes.o \
	    $(OBJDIR)/crystal_json.o $(OBJDIR)/thermodynamics.o $(OBJDIR)/castep_vib.o \
	    $(OBJDIR)/poscastep_menu.o $(OBJDIR)/main.o

fortran: $(SRCDIR)/config.f90 $(SRCDIR)/term_utils.f90 $(SRCDIR)/symmetry.f90 \
	           $(SRCDIR)/parser.f90 \
	           $(SRCDIR)/cell_writer.f90 $(SRCDIR)/param_writer.f90 \
	           $(SRCDIR)/bands_parser.f90 $(SRCDIR)/bands_plotter.f90 \
	           $(SRCDIR)/pdos_parser.f90 $(SRCDIR)/dos_compute.f90 \
	           $(SRCDIR)/dos_plotter.f90 $(SRCDIR)/cli_menu.f90 \
	           $(SRCDIR)/phonon_dos.f90 $(SRCDIR)/polarizability.f90 \
	           $(SRCDIR)/phonon_modes.f90 $(SRCDIR)/crystal_json.f90 \
	           $(SRCDIR)/thermodynamics.f90 $(SRCDIR)/castep_vib.f90 \
	           $(SRCDIR)/poscastep_menu.f90 $(SRCDIR)/main.f90 \
	           | $(OBJDIR)
	@$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $(OBJDIR)/config.o $(SRCDIR)/config.f90
	@$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $(OBJDIR)/term_utils.o $(SRCDIR)/term_utils.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/symmetry.o $(SRCDIR)/symmetry.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/parser.o $(SRCDIR)/parser.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/cell_writer.o $(SRCDIR)/cell_writer.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/param_writer.o $(SRCDIR)/param_writer.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/bands_parser.o $(SRCDIR)/bands_parser.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/bands_plotter.o $(SRCDIR)/bands_plotter.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/pdos_parser.o $(SRCDIR)/pdos_parser.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/dos_compute.o $(SRCDIR)/dos_compute.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/dos_plotter.o $(SRCDIR)/dos_plotter.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/cli_menu.o $(SRCDIR)/cli_menu.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/phonon_dos.o $(SRCDIR)/phonon_dos.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/polarizability.o $(SRCDIR)/polarizability.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/phonon_modes.o $(SRCDIR)/phonon_modes.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/crystal_json.o $(SRCDIR)/crystal_json.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/thermodynamics.o $(SRCDIR)/thermodynamics.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/castep_vib.o $(SRCDIR)/castep_vib.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/poscastep_menu.o $(SRCDIR)/poscastep_menu.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/main.o $(SRCDIR)/main.f90
	@$(FC) $(FCFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJDIR)/config.o $(OBJDIR)/term_utils.o \
	    $(OBJDIR)/symmetry.o $(OBJDIR)/parser.o $(OBJDIR)/cell_writer.o \
	    $(OBJDIR)/param_writer.o \
	    $(OBJDIR)/bands_parser.o $(OBJDIR)/bands_plotter.o $(OBJDIR)/pdos_parser.o \
	    $(OBJDIR)/dos_compute.o $(OBJDIR)/dos_plotter.o $(OBJDIR)/cli_menu.o \
	    $(OBJDIR)/phonon_dos.o $(OBJDIR)/polarizability.o $(OBJDIR)/phonon_modes.o \
	    $(OBJDIR)/crystal_json.o $(OBJDIR)/thermodynamics.o $(OBJDIR)/castep_vib.o \
	    $(OBJDIR)/poscastep_menu.o $(OBJDIR)/main.o

clean:
	rm -rf $(OBJDIR) $(TARGET)
	rm -rf $(VIEWER_DIR)/target

run: $(TARGET)
	./$(TARGET)

debug: FCFLAGS += -O0 -fcheck=all -g -fbacktrace
debug: $(TARGET)
