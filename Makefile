# CASTEP_Suite Makefile
FC      = gfortran
FCFLAGS = -std=f2008 -g -O2 -ffree-form -fimplicit-none -Wall -Wextra -Wno-unused-dummy-argument
LDFLAGS =

SRCDIR  = src
OBJDIR  = obj
TARGET  = CASTEP_Suite

.PHONY: all clean run debug

all: $(TARGET)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(TARGET): $(SRCDIR)/config.f90 $(SRCDIR)/term_utils.f90 $(SRCDIR)/parser.f90 \
           $(SRCDIR)/cell_writer.f90 $(SRCDIR)/param_writer.f90 \
           $(SRCDIR)/bands_parser.f90 $(SRCDIR)/bands_plotter.f90 \
           $(SRCDIR)/pdos_parser.f90 $(SRCDIR)/dos_compute.f90 \
           $(SRCDIR)/dos_plotter.f90 $(SRCDIR)/cli_menu.f90 \
           $(SRCDIR)/phonon_dos.f90 $(SRCDIR)/polarizability.f90 \
           $(SRCDIR)/poscastep_menu.f90 $(SRCDIR)/main.f90 \
           | $(OBJDIR)
	@$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $(OBJDIR)/config.o $(SRCDIR)/config.f90
	@$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $(OBJDIR)/term_utils.o $(SRCDIR)/term_utils.f90
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
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/poscastep_menu.o $(SRCDIR)/poscastep_menu.f90
	@$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $(OBJDIR)/main.o $(SRCDIR)/main.f90
	@$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $(OBJDIR)/config.o $(OBJDIR)/term_utils.o \
	    $(OBJDIR)/parser.o $(OBJDIR)/cell_writer.o $(OBJDIR)/param_writer.o \
	    $(OBJDIR)/bands_parser.o $(OBJDIR)/bands_plotter.o $(OBJDIR)/pdos_parser.o \
	    $(OBJDIR)/dos_compute.o $(OBJDIR)/dos_plotter.o $(OBJDIR)/cli_menu.o \
	    $(OBJDIR)/phonon_dos.o $(OBJDIR)/polarizability.o $(OBJDIR)/poscastep_menu.o \
	    $(OBJDIR)/main.o

clean:
	rm -rf $(OBJDIR) $(TARGET)

run: $(TARGET)
	./$(TARGET)

debug: FCFLAGS += -O0 -fcheck=all -g -fbacktrace
debug: $(TARGET)
