# ===========================================================================
# Makefile for SCT27 (SCFLY + ZELDA)
# NLTE atomic kinetics + Boltzmann EEDF solver for XFEL plasma experiments
#
# Reference: Cho et al., Phys. Rev. E 109, 045207 (2024)
#
# Compiler: gfortran (tested with gfortran >= 9.0)
# ===========================================================================

FC = gfortran

# ---------------------------------------------------------------------------
# Compiler flag sets
# ---------------------------------------------------------------------------

# DEBUG: full bounds checking, warnings, no optimization
FLAGS_DEBUG   = -g -c -w -O0 -fPIC -fbounds-check -fcheck=all -Wall

# STATIC: static storage (required for legacy common block code)
FLAGS_STATIC  = -g -c -w -static

# RELEASE: optimized build (use after code is validated)
FLAGS_RELEASE = -O2 -c -w -fPIC

# Default: debug mode for reproducibility validation
FLAGS_SCFLY   = $(FLAGS_DEBUG)
FLAGS_ZELDA   = $(FLAGS_STATIC)

# ---------------------------------------------------------------------------
# Source files
# ---------------------------------------------------------------------------

F_SRC = scmn_v3.f   \
        scdrv.f     \
        sclsd.f     \
        scslv.f     \
        scstr.f     \
        sctrn.f     \
        zchkxc.f    \
        zchkfly.f   \
        zmodules.f90        \
        zchk_stuff.f90      \
        zcoeff_elastic.f90  \
        zelda.f90

F_OBJ = scmn_v3.o   \
        scdrv.o     \
        sclsd.o     \
        scslv.o     \
        scstr.o     \
        sctrn.o     \
        zchkxc.o    \
        zchkfly.o   \
        zmodules.o          \
        zchk_stuff.o        \
        zcoeff_elastic.o    \
        zelda.o

# Include files (dependencies, not compiled directly)
F_INC = mainstuf  \
        popstuf   \
        runstuf   \
        timestuf  \
        xc_stuf   \
        flystuf   \
        kinstuf

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

.PHONY: all clean debug release

all: sct27

sct27: $(F_OBJ)
	$(FC) -fPIC -o sct27 $(F_OBJ)
	@echo "Build complete: ./sct27"

# Alias for debug build (default)
debug: sct27

# Optimized build
release: FLAGS_SCFLY  = $(FLAGS_RELEASE)
release: FLAGS_ZELDA  = $(FLAGS_RELEASE)
release: clean sct27

# ---------------------------------------------------------------------------
# Individual file compilation rules
# ---------------------------------------------------------------------------

# SCFLY core (NLTE population solver)
scmn_v3.o: $(F_INC) scmn_v3.f
	$(FC) $(FLAGS_SCFLY) scmn_v3.f

scdrv.o: $(F_INC) scdrv.f
	$(FC) $(FLAGS_SCFLY) scdrv.f

sclsd.o: $(F_INC) sclsd.f
	$(FC) $(FLAGS_DEBUG) sclsd.f

scslv.o: $(F_INC) scslv.f
	$(FC) $(FLAGS_SCFLY) scslv.f

scstr.o: $(F_INC) scstr.f
	$(FC) $(FLAGS_SCFLY) scstr.f

sctrn.o: $(F_INC) sctrn.f
	$(FC) $(FLAGS_SCFLY) sctrn.f

# ZELDA interface
zchkfly.o: $(F_INC) zchkfly.f
	$(FC) $(FLAGS_ZELDA) zchkfly.f

zchkxc.o: $(F_INC) zchkxc.f
	$(FC) $(FLAGS_ZELDA) zchkxc.f

# ZELDA Boltzmann solver (Fortran 90 modules - ORDER MATTERS)
zmodules.o: zmodules.f90
	$(FC) $(FLAGS_ZELDA) zmodules.f90

zchk_stuff.o: zchk_stuff.f90 zmodules.o
	$(FC) $(FLAGS_ZELDA) zchk_stuff.f90

zcoeff_elastic.o: zcoeff_elastic.f90 zmodules.o zchk_stuff.o
	$(FC) $(FLAGS_ZELDA) zcoeff_elastic.f90

zelda.o: zelda.f90 zmodules.o zchk_stuff.o
	$(FC) $(FLAGS_ZELDA) zelda.f90

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

clean:
	rm -f *.o *.mod sct27
	@echo "Cleaned object files and executable"
