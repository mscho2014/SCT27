#!/bin/bash
# =============================================================================
# run_fig1b.sh
# Reproduces Fig. 1(b) of Cho et al., Phys. Rev. E 109, 045207 (2024)
#
# Scenario: Neutral Neon gas (500 torr) irradiated by 2 keV XFEL pulse
#   - Nt = 2e19 /cc, T = 300 K
#   - hv = 2000 eV, FWHM ~ 230 fs, peak at t = 300 fs
#   - Simulates EEDF evolution through photoionization + Auger cascade
#   - Run time: < 600 fs (one simulation)
#
# Expected output: charge state distribution (CSD) and EEDF vs time
# Key result: non-Maxwellian EEDF with photoelectron peak (~1130 eV)
#             and Auger peak (~770 eV)
# =============================================================================

set -e

BINDIR="../"
INPUTDIR="."
WORKDIR="run_fig1b"

mkdir -p $WORKDIR
cd $WORKDIR

# Link input files
ln -sf ../*.dat .       2>/dev/null || true
ln -sf ../atomic.data . 2>/dev/null || true
ln -sf ../atomic_inp.10 atomic.inp 2>/dev/null || true
ln -sf ../fe_200g_2000ev fe_200g_2000ev 2>/dev/null || true
ln -sf ../hvfile hvfile 2>/dev/null || true
ln -sf ../his_ne_29_200step his_ne_29_200step 2>/dev/null || true
ln -sf ../initial_ne.dat initial.dat 2>/dev/null || true
cp ../runfile_fig1b input

echo "=== SCT27 run: Fig.1(b) - XFEL pulse on Ne gas ==="
echo "    ee_scale = 1.0 (physical e-e collisions)"
echo ""

../../sct27 < input 2>&1 | tee run.log

echo ""
echo "=== Run complete. Key output files: ==="
echo "    ne_xfel_out  : charge state distribution vs time"
echo "    zelda.out    : EEDF evolution (Te, Ne per ZELDA step)"
echo "    source_monitor: photoionization/Auger source rates"
echo ""
echo "Save the final EEDF (fe output) for use as input to run_fig1c.sh"
