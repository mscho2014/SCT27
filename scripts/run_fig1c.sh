#!/bin/bash
# =============================================================================
# run_fig1c.sh
# Reproduces Fig. 1(c) of Cho et al., Phys. Rev. E 109, 045207 (2024)
#
# This is a SEPARATE simulation from Fig.1(b).
#
# Method (from paper):
#   1. Take the final EEDF from the <600 fs simulation (run_fig1b) as input.
#   2. Run with ee_scale ~ 100 (artificially accelerated e-e collisions).
#   3. Physical simulation time is ~100x shorter (e.g., ~200 fs range).
#   4. Time axis is linearly rescaled by ee_scale to show ~20 ps thermalization.
#
# This approach was necessary because direct simulation to 20 ps was
# computationally prohibitive (time steps ~ 0.01 fs required).
#
# To use:
#   1. First run run_fig1b.sh and obtain the final EEDF output file.
#   2. Copy/rename that EEDF file to fe_thermalization (edit below).
#   3. Edit ee_scale in src/scmn_v3.f (line ~88): ee_scale = 100.0d0
#   4. Rebuild: cd .. && make clean && make
#   5. Run this script.
# =============================================================================

set -e

WORKDIR="run_fig1c"
EE_SCALE=100.0    # Must match value set in scmn_v3.f

# Input EEDF file: output from the end of the Fig.1(b) run
# Rename/copy your EEDF output here:
FE_INPUT="fe_after_pulse"   # <-- Edit this to match your actual filename

if [ ! -f "../input/neon_2keV/$FE_INPUT" ]; then
  echo "ERROR: EEDF input file not found: input/neon_2keV/$FE_INPUT"
  echo "Please run run_fig1b.sh first and copy the final EEDF output."
  echo "Then set FE_INPUT in this script accordingly."
  exit 1
fi

echo "=== SCT27 run: Fig.1(c) - EEDF thermalization study ==="
echo "    ee_scale = $EE_SCALE (accelerated e-e collisions)"
echo "    Input EEDF: $FE_INPUT"
echo "    NOTE: Time axis must be divided by $EE_SCALE to get physical time"
echo ""
echo "  IMPORTANT: verify that ee_scale = ${EE_SCALE}d0 is set in src/scmn_v3.f"
echo "  and that the code was rebuilt (make clean && make)"
echo ""

mkdir -p $WORKDIR
cd $WORKDIR

ln -sf ../../input/neon_2keV/*.dat .  2>/dev/null || true
ln -sf ../../input/neon_2keV/atomic.data . 2>/dev/null || true
ln -sf ../../input/neon_2keV/atomic_inp.29 atomic.inp 2>/dev/null || true
ln -sf "../../input/neon_2keV/$FE_INPUT" fe_thermalization 2>/dev/null || true

cat > input << EOF
z 29
outfile ne_therm_out
initial file initial_ne.dat
fe file fe_thermalization
evolve td
history historyfile his_ne_29_200step rho
time pump 1.d-18 1.d-15 1.e-14 36 34
end
EOF

../../sct27 < input 2>&1 | tee run.log

echo ""
echo "=== Thermalization run complete ==="
echo "    Divide time axis by ee_scale=$EE_SCALE to recover physical time."
echo "    At t_physical ~ 20 ps, EEDF should converge to Maxwell-Boltzmann."
