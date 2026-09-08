#!/bin/bash
# Run this from your CASE directory (e.g. oilWaterAirParticle_quicktest,
# oilWaterAirParticle_kickstart, or wherever you just ran the solver).
# It gathers everything needed to diagnose the bubble-rise issue into one
# file: diagnostics.tar.gz -- just attach that single file back to me.

set -e
OUT=diagnostics_collect
rm -rf "$OUT" diagnostics.tar.gz
mkdir -p "$OUT"

echo "=== pwd ===" > "$OUT/summary.txt"
pwd >> "$OUT/summary.txt"

echo -e "\n=== constant/g ===" >> "$OUT/summary.txt"
cat constant/g >> "$OUT/summary.txt" 2>&1 || echo "NOT FOUND" >> "$OUT/summary.txt"

echo -e "\n=== system/controlDict (key entries) ===" >> "$OUT/summary.txt"
grep -E "endTime|deltaT|writeInterval|writeControl" system/controlDict >> "$OUT/summary.txt" 2>&1

echo -e "\n=== system/setFieldsDict ===" >> "$OUT/summary.txt"
cat system/setFieldsDict >> "$OUT/summary.txt" 2>&1 || echo "NOT FOUND" >> "$OUT/summary.txt"

echo -e "\n=== constant/transportProperties ===" >> "$OUT/summary.txt"
cat constant/transportProperties >> "$OUT/summary.txt" 2>&1 || echo "NOT FOUND" >> "$OUT/summary.txt"

# Find the solver source used to build the currently-running binary.
SOLVER_DIR=$(which oilWaterAirParticleFoam 2>/dev/null || true)
echo -e "\n=== which oilWaterAirParticleFoam ===" >> "$OUT/summary.txt"
echo "$SOLVER_DIR" >> "$OUT/summary.txt"

echo -e "\n=== pEqn.H phig block (from solver source dir, if you're in it or its sibling) ===" >> "$OUT/summary.txt"
for d in ../oilWaterAirParticleFoam_final ../oilWaterAirParticleFoam_solver_final .; do
    if [ -f "$d/pEqn.H" ]; then
        echo "-- found at $d/pEqn.H --" >> "$OUT/summary.txt"
        grep -n -B2 -A5 "surfaceScalarField phig" "$d/pEqn.H" >> "$OUT/summary.txt"
        break
    fi
done

# List available result times and pick the latest.
echo -e "\n=== available time directories ===" >> "$OUT/summary.txt"
TIMES=$(ls -d [0-9]*.[0-9]* [0-9]* 2>/dev/null | grep -E "^[0-9.]+$" | sort -g)
echo "$TIMES" >> "$OUT/summary.txt"
LATEST=$(echo "$TIMES" | tail -1)
echo "Latest time: $LATEST" >> "$OUT/summary.txt"

# Copy the latest time directory's key fields + the log.
if [ -n "$LATEST" ] && [ -d "$LATEST" ]; then
    mkdir -p "$OUT/latestTime"
    cp "$LATEST"/alpha.air "$LATEST"/alpha.oil "$LATEST"/alpha.water "$LATEST"/U "$LATEST"/p_rgh "$LATEST"/rho "$OUT/latestTime/" 2>/dev/null || true
fi

if [ -f log.solver ]; then
    cp log.solver "$OUT/log.solver"
elif ls log.oilWaterAirParticleFoam* >/dev/null 2>&1; then
    cp log.oilWaterAirParticleFoam* "$OUT/log.solver"
fi

tar -czf diagnostics.tar.gz -C "$OUT" .
rm -rf "$OUT"
echo "Done. Attach: $(pwd)/diagnostics.tar.gz"
