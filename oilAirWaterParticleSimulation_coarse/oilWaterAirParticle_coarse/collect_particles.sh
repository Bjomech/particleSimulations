#!/bin/bash
# Run this from your CASE directory after a run. Pulls the Lagrangian
# particle cloud data (positions, velocities) from EVERY saved timestep so
# the actual particle trajectories can be reconstructed and checked, not
# just guessed from a ParaView animation. Attach the resulting
# particles.tar.gz back to me.

set -e
OUT=particles_collect
rm -rf "$OUT" particles.tar.gz
mkdir -p "$OUT"

TIMES=$(ls -d [0-9]*.[0-9]* [0-9]* 2>/dev/null | grep -E "^[0-9.]+$" | sort -gu)
echo "$TIMES" > "$OUT/times_list.txt"

# Find the cloud name automatically (usually kinematicCloud) instead of
# hardcoding it, in case it differs.
CLOUD_DIR=""
for t in $TIMES; do
    if [ -d "$t/lagrangian" ]; then
        CLOUD_DIR=$(ls "$t/lagrangian" | head -1)
        break
    fi
done
echo "Detected cloud directory name: $CLOUD_DIR" >> "$OUT/times_list.txt"

for t in $TIMES; do
    SRC="$t/lagrangian/$CLOUD_DIR"
    if [ -d "$SRC" ]; then
        mkdir -p "$OUT/$t"
        cp -r "$SRC"/* "$OUT/$t/" 2>/dev/null || true
    fi
done

tar -czf particles.tar.gz -C "$OUT" .
rm -rf "$OUT"
echo "Done. Attach: $(pwd)/particles.tar.gz"
