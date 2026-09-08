#!/bin/bash
# Run this from your CASE directory after a completed run (e.g. the
# kickstart case oilWaterAirParticle_v_0). Pulls alpha.air + U from EVERY
# saved timestep so the full bubble trajectory can be reconstructed, not
# just a single snapshot. Attach the resulting trajectory.tar.gz back to me.

set -e
OUT=trajectory_collect
rm -rf "$OUT" trajectory.tar.gz
mkdir -p "$OUT"

TIMES=$(ls -d [0-9]*.[0-9]* [0-9]* 2>/dev/null | grep -E "^[0-9.]+$" | sort -gu)
echo "$TIMES" > "$OUT/times_list.txt"

for t in $TIMES; do
    if [ -f "$t/alpha.air" ] && [ -f "$t/U" ]; then
        mkdir -p "$OUT/$t"
        cp "$t/alpha.air" "$t/U" "$OUT/$t/"
    fi
done

tar -czf trajectory.tar.gz -C "$OUT" .
rm -rf "$OUT"
echo "Done. Attach: $(pwd)/trajectory.tar.gz"
