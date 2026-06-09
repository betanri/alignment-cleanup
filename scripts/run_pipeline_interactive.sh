#!/bin/zsh
set -euo pipefail

BASE_DIR="${1:-$(pwd)}"
INPUT_DIR="${2:-}"
JULIA_BIN="${JULIA_BIN:-/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia}"
TAPER_SCRIPT="${TAPER_SCRIPT:-/path/to/TAPER/correction_multi.jl}"

if [[ -z "${INPUT_DIR}" ]]; then
  echo "Usage: $0 REPO_BASE INPUT_ALIGNMENT_DIR"
  exit 1
fi

echo "Step 1: custom trim"
read "?Use absolute min column occupancy taxa? Leave blank to use percent [default 60%]: " COL_ABS
if [[ -z "${COL_ABS}" ]]; then
  read "?Minimum column occupancy fraction [0.60]: " COL_FRAC
  COL_FRAC="${COL_FRAC:-0.60}"
  STEP1_COL_ARGS=(--min-col-occupancy-frac "$COL_FRAC")
else
  STEP1_COL_ARGS=(--min-col-occupancy-taxa "$COL_ABS")
fi
read "?Use absolute min retained taxa per locus? Leave blank to use percent or default absolute 6: " LOC_ABS
if [[ -z "${LOC_ABS}" ]]; then
  read "?Minimum retained taxa fraction per locus [blank keeps default absolute 6]: " LOC_FRAC
  if [[ -z "${LOC_FRAC}" ]]; then
    STEP1_LOC_ARGS=(--min-locus-taxa 6)
  else
    STEP1_LOC_ARGS=(--min-locus-frac "$LOC_FRAC")
  fi
else
  STEP1_LOC_ARGS=(--min-locus-taxa "$LOC_ABS")
fi
read "?Minimum fragment bp [100]: " MIN_BP
MIN_BP="${MIN_BP:-100}"
read "?Minimum fragment fraction of longest sequence [0.25]: " MIN_FRAC
MIN_FRAC="${MIN_FRAC:-0.25}"

STEP1_OUT="$BASE_DIR/out_step1_trimmed"
STEP2_OUT="$BASE_DIR/out_step2_taper"
STEP3_OUT="$BASE_DIR/out_step3_flankmask"
mkdir -p "$STEP1_OUT" "$STEP2_OUT" "$STEP3_OUT"

python3 "$BASE_DIR/scripts/step1_trim.py" "$INPUT_DIR" "$STEP1_OUT" \
  "${STEP1_COL_ARGS[@]}" "${STEP1_LOC_ARGS[@]}" \
  --min-fragment-bp "$MIN_BP" --min-fragment-frac "$MIN_FRAC"

echo "Step 2: TAPER"
read "?TAPER cutoff [1.5]: " TAPER_CUTOFF
TAPER_CUTOFF="${TAPER_CUTOFF:-1.5}"
LIST_FILE="$STEP2_OUT/taper_file_list.txt"
: > "$LIST_FILE"
for f in "$STEP1_OUT"/*.fasta; do
  base=$(basename "$f")
  printf '%s\n%s\n' "$f" "$STEP2_OUT/$base" >> "$LIST_FILE"
done
"$JULIA_BIN" "$TAPER_SCRIPT" -l -m N -a N -c "$TAPER_CUTOFF" "$LIST_FILE"

echo "Step 3: flank mask"
read "?Use absolute flank occupancy taxa? Leave blank to use percent [default 60%]: " FLANK_ABS
if [[ -z "${FLANK_ABS}" ]]; then
  read "?Flank occupancy fraction [0.60]: " FLANK_FRAC
  FLANK_FRAC="${FLANK_FRAC:-0.60}"
  FLANK_ARGS=(--min-col-occupancy-frac "$FLANK_FRAC")
else
  FLANK_ARGS=(--min-col-occupancy-taxa "$FLANK_ABS")
fi
read "?Required consecutive good bases before stopping flank masking [10]: " GOOD_RUN
GOOD_RUN="${GOOD_RUN:-10}"
python3 "$BASE_DIR/scripts/step3_flankmask.py" "$STEP2_OUT" "$STEP3_OUT" \
  "${FLANK_ARGS[@]}" --good-run "$GOOD_RUN"

echo "Done."
echo "Step 1: $STEP1_OUT"
echo "Step 2: $STEP2_OUT"
echo "Step 3: $STEP3_OUT"
