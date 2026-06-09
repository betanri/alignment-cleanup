#!/bin/zsh
set -euo pipefail

REPO_DIR="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
IQTREE_BIN="${IQTREE_BIN:-/Applications/iqtree-2.4.0-macOS/iqtree2}"
FASTTREE_BIN="${FASTTREE_BIN:-$(command -v FastTree || true)}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/private/tmp/mplconfig}"
export MPLCONFIGDIR
mkdir -p "$MPLCONFIGDIR"

DATA_DIR="$REPO_DIR/data/examples"
TREE_DIR="$REPO_DIR/docs/trees"
FIG_DIR="$REPO_DIR/docs/figures"
mkdir -p "$TREE_DIR" "$FIG_DIR"

for STEP in \
  01_Gene_Alignments_Non-Trimmed \
  02_Gene_Alignments_TRIMMED_CODEX \
  03_Gene_Alignments_TRIMMED_CODEX_TAPER \
  04_Gene_Alignments_TRIMMED_CODEX_TAPER_FLANKMASK
do
  CONCAT="$TREE_DIR/${STEP}.concat.fasta"
  python3 "$REPO_DIR/scripts/concatenate_alignments.py" "$DATA_DIR/$STEP" "$CONCAT"
  PREFIX="$TREE_DIR/${STEP}"
  if [[ -x "$IQTREE_BIN" ]]; then
    "$IQTREE_BIN" -s "$CONCAT" -m GTR+G -nt AUTO -seed 12345 -redo --prefix "$PREFIX"
    TREEFILE="${PREFIX}.treefile"
  elif [[ -n "$FASTTREE_BIN" ]]; then
    "$FASTTREE_BIN" -gtr -nt "$CONCAT" > "${PREFIX}.treefile"
    TREEFILE="${PREFIX}.treefile"
  else
    echo "No IQ-TREE or FastTree executable found." >&2
    exit 1
  fi
  python3 "$REPO_DIR/scripts/draw_tree.py" "$TREEFILE" "$FIG_DIR/${STEP}.png" --title "$STEP"
done
