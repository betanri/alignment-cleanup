#!/usr/bin/env python3
from __future__ import annotations

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_tree")
    ap.add_argument("output_png")
    ap.add_argument("--title", default="")
    ap.add_argument("--outgroup", default="")
    args = ap.parse_args()
    tree = Phylo.read(args.input_tree, "newick")
    if args.outgroup:
        matches = [clade for clade in tree.find_clades() if getattr(clade, "name", None) == args.outgroup]
        if not matches:
            raise SystemExit(f"Outgroup '{args.outgroup}' not found in tree.")
        tree.root_with_outgroup(matches[0])
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    ax.set_title(args.title, fontsize=14)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=180)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
