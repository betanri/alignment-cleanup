#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass

MISSING = set("-Nn?")


@dataclass
class SeqRec:
    name: str
    seq: str


def read_fasta(path: str) -> list[SeqRec]:
    seqs: list[SeqRec] = []
    name = None
    buf: list[str] = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append(SeqRec(name=name, seq="".join(buf)))
                name = line[1:]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs.append(SeqRec(name=name, seq="".join(buf)))
    return seqs


def write_fasta(path: str, seqs: list[SeqRec], wrap: int = 80) -> None:
    with open(path, "w") as handle:
        for rec in seqs:
            handle.write(f">{rec.name}\n")
            for i in range(0, len(rec.seq), wrap):
                handle.write(rec.seq[i : i + wrap] + "\n")


def occupancy(col: list[str]) -> int:
    return sum(1 for c in col if c not in MISSING)


def trim_terminal_columns(seqs: list[SeqRec], min_taxa: int) -> tuple[list[SeqRec], int, int]:
    ncol = len(seqs[0].seq)
    left = 0
    right = ncol - 1
    while left < ncol:
        if occupancy([rec.seq[left] for rec in seqs]) >= min_taxa:
            break
        left += 1
    while right >= left:
        if occupancy([rec.seq[right] for rec in seqs]) >= min_taxa:
            break
        right -= 1
    if right < left:
        return [SeqRec(rec.name, "") for rec in seqs], left, ncol
    return [SeqRec(rec.name, rec.seq[left:right + 1]) for rec in seqs], left, ncol - right - 1


def nongap_len(seq: str) -> int:
    return sum(1 for c in seq if c not in MISSING)


def mask_short_fragments(seqs: list[SeqRec], min_bp: int, min_frac: float) -> tuple[list[SeqRec], int, int]:
    lengths = [nongap_len(rec.seq) for rec in seqs]
    longest = max(lengths) if lengths else 0
    threshold = max(min_bp, math.ceil(longest * min_frac))
    masked = 0
    out: list[SeqRec] = []
    for rec in seqs:
        if nongap_len(rec.seq) < threshold:
            out.append(SeqRec(rec.name, "-" * len(rec.seq)))
            masked += 1
        else:
            out.append(rec)
    return out, masked, threshold


def drop_all_gap_columns(seqs: list[SeqRec]) -> tuple[list[SeqRec], int]:
    if not seqs or not seqs[0].seq:
        return seqs, 0
    keep = []
    for i in range(len(seqs[0].seq)):
        if any(rec.seq[i] not in MISSING for rec in seqs):
            keep.append(i)
    dropped = len(seqs[0].seq) - len(keep)
    return [SeqRec(rec.name, "".join(rec.seq[i] for i in keep)) for rec in seqs], dropped


def choose_threshold(ntaxa: int, absolute: int | None, fraction: float | None, label: str) -> int:
    if absolute is not None and fraction is not None:
        raise SystemExit(f"Choose either absolute or fractional threshold for {label}, not both.")
    if absolute is not None:
        return absolute
    if fraction is not None:
        return max(1, math.ceil(ntaxa * fraction))
    raise SystemExit(f"No threshold provided for {label}.")


def main() -> int:
    ap = argparse.ArgumentParser(description="Custom end-trim and occupancy filter for locus alignments.")
    ap.add_argument("input_dir")
    ap.add_argument("output_dir")
    ap.add_argument("--min-col-occupancy-taxa", type=int)
    ap.add_argument("--min-col-occupancy-frac", type=float, default=0.60)
    ap.add_argument("--min-locus-taxa", type=int, default=6)
    ap.add_argument("--min-locus-frac", type=float)
    ap.add_argument("--min-fragment-bp", type=int, default=100)
    ap.add_argument("--min-fragment-frac", type=float, default=0.25)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    dropped_dir = os.path.join(args.output_dir, "dropped_loci")
    os.makedirs(dropped_dir, exist_ok=True)

    rows = []
    for fname in sorted(f for f in os.listdir(args.input_dir) if f.endswith(".fasta")):
        seqs = read_fasta(os.path.join(args.input_dir, fname))
        if not seqs:
            continue
        ntaxa = len(seqs)
        min_col_taxa = choose_threshold(ntaxa, args.min_col_occupancy_taxa, args.min_col_occupancy_frac, "column occupancy")
        min_locus_taxa = choose_threshold(ntaxa, args.min_locus_taxa, args.min_locus_frac, "locus retention")
        initial_len = len(seqs[0].seq)
        seqs1, left_trim, right_trim = trim_terminal_columns(seqs, min_col_taxa)
        seqs2, masked, threshold = mask_short_fragments(seqs1, args.min_fragment_bp, args.min_fragment_frac)
        seqs3, all_gap_cols_dropped = drop_all_gap_columns(seqs2)
        kept_taxa = sum(1 for rec in seqs3 if nongap_len(rec.seq) > 0)
        final_len = len(seqs3[0].seq) if seqs3 else 0
        kept = kept_taxa >= min_locus_taxa and final_len > 0
        outdir = args.output_dir if kept else dropped_dir
        write_fasta(os.path.join(outdir, fname), seqs3)
        rows.append({
            "file": fname,
            "n_taxa": ntaxa,
            "min_col_taxa": min_col_taxa,
            "min_locus_taxa": min_locus_taxa,
            "initial_len": initial_len,
            "left_trim": left_trim,
            "right_trim": right_trim,
            "masked_short_fragments": masked,
            "fragment_threshold_bp": threshold,
            "all_gap_cols_dropped": all_gap_cols_dropped,
            "kept_taxa": kept_taxa,
            "final_len": final_len,
            "kept_locus": kept,
        })

    with open(os.path.join(args.output_dir, "trim_summary.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
