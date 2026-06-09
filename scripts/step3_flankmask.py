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


def drop_all_gap_columns(seqs: list[SeqRec]) -> tuple[list[SeqRec], int]:
    if not seqs or not seqs[0].seq:
        return seqs, 0
    keep = []
    for i in range(len(seqs[0].seq)):
        if any(rec.seq[i] not in MISSING for rec in seqs):
            keep.append(i)
    dropped = len(seqs[0].seq) - len(keep)
    return [SeqRec(rec.name, "".join(rec.seq[i] for i in keep)) for rec in seqs], dropped


def choose_threshold(ntaxa: int, absolute: int | None, fraction: float | None) -> int:
    if absolute is not None and fraction is not None:
        raise SystemExit("Choose either absolute or fractional flank occupancy threshold, not both.")
    if absolute is not None:
        return absolute
    if fraction is not None:
        return max(1, math.ceil(ntaxa * fraction))
    raise SystemExit("No flank occupancy threshold provided.")


def mask_flanks(seqs: list[SeqRec], min_col_taxa: int, good_run: int) -> tuple[list[SeqRec], int]:
    ncol = len(seqs[0].seq)
    cols = [[rec.seq[i] for rec in seqs] for i in range(ncol)]
    counts_per_col = []
    occupancies = []
    for col in cols:
        counts = {}
        occ = 0
        for c in col:
            if c in MISSING:
                continue
            counts[c] = counts.get(c, 0) + 1
            occ += 1
        counts_per_col.append(counts)
        occupancies.append(occ)

    def is_good(seq_idx: int, pos: int) -> bool:
        base = seqs[seq_idx].seq[pos]
        if base in MISSING or occupancies[pos] < min_col_taxa:
            return False
        return counts_per_col[pos].get(base, 0) > 1

    total_masked = 0
    out = []
    for seq_idx, rec in enumerate(seqs):
        seq = list(rec.seq)
        left_anchor = ncol
        streak = 0
        for pos in range(ncol):
            if is_good(seq_idx, pos):
                streak += 1
                if streak >= good_run:
                    left_anchor = pos - good_run + 1
                    break
            else:
                streak = 0
        right_anchor = -1
        streak = 0
        for pos in range(ncol - 1, -1, -1):
            if is_good(seq_idx, pos):
                streak += 1
                if streak >= good_run:
                    right_anchor = pos + good_run - 1
                    break
            else:
                streak = 0
        for pos in range(0, min(left_anchor, ncol)):
            if seq[pos] not in MISSING:
                seq[pos] = "N"
                total_masked += 1
        for pos in range(max(right_anchor + 1, 0), ncol):
            if seq[pos] not in MISSING:
                seq[pos] = "N"
                total_masked += 1
        out.append(SeqRec(rec.name, "".join(seq)))
    return out, total_masked


def main() -> int:
    ap = argparse.ArgumentParser(description="Flank-only private singleton masker.")
    ap.add_argument("input_dir")
    ap.add_argument("output_dir")
    ap.add_argument("--min-col-occupancy-taxa", type=int)
    ap.add_argument("--min-col-occupancy-frac", type=float, default=0.60)
    ap.add_argument("--good-run", type=int, default=10)
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    rows = []
    for fname in sorted(f for f in os.listdir(args.input_dir) if f.endswith(".fasta")):
        seqs = read_fasta(os.path.join(args.input_dir, fname))
        if not seqs:
            continue
        ntaxa = len(seqs)
        min_col_taxa = choose_threshold(ntaxa, args.min_col_occupancy_taxa, args.min_col_occupancy_frac)
        initial_len = len(seqs[0].seq)
        initial_bases = sum(sum(1 for c in rec.seq if c not in MISSING) for rec in seqs)
        masked, masked_bases = mask_flanks(seqs, min_col_taxa, args.good_run)
        cleaned, dropped_cols = drop_all_gap_columns(masked)
        final_len = len(cleaned[0].seq) if cleaned else 0
        final_bases = sum(sum(1 for c in rec.seq if c not in MISSING) for rec in cleaned)
        write_fasta(os.path.join(args.output_dir, fname), cleaned)
        rows.append({
            "file": fname,
            "n_taxa": ntaxa,
            "min_col_taxa": min_col_taxa,
            "good_run": args.good_run,
            "initial_len": initial_len,
            "masked_flank_bases": masked_bases,
            "all_gap_cols_dropped": dropped_cols,
            "final_len": final_len,
            "initial_real_bases": initial_bases,
            "final_real_bases": final_bases,
        })
    with open(os.path.join(args.output_dir, "flank_mask_summary.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
