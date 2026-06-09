#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs = []
    name = None
    buf: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf)))
            name = line[1:]
            buf = []
        else:
            buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def write_fasta(path: Path, seqs: list[tuple[str, str]]) -> None:
    with open(path, "w") as handle:
        for name, seq in seqs:
            handle.write(f">{name}\n{seq}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-nontrimmed", required=True)
    ap.add_argument("--source-trimmed", required=True)
    ap.add_argument("--source-taper", required=True)
    ap.add_argument("--source-flankmask", required=True)
    ap.add_argument("--repo-data-dir", required=True)
    ap.add_argument("--n-loci", type=int, default=50)
    args = ap.parse_args()

    steps = {
        "01_Gene_Alignments_Non-Trimmed": Path(args.source_nontrimmed),
        "02_Gene_Alignments_TRIMMED_CODEX": Path(args.source_trimmed),
        "03_Gene_Alignments_TRIMMED_CODEX_TAPER": Path(args.source_taper),
        "04_Gene_Alignments_TRIMMED_CODEX_TAPER_FLANKMASK": Path(args.source_flankmask),
    }
    repo_data_dir = Path(args.repo_data_dir)

    final_files = sorted(p.name for p in steps["04_Gene_Alignments_TRIMMED_CODEX_TAPER_FLANKMASK"].glob("*.fasta"))
    chosen = final_files[: args.n_loci]
    first_locus = read_fasta(steps["04_Gene_Alignments_TRIMMED_CODEX_TAPER_FLANKMASK"] / chosen[0])
    orig_names = [name for name, _ in first_locus]
    if len(orig_names) != 10:
        raise SystemExit("Expected 10 taxa in the example set; update this helper before publishing.")
    rename = {name: f"Sp{i+1}_ind1" for i, name in enumerate(orig_names)}

    for outname, srcdir in steps.items():
        outdir = repo_data_dir / outname
        outdir.mkdir(parents=True, exist_ok=True)
        for fname in chosen:
            records = read_fasta(srcdir / fname)
            anon = [(rename[name], seq) for name, seq in records]
            write_fasta(outdir / fname, anon)

    with open(repo_data_dir / "example_loci.txt", "w") as handle:
        for fname in chosen:
            handle.write(fname + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
