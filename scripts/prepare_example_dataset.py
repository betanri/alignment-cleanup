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
    rename = {
        "Galeichthys_ater_RB0513": "Sp1_ind1",
        "Galeichthys_feliceps_RB0508": "Sp2_ind1",
        "Galeichthys_feliceps_RB0509": "Sp2_ind2",
        "Galeichthys_feliceps_RB0510": "Sp2_ind3",
        "Galeichthys_feliceps_RB0512": "Sp2_ind4",
        "Galeichthys_trowi_RB0504": "Sp3_ind1",
        "Galeichthys_trowi_RB0505": "Sp3_ind2",
        "Galeichthys_trowi_RB0507": "Sp3_ind3",
        "Galeichthys_peruvianus_EPLATE_86_C03": "Sp4_ind1",
        "Neoarius_graeffei": "Sp5_ind1",
    }
    if set(orig_names) != set(rename):
        raise SystemExit("Unexpected taxon set; update anonymization mapping before publishing.")

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
