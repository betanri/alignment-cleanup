#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from collections import OrderedDict


def read_fasta(path: str) -> OrderedDict[str, str]:
    data: OrderedDict[str, str] = OrderedDict()
    name = None
    buf = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    data[name] = "".join(buf)
                name = line[1:]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        data[name] = "".join(buf)
    return data


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_dir")
    ap.add_argument("output_fasta")
    args = ap.parse_args()
    combined: OrderedDict[str, str] = OrderedDict()
    for fname in sorted(f for f in os.listdir(args.input_dir) if f.endswith(".fasta")):
        aln = read_fasta(os.path.join(args.input_dir, fname))
        locus_len = len(next(iter(aln.values()))) if aln else 0
        for k in aln:
            if k not in combined:
                combined[k] = ""
        for k in list(combined.keys()):
            combined[k] += aln.get(k, "-" * locus_len)
    with open(args.output_fasta, "w") as handle:
        for name, seq in combined.items():
            handle.write(f">{name}\n{seq}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
