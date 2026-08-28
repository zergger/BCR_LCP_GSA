#!/usr/bin/env python3
"""Generate deterministic reads with known SNP and 1--10 bp InDel events."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import random


CASES = [
    ("snp", 1, 6, 0.00),
    ("ins", 1, 12, 0.01),
    ("ins", 2, 24, 0.03),
    ("ins", 5, 12, 0.00),
    ("ins", 10, 12, 0.00),
    ("del", 1, 24, 0.01),
    ("del", 2, 6, 0.00),
    ("del", 5, 12, 0.03),
    ("del", 10, 24, 0.01),
]
DNA = "ACGT"


def deterministic_dna(label: str, length: int) -> str:
    output = []
    counter = 0
    while len(output) < length:
        digest = hashlib.sha256(f"{label}:{counter}".encode()).digest()
        output.extend(DNA[byte & 3] for byte in digest)
        counter += 1
    return "".join(output[:length])


def add_errors(sequence: str, rate: float, seed: int) -> str:
    rng = random.Random(seed)
    output = []
    for symbol in sequence:
        if rng.random() < rate:
            output.append(rng.choice(DNA.replace(symbol, "")))
        else:
            output.append(symbol)
    return "".join(output)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)

    reference_records = []
    alternate_records = []
    truth = ["event\ttype\tlength\tdepth\terror_rate\tposition\tref\talt"]
    for index, (kind, length, depth, error_rate) in enumerate(CASES, start=1):
        event = f"event_{index:02d}_{kind}{length}"
        reference = deterministic_dna(event, 120)
        position = 60
        if kind == "snp":
            ref_allele = reference[position]
            alt_allele = DNA[(DNA.index(ref_allele) + 1) % len(DNA)]
            alternate = reference[:position] + alt_allele + reference[position + 1 :]
        elif kind == "ins":
            ref_allele = "-"
            alt_allele = deterministic_dna(event + ":insert", length)
            alternate = reference[:position] + alt_allele + reference[position:]
        else:
            ref_allele = reference[position : position + length]
            alt_allele = "-"
            alternate = reference[:position] + reference[position + length :]

        truth.append(
            f"{event}\t{kind}\t{length}\t{depth}\t{error_rate:.2f}\t"
            f"{position}\t{ref_allele}\t{alt_allele}"
        )
        for copy in range(depth):
            seed = index * 100000 + copy
            reference_records.append(
                (f"{event}_ref_{copy:03d}", add_errors(reference, error_rate, seed))
            )
            alternate_records.append(
                (f"{event}_alt_{copy:03d}", add_errors(alternate, error_rate, seed + 50000))
            )

    write_fasta(args.output_dir / "reference.fa", reference_records)
    write_fasta(args.output_dir / "alternate.fa", alternate_records)
    (args.output_dir / "truth.tsv").write_text("\n".join(truth) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
