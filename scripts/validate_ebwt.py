#!/usr/bin/env python3
"""Validate a FASTA/FASTQ input and its BCR eBWT for ebwt2InDel."""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
from pathlib import Path
import sys

from bwt_format import BwtReader, read_all


DNA = frozenset("ACGT")


def open_text(path: Path):
    with path.open("rb") as raw:
        is_gzip = raw.read(2) == b"\x1f\x8b"
    if is_gzip:
        return gzip.open(path, "rt", encoding="ascii", newline=None)
    return path.open("rt", encoding="ascii", newline=None)


def detect_format(path: Path) -> str:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            raise ValueError("input is neither FASTA nor FASTQ")
    raise ValueError("input contains no records")


def checked_sequence_line(line: str, line_number: int) -> bytes:
    invalid = sorted(set(line) - DNA)
    if invalid:
        shown = ", ".join(repr(char) for char in invalid[:8])
        raise ValueError(
            f"ebwt2InDel accepts uppercase A/C/G/T only; line {line_number} contains {shown}"
        )
    return line.encode("ascii")


def sequence_records(path: Path, input_format: str):
    with open_text(path) as handle:
        lines = iter(enumerate(handle, start=1))

        first = None
        for item in lines:
            if item[1].rstrip("\r\n"):
                first = item
                break
        if first is None:
            raise ValueError("input contains no records")

        if input_format == "fasta":
            if not first[1].startswith(">"):
                raise ValueError("FASTA input does not begin with '>'")
            current = bytearray()
            for line_number, raw_line in lines:
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if not current:
                        raise ValueError(f"empty FASTA sequence before line {line_number}")
                    yield bytes(current)
                    current.clear()
                else:
                    current.extend(checked_sequence_line(line, line_number))
            if not current:
                raise ValueError("the final FASTA sequence is empty")
            yield bytes(current)
            return

        header = first
        while header is not None:
            header_number, header_line = header
            if not header_line.startswith("@"):
                raise ValueError(f"FASTQ record at line {header_number} does not begin with '@'")

            sequence = bytearray()
            plus_found = False
            for line_number, raw_line in lines:
                line = raw_line.rstrip("\r\n")
                if line.startswith("+"):
                    plus_found = True
                    break
                if line:
                    sequence.extend(checked_sequence_line(line, line_number))
            if not plus_found:
                raise ValueError(f"FASTQ record at line {header_number} has no '+' separator")
            if not sequence:
                raise ValueError(f"FASTQ record at line {header_number} has an empty sequence")

            quality_length = 0
            for line_number, raw_line in lines:
                quality_length += len(raw_line.rstrip("\r\n"))
                if quality_length >= len(sequence):
                    break
            if quality_length != len(sequence):
                raise ValueError(
                    f"FASTQ record at line {header_number} has sequence/quality lengths "
                    f"{len(sequence)}/{quality_length}"
                )
            yield bytes(sequence)

            header = None
            for item in lines:
                if item[1].rstrip("\r\n"):
                    header = item
                    break


def sequence_stats(path: Path, max_read_length: int) -> tuple[str, int, int, int]:
    input_format = detect_format(path)
    sequences = 0
    total_bases = 0
    longest = 0
    for sequence in sequence_records(path, input_format):
        length = len(sequence)
        if length > max_read_length:
            raise ValueError(f"read length {length} exceeds this build's {max_read_length}-base limit")
        sequences += 1
        total_bases += length
        longest = max(longest, length)
    if sequences == 0:
        raise ValueError("input contains no sequences")
    return input_format, sequences, total_bases, longest


def ebwt_stats(path: Path, terminator: str) -> tuple[int, int, str, int]:
    allowed = {ord(base) for base in DNA}
    allowed.add(ord(terminator))
    size = 0
    terminators = 0

    with BwtReader(path, terminator) as handle:
        encoding = handle.encoding
        stored_bytes = path.stat().st_size
        for chunk in handle.chunks():
            size += len(chunk)
            invalid = sorted(set(chunk) - allowed)
            if invalid:
                shown = ", ".join(str(value) for value in invalid[:8])
                raise ValueError(f"eBWT contains unsupported byte value(s): {shown}")
            terminators += chunk.count(ord(terminator))
    return size, terminators, encoding, stored_bytes


def input_multiset(path: Path) -> Counter[bytes]:
    return Counter(sequence_records(path, detect_format(path)))


def invert_ebwt(path: Path, terminator: str, max_bytes: int) -> Counter[bytes]:
    bwt = read_all(path, terminator, max_bytes)
    term = ord(terminator)
    counts = [0] * 256
    for symbol in bwt:
        counts[symbol] += 1

    starts = [0] * 256
    cumulative = 0
    for symbol, count in enumerate(counts):
        starts[symbol] = cumulative
        cumulative += count

    seen = [0] * 256
    lf = [0] * len(bwt)
    terminator_rows = []
    for position, symbol in enumerate(bwt):
        lf[position] = starts[symbol] + seen[symbol]
        seen[symbol] += 1
        if symbol == term:
            terminator_rows.append(position)

    reconstructed: Counter[bytes] = Counter()
    for start in terminator_rows:
        position = start
        reversed_sequence = bytearray()
        first_terminator = True
        for _ in range(len(bwt) + 1):
            symbol = bwt[position]
            position = lf[position]
            if symbol == term:
                if first_terminator:
                    first_terminator = False
                    continue
                reconstructed[bytes(reversed(reversed_sequence))] += 1
                break
            reversed_sequence.append(symbol)
        else:
            raise ValueError("eBWT inversion did not close a terminator cycle")
    return reconstructed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path, help="input FASTA/FASTQ, optionally gzip-compressed")
    parser.add_argument("--ebwt", type=Path, help="BCR eBWT to validate")
    parser.add_argument("--terminator", default="$", help="single-byte BCR terminator (default: $)")
    parser.add_argument("--max-read-length", type=int, default=254)
    parser.add_argument("--verify-read-multiset", action="store_true")
    parser.add_argument("--max-inversion-bytes", type=int, default=16 * 1024 * 1024)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if len(args.terminator.encode("ascii", errors="strict")) != 1:
        raise ValueError("--terminator must be one ASCII byte")
    if args.max_read_length < 1:
        raise ValueError("--max-read-length must be positive")

    input_format, sequences, bases, longest = sequence_stats(args.input, args.max_read_length)
    fields = [
        f"input_format={input_format}",
        f"input_sequences={sequences}",
        f"input_bases={bases}",
        f"max_read_length={longest}",
    ]

    if args.ebwt is not None:
        ebwt_bytes, terminators, encoding, stored_bytes = ebwt_stats(
            args.ebwt, args.terminator
        )
        expected_bytes = bases + sequences
        if ebwt_bytes != expected_bytes:
            raise ValueError(f"eBWT size is {ebwt_bytes}, expected {expected_bytes}")
        if terminators != sequences:
            raise ValueError(f"eBWT has {terminators} terminators, expected {sequences}")
        fields.extend(
            [
                f"ebwt_bytes={ebwt_bytes}",
                f"ebwt_encoding={encoding}",
                f"ebwt_stored_bytes={stored_bytes}",
                f"terminator_ascii={ord(args.terminator)}",
                f"terminators={terminators}",
            ]
        )
        if args.verify_read_multiset:
            if args.max_inversion_bytes < 1:
                raise ValueError("--max-inversion-bytes must be positive")
            expected = input_multiset(args.input)
            observed = invert_ebwt(args.ebwt, args.terminator, args.max_inversion_bytes)
            if observed != expected:
                missing = sum((expected - observed).values())
                extra = sum((observed - expected).values())
                raise ValueError(f"inverted read multiset differs: missing={missing}, extra={extra}")
            fields.append("read_multiset_verified=1")

    print("\t".join(fields))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, UnicodeError, ValueError) as error:
        print(f"validate_ebwt.py: {error}", file=sys.stderr)
        raise SystemExit(1)
