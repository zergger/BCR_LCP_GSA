#!/usr/bin/env bash
set -Eeuo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
runner="$repo_dir/scripts/run_bcr_for_ebwt2indel.sh"
validator="$repo_dir/scripts/validate_ebwt.py"
fixture="$repo_dir/test/ebwt2indel_smoke.fa"
expected_sha256=ef62854100af8d59f557dd153ec79007e218582d31102853c0843b2b8dc7f0d4
expected_fastq_sha256=c1bdfb98d5dd64bd5722f412653781eb9088aa22a2952195f0f12d5c395ffdea
fastq_fixture="$repo_dir/test/test.fq.gz"

work_dir=$(mktemp -d -p /tmp bcr-ebwt-test.XXXXXX)
cleanup() {
    find "$work_dir" -type f -delete
    find "$work_dir" -depth -type d -empty -delete
}
trap cleanup EXIT

expect_failure() {
    local label=$1
    shift
    if "$@" >"$work_dir/$label.stdout" 2>"$work_dir/$label.stderr"; then
        echo "Expected failure did not occur: $label" >&2
        exit 1
    fi
}

"$runner" "$fixture" "$work_dir/valid"
python3 "$validator" --input "$fixture" --ebwt "$work_dir/valid.rl_bwt" \
    --terminator '$' --verify-read-multiset
"$runner" "$fixture" "$work_dir/plain" plain
"$runner" "$fixture" "$work_dir/packed" packed
python3 "$validator" --input "$fixture" --ebwt "$work_dir/packed.pck_bwt" \
    --terminator '$' --verify-read-multiset
actual_sha256=$(sha256sum "$work_dir/plain.ebwt" | awk '{print $1}')
if [[ "$actual_sha256" != "$expected_sha256" ]]; then
    echo "Unexpected eBWT SHA-256: $actual_sha256" >&2
    exit 1
fi
test ! -e "$work_dir/valid.len"
test ! -e "$work_dir/valid.info"
awk -F '\t' '$1 == "peak" { found = 1; if ($2 <= 0) exit 1 } END { if (!found) exit 1 }' \
    "$work_dir/valid.bcr.space.tsv"

PYTHONPATH="$repo_dir/scripts" python3 - "$work_dir" <<'PY'
from pathlib import Path
import sys
from bwt_format import read_all

root = Path(sys.argv[1])
plain = read_all(root / "plain.ebwt", "$", 16 * 1024 * 1024)
assert read_all(root / "valid.rl_bwt", "$", 16 * 1024 * 1024) == plain
assert read_all(root / "packed.pck_bwt", "$", 16 * 1024 * 1024) == plain

rle = bytearray((root / "valid.rl_bwt").read_bytes())
(root / "truncated.rl_bwt").write_bytes(rle[:-1])
rle[24] ^= 1
(root / "bad_checksum.rl_bwt").write_bytes(rle)

packed = bytearray((root / "packed.pck_bwt").read_bytes())
packed[32] = (packed[32] & 0xF8) | 0x07
(root / "bad_symbol.pck_bwt").write_bytes(packed)
PY

expect_failure truncated_rle \
    python3 "$validator" --input "$fixture" --ebwt "$work_dir/truncated.rl_bwt" --terminator '$'
expect_failure bad_checksum \
    python3 "$validator" --input "$fixture" --ebwt "$work_dir/bad_checksum.rl_bwt" --terminator '$'
expect_failure bad_packed_symbol \
    python3 "$validator" --input "$fixture" --ebwt "$work_dir/bad_symbol.pck_bwt" --terminator '$'

"$runner" "$fastq_fixture" "$work_dir/fastq"
"$runner" "$fastq_fixture" "$work_dir/fastq_plain" plain
python3 "$validator" --input "$fastq_fixture" --ebwt "$work_dir/fastq.rl_bwt" \
    --terminator '$' --verify-read-multiset
actual_fastq_sha256=$(sha256sum "$work_dir/fastq_plain.ebwt" | awk '{print $1}')
if [[ "$actual_fastq_sha256" != "$expected_fastq_sha256" ]]; then
    echo "Unexpected FASTQ eBWT SHA-256: $actual_fastq_sha256" >&2
    exit 1
fi

expect_failure wrong_terminator \
    python3 "$validator" --input "$fixture" --ebwt "$work_dir/valid.rl_bwt" --terminator '#'
expect_failure overwrite "$runner" "$fixture" "$work_dir/valid"

printf '>invalid_n\nACNT\n' >"$work_dir/invalid_n.fa"
expect_failure invalid_n python3 "$validator" --input "$work_dir/invalid_n.fa"

printf '>invalid_lower\nACgT\n' >"$work_dir/invalid_lower.fa"
expect_failure invalid_lower python3 "$validator" --input "$work_dir/invalid_lower.fa"

printf '@bad_fastq\nACGT\n+\n!!!\n' >"$work_dir/invalid_quality.fq"
expect_failure invalid_quality python3 "$validator" --input "$work_dir/invalid_quality.fq"

printf '>too_long\n' >"$work_dir/too_long.fa"
for _ in $(seq 1 255); do
    printf 'A' >>"$work_dir/too_long.fa"
done
printf '\n' >>"$work_dir/too_long.fa"
expect_failure too_long python3 "$validator" --input "$work_dir/too_long.fa"

printf '>max_length_1\n' >"$work_dir/max_length.fa"
for _ in $(seq 1 254); do
    printf 'A' >>"$work_dir/max_length.fa"
done
printf '\n>max_length_2\n' >>"$work_dir/max_length.fa"
for _ in $(seq 1 254); do
    printf 'C' >>"$work_dir/max_length.fa"
done
printf '\n' >>"$work_dir/max_length.fa"
"$runner" "$work_dir/max_length.fa" "$work_dir/max_length"

if [[ -n "${EBWT2INDEL_BIN:-}" ]]; then
    "$EBWT2INDEL_BIN" -1 "$work_dir/plain.ebwt" -o "$work_dir/plain.snp" -t 36 \
        >"$work_dir/ebwt2indel_plain.log" 2>&1
    "$EBWT2INDEL_BIN" -1 "$work_dir/valid.rl_bwt" -o "$work_dir/rle.snp" -t 36 \
        >"$work_dir/ebwt2indel_rle.log" 2>&1
    "$EBWT2INDEL_BIN" -1 "$work_dir/packed.pck_bwt" -o "$work_dir/packed.snp" -t 36 \
        >"$work_dir/ebwt2indel_packed.log" 2>&1
    cmp "$work_dir/plain.snp" "$work_dir/rle.snp"
    cmp "$work_dir/plain.snp" "$work_dir/packed.snp"
    expect_failure ebwt2indel_truncated "$EBWT2INDEL_BIN" \
        -1 "$work_dir/truncated.rl_bwt" -o "$work_dir/truncated.snp" -t 36
    expect_failure ebwt2indel_checksum "$EBWT2INDEL_BIN" \
        -1 "$work_dir/bad_checksum.rl_bwt" -o "$work_dir/bad_checksum.snp" -t 36
    expect_failure ebwt2indel_bad_symbol "$EBWT2INDEL_BIN" \
        -1 "$work_dir/bad_symbol.pck_bwt" -o "$work_dir/bad_symbol.snp" -t 36

    python3 "$repo_dir/test/generate_variant_fixture.py" \
        --output-dir "$work_dir/variant_fixture"
    truth="$work_dir/variant_fixture/truth.tsv"
    test "$(awk -F '\t' 'NR > 1 && $2 == "snp" { n++ } END { print n + 0 }' "$truth")" -eq 1
    test "$(awk -F '\t' 'NR > 1 && $2 == "ins" { n++ } END { print n + 0 }' "$truth")" -eq 4
    test "$(awk -F '\t' 'NR > 1 && $2 == "del" { n++ } END { print n + 0 }' "$truth")" -eq 4
    test "$(awk -F '\t' 'NR > 1 { seen[$4] = 1 } END { print length(seen) }' "$truth")" -eq 3
    test "$(awk -F '\t' 'NR > 1 { seen[$5] = 1 } END { print length(seen) }' "$truth")" -eq 3

    for collection in reference alternate; do
        input="$work_dir/variant_fixture/$collection.fa"
        "$runner" "$input" "$work_dir/${collection}_rle"
        "$runner" "$input" "$work_dir/${collection}_plain" plain
        "$runner" "$input" "$work_dir/${collection}_packed" packed
        PYTHONPATH="$repo_dir/scripts" python3 - "$work_dir" "$collection" <<'PY'
from pathlib import Path
import sys
from bwt_format import read_all

root = Path(sys.argv[1])
name = sys.argv[2]
plain = read_all(root / f"{name}_plain.ebwt", "$", 64 * 1024 * 1024)
assert read_all(root / f"{name}_rle.rl_bwt", "$", 64 * 1024 * 1024) == plain
assert read_all(root / f"{name}_packed.pck_bwt", "$", 64 * 1024 * 1024) == plain
PY
    done

    "$EBWT2INDEL_BIN" -1 "$work_dir/reference_plain.ebwt" \
        -2 "$work_dir/alternate_plain.ebwt" -o "$work_dir/variants_plain.snp" -t 36 \
        >"$work_dir/variants_plain.log" 2>&1
    "$EBWT2INDEL_BIN" -1 "$work_dir/reference_rle.rl_bwt" \
        -2 "$work_dir/alternate_rle.rl_bwt" -o "$work_dir/variants_rle.snp" -t 36 \
        >"$work_dir/variants_rle.log" 2>&1
    "$EBWT2INDEL_BIN" -1 "$work_dir/reference_packed.pck_bwt" \
        -2 "$work_dir/alternate_packed.pck_bwt" -o "$work_dir/variants_packed.snp" -t 36 \
        >"$work_dir/variants_packed.log" 2>&1
    python3 - "$truth" "$work_dir/variants_plain.snp" <<'PY'
from pathlib import Path
import csv
import re
import sys

truth_path = Path(sys.argv[1])
calls_path = Path(sys.argv[2])
called = set(re.findall(r"type:_(?:SNP|INDEL)_event:([ACGT]*/[ACGT]*)", calls_path.read_text()))
missing = []
with truth_path.open(newline="") as handle:
    for event in csv.DictReader(handle, delimiter="\t"):
        ref = "" if event["ref"] == "-" else event["ref"]
        alt = "" if event["alt"] == "-" else event["alt"]
        allele = alt if not ref else ref
        rotations = {allele[index:] + allele[:index] for index in range(len(allele))}
        expected = {f"{ref}/{alt}", f"{alt}/{ref}"}
        if event["type"] != "snp":
            expected |= {f"/{rotation}" for rotation in rotations}
            expected |= {f"{rotation}/" for rotation in rotations}
        if called.isdisjoint(expected):
            missing.append(event["event"])
if missing:
    raise SystemExit("known variants not recovered: " + ", ".join(missing))
PY
    cmp "$work_dir/variants_plain.snp" "$work_dir/variants_rle.snp"
    cmp "$work_dir/variants_plain.snp" "$work_dir/variants_packed.snp"
fi

echo "PASS: ebwt2InDel BCR path"
