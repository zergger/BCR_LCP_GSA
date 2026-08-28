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
python3 "$validator" --input "$fixture" --ebwt "$work_dir/valid.ebwt" \
    --terminator '$' --verify-read-multiset
actual_sha256=$(sha256sum "$work_dir/valid.ebwt" | awk '{print $1}')
if [[ "$actual_sha256" != "$expected_sha256" ]]; then
    echo "Unexpected eBWT SHA-256: $actual_sha256" >&2
    exit 1
fi
test ! -e "$work_dir/valid.len"
test ! -e "$work_dir/valid.info"

"$runner" "$fastq_fixture" "$work_dir/fastq"
python3 "$validator" --input "$fastq_fixture" --ebwt "$work_dir/fastq.ebwt" \
    --terminator '$' --verify-read-multiset
actual_fastq_sha256=$(sha256sum "$work_dir/fastq.ebwt" | awk '{print $1}')
if [[ "$actual_fastq_sha256" != "$expected_fastq_sha256" ]]; then
    echo "Unexpected FASTQ eBWT SHA-256: $actual_fastq_sha256" >&2
    exit 1
fi

expect_failure wrong_terminator \
    python3 "$validator" --input "$fixture" --ebwt "$work_dir/valid.ebwt" --terminator '#'
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
    "$EBWT2INDEL_BIN" -1 "$work_dir/valid.ebwt" -o "$work_dir/valid.snp" -t 36 \
        >"$work_dir/ebwt2indel.log" 2>&1
fi

echo "PASS: ebwt2InDel BCR path"
