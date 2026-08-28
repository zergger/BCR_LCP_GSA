#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    echo "Usage: $0 INPUT_FASTA_OR_FASTQ OUTPUT_PREFIX" >&2
    echo "Builds OUTPUT_PREFIX.ebwt; use ebwt2InDel with -t 36." >&2
}

if [[ $# -ne 2 ]]; then
    usage
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(cd -- "$script_dir/.." && pwd -P)
validator="$script_dir/validate_ebwt.py"
bcr_bin=${BCR_BIN:-"$repo_dir/build/ebwt2indel/BCR_LCP_GSA"}

if [[ ! -x "$bcr_bin" ]]; then
    echo "BCR binary not found or not executable: $bcr_bin" >&2
    echo "Build it with: make -C $repo_dir ebwt2indel" >&2
    exit 1
fi
bcr_bin=$(realpath -e -- "$bcr_bin")

input=$(realpath -e -- "$1")
output_arg=$2
output_parent=$(dirname -- "$output_arg")
output_name=$(basename -- "$output_arg")
mkdir -p -- "$output_parent"
output_dir=$(cd -- "$output_parent" && pwd -P)
output_prefix="$output_dir/$output_name"
final_ebwt="$output_prefix.ebwt"
run_log="$output_prefix.bcr.log"
time_log="$output_prefix.bcr.time"

for target in "$final_ebwt" "$run_log" "$time_log"; do
    if [[ -e "$target" ]]; then
        echo "Refusing to overwrite existing output: $target" >&2
        exit 1
    fi
done

python3 "$validator" --input "$input"

stage=$(mktemp -d -- "$output_dir/.bcr-ebwt.XXXXXX")
completed=0
finish() {
    if [[ $completed -eq 0 ]]; then
        echo "BCR stage retained for diagnosis: $stage" >&2
    fi
}
trap finish EXIT

(
    cd -- "$stage"
    /usr/bin/time -v -o "$time_log" "$bcr_bin" "$input" bcr
) >"$run_log" 2>&1

python3 "$validator" --input "$input" --ebwt "$stage/bcr.ebwt" --terminator '$'

unexpected=$(find "$stage" -mindepth 1 -maxdepth 1 \
    ! -name 'bcr.ebwt' ! -name 'bcr.info' ! -name 'bcr.len' -print -quit)
if [[ -n "$unexpected" ]]; then
    echo "Unexpected BCR stage artifact: $unexpected" >&2
    exit 1
fi

rm -f -- "$stage/bcr.info" "$stage/bcr.len"
mv -- "$stage/bcr.ebwt" "$final_ebwt"
rmdir -- "$stage"
completed=1
trap - EXIT

echo "Published eBWT: $final_ebwt"
echo "ebwt2InDel terminator option: -t 36"
