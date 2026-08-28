#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    echo "Usage: $0 INPUT_FASTA_OR_FASTQ OUTPUT_PREFIX [rle|packed|plain]" >&2
    echo "Defaults to OUTPUT_PREFIX.rl_bwt; use ebwt2InDel with -t 36." >&2
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
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
format=${3:-${BWT_FORMAT:-rle}}
case "$format" in
    rle) suffix=.rl_bwt ;;
    packed) suffix=.pck_bwt ;;
    plain) suffix=.ebwt ;;
    *)
        echo "Unsupported BWT format: $format (expected rle, packed, or plain)" >&2
        exit 2
        ;;
esac
final_ebwt="$output_prefix$suffix"
run_log="$output_prefix.bcr.log"
time_log="$output_prefix.bcr.time"
space_log="$output_prefix.bcr.space.tsv"

for target in "$final_ebwt" "$run_log" "$time_log" "$space_log"; do
    if [[ -e "$target" ]]; then
        echo "Refusing to overwrite existing output: $target" >&2
        exit 1
    fi
done

python3 "$validator" --input "$input"

stage=$(mktemp -d -- "$output_dir/.bcr-ebwt.XXXXXX")
completed=0
bcr_pid=
finish() {
    if [[ -n "$bcr_pid" ]] && kill -0 "$bcr_pid" 2>/dev/null; then
        kill "$bcr_pid" 2>/dev/null || true
        wait "$bcr_pid" 2>/dev/null || true
    fi
    if [[ $completed -eq 0 ]]; then
        echo "BCR stage retained for diagnosis: $stage" >&2
    fi
}
trap finish EXIT

printf 'elapsed_ms\tstage_bytes\n' >"$space_log"
start_ns=$(date +%s%N)
(
    cd -- "$stage"
    BCR_BWT_FORMAT="$format" /usr/bin/time -v -o "$time_log" "$bcr_bin" "$input" bcr
) >"$run_log" 2>&1 &
bcr_pid=$!
peak_bytes=0
while kill -0 "$bcr_pid" 2>/dev/null; do
    current_bytes=$(du -sb -- "$stage" 2>/dev/null | awk 'END { print $1 + 0 }' || true)
    [[ "$current_bytes" =~ ^[0-9]+$ ]] || current_bytes=0
    elapsed_ms=$(( ($(date +%s%N) - start_ns) / 1000000 ))
    printf '%s\t%s\n' "$elapsed_ms" "$current_bytes" >>"$space_log"
    if (( current_bytes > peak_bytes )); then
        peak_bytes=$current_bytes
    fi
    sleep 0.1
done
set +e
wait "$bcr_pid"
bcr_status=$?
bcr_pid=
set -e
current_bytes=$(du -sb -- "$stage" 2>/dev/null | awk 'END { print $1 + 0 }' || true)
[[ "$current_bytes" =~ ^[0-9]+$ ]] || current_bytes=0
elapsed_ms=$(( ($(date +%s%N) - start_ns) / 1000000 ))
printf '%s\t%s\n' "$elapsed_ms" "$current_bytes" >>"$space_log"
if (( current_bytes > peak_bytes )); then
    peak_bytes=$current_bytes
fi
printf 'peak\t%s\n' "$peak_bytes" >>"$space_log"
if (( bcr_status != 0 )); then
    echo "BCR failed with exit status $bcr_status" >&2
    exit "$bcr_status"
fi

stage_ebwt="$stage/bcr$suffix"
python3 "$validator" --input "$input" --ebwt "$stage_ebwt" --terminator '$'

unexpected=$(find "$stage" -mindepth 1 -maxdepth 1 \
    ! -name "bcr$suffix" ! -name 'bcr.info' ! -name 'bcr.len' -print -quit)
if [[ -n "$unexpected" ]]; then
    echo "Unexpected BCR stage artifact: $unexpected" >&2
    exit 1
fi

rm -f -- "$stage/bcr.info" "$stage/bcr.len"
mv -- "$stage_ebwt" "$final_ebwt"
rmdir -- "$stage"
completed=1
trap - EXIT

echo "Published $format eBWT: $final_ebwt"
echo "ebwt2InDel terminator option: -t 36"
