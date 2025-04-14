#!/bin/bash
#
# avg_depth_over_windows.sh
#
# This script computes the average depth across specified genome windows.
# It takes as input:
#   - A depth file with columns: chromosome, position, depth
#   - A windows file with columns: chromosome, starting position, ending position
#
# It outputs a file with five columns:
#   chromosome, starting position, ending position, average depth, number of sites
#
# The script splits the depth file into chunks, processes each chunk in parallel,
# and finally merges results from all chunks.
#
# Usage:
#   avg_depth_over_windows.sh -d <DEPTH_FILE> -w <WIN_FILE> -a <AVG_OUTPUT_FILE> -t <THREADS>
#
# Example sbatch call:
#   sbatch --account=mcnew --job-name=fstavgdepthfastparallel \
#          --partition=standard --mail-type=ALL \
#          --output=slurm_output/output.%j \
#          --nodes=1 --ntasks-per-node=94 --time=24:00:00 \
#          avg_depth_over_windows.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVG_OUTPUT_FILE} -t 94

set -euo pipefail

# Parse command-line arguments
usage() {
    cat <<EOF
Usage: $(basename "$0") -d <depth_file> -w <win_file> -a <avg_output_file> -t <threads>

This script computes the average of depth values across genomic windows.
Required arguments:
  -d   Path to the depth file (chromosome, position, depth)
  -w   Path to the windows file (chromosome, start, end)
  -a   Path to the output file (chromosome, start, end, average depth, count)
  -t   Number of threads to use

EOF
    exit 1
}

while getopts "d:w:a:t:" opt; do
    case "${opt}" in
        d) DEPTH_FILE=${OPTARG} ;;
        w) WIN_FILE=${OPTARG} ;;
        a) AVG_OUTPUT_FILE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        *) usage ;;
    esac
done

if [[ -z "${DEPTH_FILE:-}" || -z "${WIN_FILE:-}" || -z "${AVG_OUTPUT_FILE:-}" || -z "${THREADS:-}" ]]; then
    usage
fi

# Create a temporary working directory (will be removed on exit)
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Splitting depth file into chunks..."
# Split the depth file into chunks of 100,000 lines (adjust if needed)
split -l 100000 "$DEPTH_FILE" "$TEMP_DIR/depth_chunk_"

# Function to process a single depth chunk
process_chunk() {
    local chunk_file=$1
    local win_file=$2
    local out_dir=$3
    local out_file
    out_file=$(mktemp "$out_dir/chunk_result.XXXXXX")
    
    # Using AWK to:
    #  1. Read the windows file (assumes it is sorted by chromosome and start).
    #  2. For each depth record in the chunk, check (in order) which windows it falls into.
    #  3. Accumulate the sum of depths and count for each window.
    #
    # Output columns for each window found in the chunk:
    #   chromosome, start, end, sum_of_depths, count_of_sites
    awk -v win_file="$win_file" -v out_file="$out_file" '
    BEGIN {
        OFS="\t";
        FS = "[ \t]+";  # works if files are whitespace or tab separated
        # Read windows file and store by chromosome.
        while ((getline line < win_file) > 0) {
            split(line, fields, FS);
            chr = fields[1];
            start = fields[2] + 0;
            end = fields[3] + 0;
            win_count[chr]++; 
            i = win_count[chr];
            win_start[chr,i] = start;
            win_end[chr,i] = end;
        }
        close(win_file);
        # Set FS for depth file input (assumed whitespace separated)
        FS = "[ \t]+";
    }
    {
        chr = $1;
        pos = $2 + 0;
        depth = $3 + 0;
        if (!(chr in win_count)) next;
        for (i = 1; i <= win_count[chr]; i++) {
            # Because windows are sorted by start, if pos is less than the current window start, 
            # no later window can contain it.
            if (pos < win_start[chr,i])
                break;
            if (pos >= win_start[chr,i] && pos <= win_end[chr,i]) {
                key = chr "\t" win_start[chr,i] "\t" win_end[chr,i];
                sum[key] += depth;
                count[key]++;
            }
        }
    }
    END {
        for (key in sum) {
            print key, sum[key], count[key] > out_file;
        }
    }
    ' "$chunk_file"
}
export -f process_chunk

echo "Processing chunks in parallel..."
# Process each chunk in parallel
find "$TEMP_DIR" -type f -name 'depth_chunk_*' | parallel -j "$THREADS" process_chunk {} "$WIN_FILE" "$TEMP_DIR"

echo "Merging intermediate results..."
# Merge results from all chunks. In case the same window was processed in multiple chunks,
# we sum the depths and counts, then compute the final average depth.
awk 'BEGIN { OFS="\t" }
{
    key = $1 "\t" $2 "\t" $3;
    total_sum[key] += $4;
    total_count[key] += $5;
}
END {
    for (k in total_sum) {
        avg_depth = total_sum[k] / total_count[k];
        print k, avg_depth, total_count[k]
    }
}' "$TEMP_DIR"/chunk_result.* > "$AVG_OUTPUT_FILE"

echo "Done. Final output written to ${AVG_OUTPUT_FILE}"
