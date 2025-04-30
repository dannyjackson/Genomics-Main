#!/bin/bash
#
# avg_depth_with_progress.sh
#
# This script computes the average depth across specified genome windows.
# It takes as input:
#   - A depth file with columns: chromosome, position, depth
#   - A windows file with columns: chromosome, starting position, ending position
#
# The final output has five columns:
#   chromosome, window start, window end, average depth, number of sites (denom)
#
# In addition, the script periodically prints progress updates (the total
# number of depth sites processed) to standard output.
#
# Usage:
#   avg_depth_with_progress.sh -d <DEPTH_FILE> -w <WIN_FILE> -a <AVG_OUTPUT_FILE> -t <THREADS>
#
# Example sbatch call:
#   sbatch --account=mcnew --job-name=fstavgdepthfastparallel \
#          --partition=standard --mail-type=ALL \
#          --output=slurm_output/output.%j \
#          --nodes=1 --ntasks-per-node=94 --time=24:00:00 \
#          avg_depth_with_progress.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVG_OUTPUT_FILE} -t 94

module load parallel

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -d <DEPTH_FILE> -w <WIN_FILE> -a <AVG_OUTPUT_FILE> -t <THREADS>

This script computes the average of depth values across genomic windows.
Required arguments:
  -d   Path to the depth file (chromosome, position, depth)
  -w   Path to the windows file (chromosome, start, end)
  -a   Path to the output file (chromosome, start, end, average depth, count)
  -t   Number of threads to use

The script will send periodic updates (sites processed) to standard output.
EOF
    exit 1
}

# Parse command-line arguments
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

# Create temporary working directory (will be removed on exit)
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Create a progress log file in the temporary directory.
PROGRESS_LOG="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/progress.log"
touch "$PROGRESS_LOG"

# Get the total number of sites (lines) in the depth file for reference.
total_sites=$(wc -l < "$DEPTH_FILE")
echo "Total number of depth sites: $total_sites"

echo "Splitting depth file into chunks..."
# Split the depth file into chunks of 100,000 lines (adjust as needed)
split -l 100000 "$DEPTH_FILE" "$TEMP_DIR/depth_chunk_"

# Function to process a single depth chunk.
# It runs the AWK command to compare depth entries to windows and then writes results.
# When finished it also writes the count of lines processed (sites) to the shared progress log.
process_chunk() {
    local chunk_file=$1
    local win_file=$2
    local out_dir=$3
    local out_file
    out_file=$(mktemp "$out_dir/chunk_result.XXXXXX")
    
    awk -v win_file="$win_file" -v out_file="$out_file" '
    BEGIN {
        OFS="\t";
        FS = "[ \t]+";
        # Load windows into arrays keyed by chromosome
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
        # Use whitespace/tab as field separator for the depth file as well.
        FS = "[ \t]+";
    }
    {
        chr = $1;
        pos = $2 + 0;
        depth = $3 + 0;
        if (!(chr in win_count)) next;
        for (i = 1; i <= win_count[chr]; i++) {
            # Because windows are sorted by start, if pos is less than the current window start,
            # no later window for that chromosome can include the position.
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
    
    # Append the number of sites processed in this chunk (i.e. the number of lines in the chunk)
    sites_in_chunk=$(wc -l < "$chunk_file")
    echo "$sites_in_chunk" >> "$PROGRESS_LOG"
}
export -f process_chunk

# Start a background progress reporter that prints updates every 10 seconds.
(
  while true; do
    # Use a flag file to determine when to exit this loop.
    if [ -f "$TEMP_DIR/parallel_done.flag" ]; then
      break
    fi
    if [ -s "$PROGRESS_LOG" ]; then
        processed=$(awk '{s+=$1} END {print s}' "$PROGRESS_LOG")
    else
        processed=0
    fi
    percent=$(awk -v proc="$processed" -v tot="$total_sites" 'BEGIN { 
         if (tot==0) { print "0.00" } else { printf "%.2f", (proc/tot)*100 } }')
    echo "Progress: Processed $processed of $total_sites sites ($percent% complete)"
    sleep 10
  done
) &
progress_pid=$!

echo "Processing chunks in parallel..."
find "$TEMP_DIR" -type f -name 'depth_chunk_*' | parallel -j "$THREADS" process_chunk {} "$WIN_FILE" "$TEMP_DIR"

# Signal the progress reporter to exit.
touch "$TEMP_DIR/parallel_done.flag"
wait $progress_pid

# Final update: sum the total sites processed.
final_processed=$(awk '{s+=$1} END {print s}' "$PROGRESS_LOG")
echo "Final progress: Processed $final_processed of $total_sites sites."

echo "Merging intermediate results..."
# Merge all per-chunk results. For windows processed in more than one chunk,
# their depths (sums) and counts are summed and the final average computed.
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
