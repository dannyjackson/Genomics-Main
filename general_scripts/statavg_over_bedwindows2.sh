#!/bin/bash

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -d <stat_file> -w <win_file> -o <site_output_file> -a <avg_output_file> -t <threads>

This script takes a bed file of windows and a file with a statistic per site and computes the average of the statistic across each window.

REQUIRED ARGUMENTS:
    -d  Path to a file containing a statistic (e.g. average depth) at each site
    -w  Path to a bed file containing windows
    -o  Path to output file for window-labeled sites (optional; can be left empty)
    -a  Path to output file for window averages
    -t  Number of threads to use
EOF
    exit 1
fi

# Parse command-line arguments
while getopts "d:w:o:a:t:" option; do
    case "${option}" in
        d) STAT_FILE=${OPTARG} ;;
        w) WIN_FILE=${OPTARG} ;;
        o) OUTPUT_FILE=${OPTARG} ;;  # Optional
        a) AVGSTAT_FILE=${OPTARG} ;;
        t) THREAD=${OPTARG} ;;
        *) echo "Error: Invalid option '-${OPTARG}'" >&2; exit 1 ;;
    esac
done

# Load parallel 
module load parallel


# Check required arguments
if [[ -z "$STAT_FILE" || -z "$WIN_FILE" || -z "$AVGSTAT_FILE" ]]; then
    echo "Missing required argument(s)."
    exit 1
fi

# Temporary directory for chunks
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Splitting stat file into chunks..."
split -l 100000 "$STAT_FILE" "$TEMP_DIR/stat_chunk_"

# Function to process each chunk
process_chunk() {
    local chunk_file=$1
    local win_file=$2
    local temp_dir=$3

    avg_output_file=$(mktemp "$temp_dir/avgstat.XXXXXX")
    
    awk -v win_file="$win_file" -v avg_output_file="$avg_output_file" '
    BEGIN {
        FS = "\t"
        # Load windows
        while ((getline line < win_file) > 0) {
            split(line, fields, FS)
            chr = fields[1]
            start = fields[2] + 0
            end = fields[3] + 0
            label = chr "_" start "_" end

            win_count[chr]++
            i = win_count[chr]
            win_start[chr, i] = start
            win_end[chr, i] = end
            win_label[chr, i] = label
        }
        close(win_file)
        FS = " "
    }
    {
        chr = $1
        pos = $2 + 0
        val = $3 + 0  # assume stat is in column 3

        if (!(chr in win_count)) next

        for (i = 1; i <= win_count[chr]; i++) {
            wstart = win_start[chr, i]
            wend = win_end[chr, i]
            label = win_label[chr, i]

            if (pos < wstart) break
            if (pos >= wstart && pos <= wend) {
                sum[label] += val
                count[label]++
            }
        }
    }
    END {
        for (label in sum) {
            avg = sum[label] / count[label]
            print label, avg, count[label] >> avg_output_file
        }
    }' "$chunk_file"
}

export -f process_chunk

# Run in parallel
echo "Processing chunks in parallel..."
find "$TEMP_DIR" -name 'stat_chunk_*' | parallel -j "$THREAD" process_chunk {} "$WIN_FILE" "$TEMP_DIR"

# Combine all average stat files
cat "$TEMP_DIR"/avgstat.* > "$AVGSTAT_FILE"

# Optional: collect site-level output
if [[ -n "$OUTPUT_FILE" ]]; then
    cat "$TEMP_DIR"/statavg.* > "$OUTPUT_FILE"
fi
