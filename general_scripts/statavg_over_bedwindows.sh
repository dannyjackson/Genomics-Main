#!/bin/bash

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <parameter_file> -i <input_file> -n <population_name> -m <metric> -w <window_size>

This script creates a gene list given a list of sites in a genome and a gff file.

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -d  Path to a file containing a statistic (e.g. average depth, mapability, etc) at each site over a reference genome
    -w  Path to a bed file containing regions corresponding to the output of a windowed analysis (e.g. fst, RAiSD, etc)
    -o  Path to output file with window assigned to each site (each site can match multiple windows, as is the case in windowed analyses)
    -a  Path to output file for the average of the statistic over each window
EOF
    exit 1
fi


# Parse command-line arguments
while getopts "d:w:o:a:" option; do
    case "${option}" in
        d) STAT_FILE=${OPTARG} ;;
        w) WIN_FILE=${OPTARG} ;;
        o) OUTPUT_FILE=${OPTARG} ;;
        a) AVGSTAT_FILE=${OPTARG} ;;
        *) echo "Error: Invalid option '-${OPTARG}'" >&2; exit 1 ;;
    esac
done


# Load parallel (if available, or use another parallelization tool)
module load parallel


# Check if required arguments are passed
if [[ -z "$STAT_FILE" || -z "$WIN_FILE" || -z "$OUTPUT_FILE" || -z "$AVGSTAT_FILE" ]]; then
    echo "Usage: $0 <stat_file> <win_file> <output_file> <avgstat_file>"
    exit 1
fi

# Temporary directory to store chunked files
TEMP_DIR=$(mktemp -d)

echo ${TEMP_DIR}

trap "rm -rf $TEMP_DIR" EXIT

echo 'splitting stat file'

# Split the stat file into chunks of 100,000 lines each
split -l 100000 "$STAT_FILE" "$TEMP_DIR/stat_chunk_"

echo 'done splitting stat file'

# Function to process each chunk
process_chunk() {
    local chunk_file=$1
    local win_file=$2
    local output_file=$3
    
    awk -v win_file="$win_file" -v output_file="$output_file" '
    BEGIN {
        FS = "\t"
        print "Reading WIN_FILE..."
        
        # Read window file and organize by chromosome
        while ((getline line < win_file) > 0) {
            split(line, fields, FS)
            chr = fields[1]
            start = fields[2] + 0
            end = fields[3] + 0
            label = chr "_" start "_" end

            # Store windows by chromosome
            win_count[chr]++
            i = win_count[chr]
            win_start[chr, i] = start
            win_end[chr, i] = end
            win_label[chr, i] = label
        }
        close(win_file)

        print "Finished reading and indexing WIN_FILE\n"

        # Sort windows by start position for each chromosome
        print "Sorting windows by chromosome..."
        for (chr in win_count) {
            n = win_count[chr]
            # Bubble sort windows by start position
            for (i = 1; i < n; i++) {
                for (j = i + 1; j <= n; j++) {
                    if (win_start[chr, i] > win_start[chr, j]) {
                        # Swap the windows if out of order
                        tmp_start = win_start[chr, i]
                        tmp_end = win_end[chr, i]
                        tmp_label = win_label[chr, i]
                        win_start[chr, i] = win_start[chr, j]
                        win_end[chr, i] = win_end[chr, j]
                        win_label[chr, i] = win_label[chr, j]
                        win_start[chr, j] = tmp_start
                        win_end[chr, j] = tmp_end
                        win_label[chr, j] = tmp_label
                    }
                }
            }
        }

        print "Finished sorting windows by chromosome\n"
        FS = " "
    }
    {
        chr = $1
        pos = $2 + 0
        print "STAT_FILE line: chr=" chr ", pos=" pos

        # Skip stat entries with no matching chromosome in the window file
        if (!(chr in win_count)) next

        # Check all windows for the given chromosome
        for (i = 1; i <= win_count[chr]; i++) {
            wstart = win_start[chr, i]
            wend = win_end[chr, i]
            label = win_label[chr, i]

            # Stop early if position is smaller than the start of the current window (because sorted)
            if (pos < wstart) break

            # Match the position with the window range
            if (pos >= wstart && pos <= wend) {
                print "  Match found with window " label
                print $0, label >> output_file
            }
        }
    }
    ' "$chunk_file"
}

# Export the function to be used by parallel
export -f process_chunk

# Process all chunks in parallel
echo 'processing chunks in parallel'

find "$TEMP_DIR" -name 'stat_chunk_*' | parallel -j 16 process_chunk {} "$WIN_FILE" "$OUTPUT_FILE"

# compute stat average per window
echo 'computing average per window'

awk '{sum[$4] += $3; count[$4]++} END {for (w in sum) print w, sum[w]/count[w]}' "$OUTPUT_FILE" > "$AVGSTAT_FILE"


# Clean up
rm -rf "$TEMP_DIR"
rm "$OUTPUT_FILE"

