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
        FS = " "
    }
    {
        chr = $1
        pos = $2 + 0

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
                print $0, label >> output_file
            }
        }
    }
    ' "$chunk_file"
}
