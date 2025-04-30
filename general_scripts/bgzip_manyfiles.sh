Gotta bgzip some files

#!/bin/bash

module load samtools

# Loop through all .pestPG files in sub-subdirectories
for file in /xdisk/mcnew/dannyjackson/cardinals/analyses/thetas/*/*/*.pestPG; do
    if [ -f "$file" ]; then
        echo "Compressing: $file"
        bgzip "$file"
    fi
done

echo "Compression complete."
