
# max win fst 

# option 2
library(data.table)

# Define window size
window <- 25000  # Change this to your desired window size
win <- window/2  # Change this to your desired window size

# Read files efficiently
file1 <- fread("/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/25000/pyrrurban_pyrrrural.25000.fst", sep = "\t", header = TRUE)
file2 <- fread("/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/1/pyrrurban_pyrrrural.1.fst", sep = "\t", header = TRUE)


# Rename columns in file1 and file2 to avoid conflicts
setnames(file1, c("region1", "chr1", "midPos1", "Nsites1", "fst1"))
setnames(file2, c("region2", "chr2", "midPos2", "Nsites2", "fst2"))

# Set keys for fast searching
setkey(file2, chr2, midPos2)

# Function to find the maximum fst within the window
get_max_fst <- function(chr_val, pos_val) {
  max(file2[chr2 == chr_val & midPos2 >= (pos_val - win) & midPos2 <= (pos_val + win), fst2], na.rm = TRUE)
}

# Apply the function row-wise using a vectorized approach
file1[, fst_max := get_max_fst(chr1, midPos1), by = .(chr1, midPos1)]

# Rename columns back for output consistency
setnames(file1, c("region1", "chr1", "midPos1", "Nsites1", "fst1", "fst_max"),
                  c("region", "chr", "midPos", "Nsites", "fst", "fst_max"))

# Save the result
fwrite(file1, "/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.25000.fst.max", sep = "\t")




#!/usr/bin/env Rscript

# max win fst 

# option 2
library(data.table)

# Define window size
window <- 25000  # Change this to your desired window size
win <- window/2  # Change this to your desired window size

# Read files efficiently
file1 <- fread("/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/25000/pyrrurban_pyrrrural.25000.fst", sep = "\t", header = TRUE)
file2 <- fread("/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/1/pyrrurban_pyrrrural.1.fst", sep = "\t", header = TRUE)


# Rename columns in file1 and file2 to avoid conflicts
setnames(file1, c("region1", "chr1", "midPos1", "Nsites1", "fst1"))
setnames(file2, c("region2", "chr2", "midPos2", "Nsites2", "fst2"))

# Set keys for fast searching
setkey(file2, chr2, midPos2)

# Function to find the maximum fst within the window
get_max_fst <- function(chr_val, pos_val) {
  max(file2[chr2 == chr_val & midPos2 >= (pos_val - win) & midPos2 <= (pos_val + win), fst2], na.rm = TRUE)
}

# Apply the function row-wise using a vectorized approach
file1[, fst_max := get_max_fst(chr1, midPos1), by = .(chr1, midPos1)]

# Rename columns back for output consistency
setnames(file1, c("region1", "chr1", "midPos1", "Nsites1", "fst1", "fst_max"),
                  c("region", "chr", "midPos", "Nsites", "fst", "fst_max"))

# Save the result
fwrite(file1, "/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.25000.fst.max", sep = "\t")

[dannyjackson@cpu37 pyrrurban_pyrrrural]$ cat fstmaxwin.test.sh 
#!/bin/sh

module load R

Rscript /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/fstmaxwin.test.R