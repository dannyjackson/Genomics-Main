import pandas as pd
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Compute average dxy in windows.")
parser.add_argument("--outdir", required=True, help="Output directory")
parser.add_argument("--pop1", required=True, help="First population name")
parser.add_argument("--pop2", required=True, help="Second population name")
parser.add_argument("--win", type=int, required=True, help="Window size")

args = parser.parse_args()

# Define output file prefix
output_prefix = f"{args.outdir}/analyses/dxy/{args.pop1}_{args.pop2}/{args.win}/{args.pop1}_{args.pop2}"

# Load the chromosome lengths

chroms_len = pd.read_csv(f"{args.outdir}/referencelists/autosomes_lengths.txt", sep="\t", header=None, names=["chromosome", "length"])

# Load the dxy data
dxy_data = pd.read_csv(f"{args.outdir}/analyses/dxy/{args.pop1}_{args.pop2}/snps/Dxy_persite_{args.pop1}_{args.pop2}.autosomes.txt", sep="\t", dtype={"chromo": str, "position": int, "dxy": float})

# Initialize an empty list to store the window data
window_data = []

# Define the window size
window_size = args.win

# Process each chromosome
for _, chrom_row in chroms_len.iterrows():
    chrom = chrom_row["chromosome"]
    chrom_length = chrom_row["length"]
    
    # Filter the dxy data for the current chromosome
    chrom_dxy = dxy_data[dxy_data["chromo"] == chrom]
    
    # Divide the chromosome into windows
    for start in range(0, chrom_length, window_size):
        end = min(start + window_size, chrom_length)
        mid = (start + end) // 2
        win_len = end - start
        # Filter dxy data for the current window
        window_sites = chrom_dxy[(chrom_dxy["position"] >= start) & (chrom_dxy["position"] < end)]
        
        # Calculate statistics
        total_sites = len(window_sites)
        average_dxy = window_sites["dxy"].mean() if total_sites > 3 else None
        
        # Append the result
        window_data.append([chrom, start, end, mid, win_len, total_sites, average_dxy])

# Convert the result into a DataFrame
output_df = pd.DataFrame(window_data, columns=["chromosome", "start", "end", "mid", "win_size", "total_sites", "dxy"])

# Ensure output directory exists
os.makedirs(f"{args.outdir}/analyses/dxy", exist_ok=True)

# Save to a new file
output_filename = f"{output_prefix}_average_dxy_{window_size}bp_windows.txt"
output_df.to_csv(output_filename, sep="\t", index=False, na_rep="NA")
