import sys
import pandas as pd
pd.set_option('compute.use_bottleneck', False)

# Read input data and optional window size from command-line arguments
input_data = sys.argv[1]
window = int(sys.argv[2]) if len(sys.argv) > 2 else None
output_file = sys.argv[3]
chromconversion = sys.argv[4]

# Read the data into a DataFrame (tab-separated, no quotes)
df = pd.read_csv(input_data, sep='\t', header=0, dtype=str)
df['position'] = df['midPos'].astype(int)  # Ensure 'position' column is integer

# Compute start and end positions based on window size
if 'pos' not in df.columns:
    df['pos'] = df['position'] - (window // 2) if window else df['position'] - 1

# Create a BED file DataFrame
bed_df = df[['chr', 'pos']].copy()
bed_df['pos'] = bed_df['pos'].astype(int)  # Ensure correct formatting


# Read chromosome conversion file (tab-separated, no quotes)
file2 = pd.read_csv(chromconversion, sep=',', header=None, dtype=str)

# Ensure file2 has at least two columns
if file2.shape[1] < 2:
    raise ValueError("Chromosome conversion file must have at least two columns.")

# Merge to replace values
mapping = dict(zip(file2[0], file2[1]))  # Create a dictionary for mapping
bed_df['chr'] = bed_df['chr'].map(mapping)

# Save to a BED file (tab-separated, no header, no index)
bed_df.to_csv(output_file, sep='\t', header=False, index=False)

print("BED file created successfully.")
