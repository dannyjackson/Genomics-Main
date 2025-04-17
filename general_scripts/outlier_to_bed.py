import sys
import numpy
import pandas as pd
pd.set_option('compute.use_bottleneck', False)

# Read input data and optional window size from command-line arguments
input_data = sys.argv[1]
window = int(sys.argv[2]) if len(sys.argv) > 2 else None
output_file = sys.argv[3]
chromconversion = sys.argv[4]


# Read the data into a DataFrame
df = pd.read_csv(input_data)

# Compute start and end positions based on window size
if 'start' not in df.columns:
    df['start'] = df['position'] - (window // 2) if window else df['position'] - 1
if 'end' not in df.columns:
    df['end'] = df['position'] + (window // 2) if window else df['position']


# Create a BED file DataFrame
bed_df = df[['chromo', 'start', 'end']].copy()  # Copy to avoid modifying df
bed_df['chromo'] = bed_df['chromo'].astype(str).str.strip()

# Read chromosome conversion file
file2 = pd.read_csv(chromconversion, sep=',', header=None)

# Ensure file2 has at least two columns
if file2.shape[1] < 2:
    raise ValueError("file2 does not have at least two columns.")

# Merge to replace values
mapping = dict(zip(file2[0], file2[1]))  # Create a dictionary for mapping
bed_df.loc[:, 'chromo'] = bed_df['chromo'].map(mapping)

# Save to a BED file
bed_df.to_csv(output_file, sep='\t', header=False, index=False)

print("BED file created successfully.")