import pandas as pd
import sys

# Read input data and optional window size from command-line arguments
input_data = sys.argv[1]
window = int(sys.argv[2]) if len(sys.argv) > 2 else None
output_file = sys.argv[3]
chromconversion = sys.argv[4]


# Read the data into a DataFrame
df = pd.read_csv(input_data)

# Compute start and stop positions based on window size
df['start'] = df['position'] - (window // 2) if window else df['position'] - 1
df['stop'] = df['position'] + (window // 2) if window else df['position']

# Create a BED file DataFrame
bed_df = df[['chromo', 'start', 'stop']]


# Read chromosome conversion file
file2 = pd.read_csv(chromconversion, sep=',', header=None)

# Merge to replace values
bed_df[0] = bed_df[0].map(file2.set_index(0)[1])

# Save to a BED file
bed_df.to_csv(output_file, sep='\t', header=False, index=False)

print("BED file created successfully.")