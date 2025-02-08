import pandas as pd
import sys

# Read input data and optional window size from command-line arguments
input_data = sys.argv[1]
window = int(sys.argv[2]) if len(sys.argv) > 2 else None
output_file = sys.argv[3]

# Read the data into a DataFrame
df = pd.read_csv(input_data)

# Compute start and stop positions based on window size
df['start'] = df['position'] - (window // 2) if window else df['position'] - 1
df['stop'] = df['position'] + (window // 2) if window else df['position']

# Create a BED file DataFrame
bed_df = df[['chromo', 'start', 'stop']]

# Save to a BED file
bed_df.to_csv(output_file, sep='\t', header=False, index=False)

print("BED file created successfully.")