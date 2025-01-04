
import pandas as pd
raw_counts = pd.read_csv('deseq2_results.csv', index_col=0)
genes, samples = raw_counts.shape
header = f"#1.2\n{genes}\t{samples}\n"
description_col=['metadata unavailable' for i in range(len(raw_counts))]
raw_counts.insert(0, 'Description', description_col)
gct_data = raw_counts.reset_index()
with open('final_deq_data.gct', 'w') as f:
    f.write(header)
    gct_data.to_csv(f, sep=',', index=False)