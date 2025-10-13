#!/usr/bin/env python3

import sys, argparse, os
import tskit
import numpy as np
import math, random
import pandas as pd

from utils import generateComparisonIntervals, load_local_trees

def main():
    parser = argparse.ArgumentParser(
        description="Create local trees in newick format to compare 2 sets of trees"
    )
    
    # Required positional args
    parser.add_argument("ts1", help="Path to .trees file for file 1")
    parser.add_argument("ts2", help="Path to .trees file for file 2")
    parser.add_argument("outfile", help="baseoutput name")

    # Optional args (have defaults)
    parser.add_argument("--output_dir", help="Output directory for newick files", default = ".")
    parser.add_argument("--ts1_offset", type=float, default=0,
                        help="Offset genomic position for ts1 (default: 0)")
    parser.add_argument("--ts2_offset", type=float, default=0,
                        help="Offset genomic position for ts2 (default: 0)")
    parser.add_argument("--file_splits", type=int, default=os.cpu_count() - 1,
                        help="Number of files to split into (default: N - 1 cores)")
    parser.add_argument("--max_trees", type=int, default=2000,
                        help="Maximum number of local trees to compare (default: 2000)")

    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output_dir), exist_ok=True)
    
    #read in tskit trees
    ts1 = tskit.load(args.ts1)
    ts2 = tskit.load(args.ts2)

    ts1_offset = args.ts1_offset
    ts2_offset = args.ts2_offset
    
    intervals = generateComparisonIntervals(ts1, ts2, ts1_offset, ts2_offset)[1:-1] #exclude the first and last interval
    if len(intervals) > args.max_trees:
        splits_of_4 = args.max_trees // 4
        
        split_intervals = np.array_split(np.array(intervals), len(intervals) // 4, axis = 0)
        
        sampled_intervals = random.sample(split_intervals, k = splits_of_4)
        intervals = np.vstack(sampled_intervals)
    print(f"Extracting {len(intervals)} local trees")
    
    ts1_local_trees = load_local_trees(
        ts=ts1,
        sampled_intervals=intervals,
        ts_offset=ts1_offset
    )

    ts2_local_trees = load_local_trees(
        ts=ts2,
        sampled_intervals=intervals,
        ts_offset=ts2_offset
    )
    
    merged_df = pd.merge(ts1_local_trees, ts2_local_trees, on = ["start", "end"])
    chunk_size = merged_df.shape[0] // args.file_splits
    
    for i in range(args.file_splits):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, merged_df.shape[0])
        if start_idx >= merged_df.shape[0]:
            break  # No more rows left
        chunk = merged_df.iloc[start_idx:end_idx]
        out_file = f"{args.output_dir}/{args.outfile}.newick_trees_part{i+1}.txt"

        # Write each row as: start end newick_ts1 newick_ts2
        chunk.to_csv(out_file, sep="\t", index=False, header = False)
        print(f"Wrote {chunk.shape[0]} rows to {out_file}")
    
if __name__ == "__main__":
    print("Running Main")
    main()
        
