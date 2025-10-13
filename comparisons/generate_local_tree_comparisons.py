#!/usr/bin/env python3
import argparse
import random
from itertools import combinations

import numpy as np
import pandas as pd
import tskit

from utils import getOffsetBreakpoints, generateComparisonIntervals, load_local_trees, compute_tree_metrics

def generate_tree_comparisons(intervals, ts1, ts2, ts1_offset = 0, ts2_offset = 0):
    # Precompute pairs of common samples
    common_samples = sorted(set(ts1.samples()).intersection(ts2.samples()))
    pairs = list(combinations(common_samples, 2))

    # Offset intervals
    offset_intervals = [
        (interval + ts1_offset, interval + ts2_offset, interval[0][0], interval[0][1])
        for interval in intervals
    ]

    results = []
    for interval_1, interval_2, interval_start, interval_end in offset_intervals:
        subset_ts1 = ts1.keep_intervals(interval_1, simplify=False, record_provenance=True)
        subset_ts2 = ts2.keep_intervals(interval_2, simplify=False, record_provenance=True)

        for sub_tree1 in subset_ts1.trees(sample_lists=True):
            for sub_tree2 in subset_ts2.trees(sample_lists=True):
                r = compute_tree_metrics(sub_tree1, sub_tree2, interval_start, interval_end, pairs)
                if r is not None:
                    results.append(r)

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Compare local trees between two tskit TreeSequences."
    )
    parser.add_argument("--ts1", required=True, help="Path to first .trees file")
    parser.add_argument("--ts2", required=True, help="Path to second .trees file")
    parser.add_argument("--out", required=True, help="Output CSV path")
    
    parser.add_argument("--start-position", type=float, help="Genomic start coordinates")
    parser.add_argument("--end-position", type=float, help="Genomic end coordinates")
    
    parser.add_argument("--comparison_name", required=True, help="Name of comparison")
    parser.add_argument("--ts1-tree", required=True, help="Type of tree ts1")
    parser.add_argument("--ts2-tree", required=True, help="Type of tree ts2")

    parser.add_argument("--max-trees", type=int, default=2000, help="Subsample cap")
    parser.add_argument("--drop-head-tail", type=int, default=1, help="Drop N intervals from head and tail (default: 2)")
    parser.add_argument("--seed", type=int, default=13, help="Random seed for subsampling")

    args = parser.parse_args()
    random.seed(args.seed)

    ts1 = tskit.load(args.ts1).trim()
    ts2 = tskit.load(args.ts2).trim()

    ts1_offset, ts1_breakpoints = getOffsetBreakpoints(ts1, args.ts1_tree)
    ts2_offset, ts2_breakpoints = getOffsetBreakpoints(ts2, args.ts2_tree)

    intervals = generateComparisonIntervals(ts1_intervals = ts1_breakpoints, ts2_intervals = ts2_breakpoints)
    # drop first/last N
    if args.drop_head_tail > 0 and len(intervals) > 2 * args.drop_head_tail:
        intervals = intervals[args.drop_head_tail : -args.drop_head_tail]

    if len(intervals) > args.max_trees:
        # sample evenly across genome by chunking into groups of 4
        splits_of_4 = max(args.max_trees // 4, 1)
        chunks = np.array_split(np.array(intervals, dtype=object), max(len(intervals) // 4, 1), axis=0)
        sampled = random.sample(list(chunks), k=min(splits_of_4, len(chunks)))
        intervals = np.vstack(sampled)

    print(f"Extracting {len(intervals)} local trees")

    df = generate_tree_comparisons(intervals = intervals, ts1 = ts1, ts2 = ts2,
                                   ts1_offset = ts1_offset, ts2_offset = ts2_offset)
    df["comparison_version"] = args.comparison_name
    df.to_csv(args.out, index=False)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
