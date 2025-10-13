#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import tskit


def getAverageTimeOfRoot(ts):
    root_times = []
    for tree in ts.trees():
        if tree.num_roots != 1:
            continue
        root_node = tree.root
        time_of_root = tree.time(root_node)
        root_times.append(time_of_root)
    return root_times

def summarize_tree_sequence(ts, label, window_length=0.25e4):
    windows = np.arange(0, ts.sequence_length, window_length)
    windows = np.append(windows, ts.sequence_length)

    starts = []
    ends = []
    num_vars = []
    average_allele_age = []
    avg_span = []
    avg_root_time = []
    number_of_subtrees = []
    avg_total_branch_length = []

    for i in range(0, len(windows) - 1, 2):
        start = windows[i]
        end = windows[i + 1] if i + 1 < len(windows) else windows[-1]
        starts.append(start)
        ends.append(end)
        tree_subset = ts.keep_intervals(
            np.array([[start, end]]),
            simplify=True,
            record_provenance=True
        ).trim()

        interval_sizes = [tree.interval.right - tree.interval.left for tree in tree_subset.trees() if tree.num_roots == 1] #get size of interval, using trim it will be based on window coordinates
        
        num_vars.append(len([x for x in tree_subset.sites()])) #number of variants
        average_allele_age.append(np.average([var.site.mutations[0].time for var in tree_subset.variants()])) #average allele age

        total_branch_lengths = [tree.total_branch_length for tree in tree_subset.trees() if tree.num_roots == 1] #get total branch length
        avg_total_branch_length.append(np.average(total_branch_lengths, weights=interval_sizes)) #average of total branch length
        
        avg_root_time.append(np.average(getAverageTimeOfRoot(tree_subset), weights = interval_sizes)) #average time to root
        
        breakpoints = [x for x in tree_subset.breakpoints()] #number of breakpoints
        
        number_of_subtrees.append(len(breakpoints) - 1) #number of subtrees
        #average size of a local tree
        if len(breakpoints) > 2:  # at least one internal breakpoint
            avg_span.append(np.mean(np.diff(breakpoints[1:-1])))
        else:
            avg_span.append(np.nan)  # or 0, depending on what makes sense 

    df = pd.DataFrame({
        "start" : starts,
        "end" : ends,
        "num_vars": num_vars,
        "avg_allele_age": average_allele_age,
        "avg_span": avg_span,
        "avg_root_time": avg_root_time,
        "avg_total_branch_length" : avg_total_branch_length,
        "num_subtrees": number_of_subtrees,
        "label": label
    })
    return df


def main():
    parser = argparse.ArgumentParser(description="Summarize tree sequence windows.")
    parser.add_argument("tree_file", help="Input tree sequence file (.trees)")
    parser.add_argument("--label", required=True, help="Label for tree sequence")
    parser.add_argument("--window_length", type=float, default=0.25e4,
                        help="Window length (default: 0.25e4)")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")

    args = parser.parse_args()

    ts = tskit.load(args.tree_file).trim().simplify()
    df = summarize_tree_sequence(ts, args.label, args.window_length)

    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, f"{args.label}.summary.csv")

    df.to_csv(out_file, index=False)
    print(f"Summary written to {out_file}")


if __name__ == "__main__":
    main()
