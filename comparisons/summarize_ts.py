#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import tskit

from utils import getTimeOfRootWindow, getMedianCoalescence

def summarize_tree_sequence(ts, label, window_length=0.25e4, offset = 0, verbose = 0):
    
    if window_length == 0: #if windows == 0, means use local tree as each start end
        windows = [ x for x in ts.breakpoints() ][1:-1] #remove first and last point 
    else:
        windows = np.arange(0, ts.sequence_length, window_length)
        windows = np.append(windows, ts.sequence_length)

    starts = []
    ends = []
    number_of_variants = []
    average_allele_age = []
    average_span = []
    average_root_time = []
    number_of_subtrees = []
    average_total_branch_length = []
    median_coalescences = []
    
    for i in range(0, len(windows) - 1):
        genomic_start = windows[i]
        genomic_end = windows[i + 1] if i + 1 < len(windows) else windows[-1]
        start = genomic_start - offset if genomic_start - offset > 0 else 0
        end = genomic_end - offset
        starts.append(start)
        ends.append(end)
        tree_subset = ts.keep_intervals(
            np.array([[genomic_start, genomic_end]]),
            simplify=True,
            record_provenance=True
        ).trim()

        interval_sizes = [tree.interval.right - tree.interval.left for tree in tree_subset.trees() if tree.num_roots == 1] #get size of interval, using trim it will be based on window coordinates

        num_vars = len([x for x in tree_subset.sites()]) #number of variants
        number_of_variants.append(num_vars)
        
        avg_allele_age = np.average([var.site.mutations[0].time for var in tree_subset.variants()]) #average allele age
        average_allele_age.append(avg_allele_age)

        avg_branch_lengths = np.average([tree.total_branch_length for tree in tree_subset.trees() if tree.num_roots == 1], weights=interval_sizes) #average of total branch length
        average_total_branch_length.append(avg_branch_lengths)
        
        avg_root_time = np.average(getTimeOfRootWindow(tree_subset), weights = interval_sizes) #average time to root
        average_root_time.append(avg_root_time)
        
        breakpoints = [x for x in tree_subset.breakpoints()] #number of breakpoints
        num_breakpoints = len(breakpoints)
        number_of_subtrees.append(num_breakpoints - 1) #number of subtrees
        
        #get median coalescence time
        median_coalescence = np.median([getMedianCoalescence(tree) for tree in tree_subset.trees() if tree.num_roots == 1]) #median coalescence
        median_coalescences.append(median_coalescence)
        
        span = 0
        if len(breakpoints) > 2:  # at least one internal breakpoint
            span = np.mean(np.diff(breakpoints))
        elif window_length == 0: #if down to local trees
            span = end - start
        elif window_length <= breakpoints[-1] - breakpoints[0]: #if window length is smaller than the local tree
            span = end - start
        else:
            raise ValueError("Uncaught Average Span")
        average_span.append(span)
        
        if verbose == 1:
            print(f"{start=}")
            print(f"{end=}")
            print(f"{genomic_start=}")
            print(f"{genomic_end=}")
            print(f"{num_vars=}")
            print(f"{num_breakpoints=}")
            print(f"{avg_allele_age=}")
            print(f"{avg_branch_lengths=}")
            print(f"{avg_root_time=}")
            print(f"{median_coalescence=}")
            print(f"{span=}")
            
    df = pd.DataFrame({
        "start" : starts,
        "end" : ends,
        "num_vars": number_of_variants,
        "avg_allele_age": average_allele_age,
        "avg_span": average_span,
        "avg_root_time": average_root_time,
        "avg_total_branch_length" : average_total_branch_length,
        "num_subtrees": number_of_subtrees,
        "median_coalescence_time": median_coalescences,
        "label": label
    })
    return df

def main():
    parser = argparse.ArgumentParser(description="Summarize tree sequence windows.")
    parser.add_argument("tree_file", help="Input tree sequence file (.trees)")
    parser.add_argument("--label", required=True, help="Label for tree sequence")
    parser.add_argument("--window_length", type=float, default=0.25e4,
                        help="Window length (default: 0.25e4), to summarize over all local trees set this to zero")
    parser.add_argument("--offset", type=float, default=0,
                        help="offest to subtract from coordinates (default: 0")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
    parser.add_argument("--verbose", type=int, default=0)

    args = parser.parse_args()

    ts = tskit.load(args.tree_file).trim().simplify()
    df = summarize_tree_sequence(ts, args.label, args.window_length, offset = args.offset, verbose = args.verbose)

    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, f"{args.label}.summary.csv")

    df.to_csv(out_file, index=False)
    print(f"Summary written to {out_file}")


if __name__ == "__main__":
    main()
