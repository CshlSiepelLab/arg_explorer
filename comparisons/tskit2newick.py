#!/usr/bin/env python3

import argparse
import tskit

def main():
    # -----------------------------
    # Parse command-line arguments
    # -----------------------------
    parser = argparse.ArgumentParser(
        description="Export single-root trees from a tskit tree sequence to a Newick file."
    )
    parser.add_argument(
        "--tree", 
        required=True, 
        help="Path to the input .trees file (tskit tree sequence)."
    )
    parser.add_argument(
        "--out", 
        required=True, 
        help="Path to the output file where Newick strings will be written."
    )
    args = parser.parse_args()

    # -----------------------------
    # Load tree sequence
    # -----------------------------
    ts = tskit.load(args.tree)

    # -----------------------------
    # Write Newick trees to output
    # -----------------------------
    with open(args.out, "w") as out_f:
        for tree in ts.trees():
            if tree.num_roots != 1:
                continue
            newick_string = tree.as_newick()
            out_f.write("\t".join([str(tree.interval.left), str(tree.interval.right), newick_string]) + "\n")

    print(f"Exported Newick trees to: {args.out}")

if __name__ == "__main__":
    main()
