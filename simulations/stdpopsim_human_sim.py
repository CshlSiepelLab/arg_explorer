#!/usr/bin/env python3

import argparse
import json
import sys

import stdpopsim
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Generic functions to simulate human chroms using stdpopsim.")
    p.add_argument("--model", required=True, help="Demographic model ID (e.g., OutOfAfrica_3G09).")
    
    p.add_argument("--species", default="HomSap", help="Species ID (default: HomSap for humans).")
    
    p.add_argument("--chrom", required=True, help="Chromosome/contig name (e.g., chr22).")
    
    p.add_argument("--left", type=float, default=None, help="Left coordinate (bp) of subregion (optional).")
    
    p.add_argument("--right", type=float, default=None, help="Right coordinate (bp) of subregion (optional).")
    
    p.add_argument("--samples", default=None,
                   help=("JSON dict of haploid sample counts per population ID, e.g. "'\'{"YRI": 5, "CHB": 5, "CEU": 0}\' '))
    
    p.add_argument("--seed", type=int, default=None,help="Random seed (optional).")
    
    p.add_argument("--out", required=True, help="Output .trees filename.")
    return p.parse_args()

def build_samples_list(model, samples_json):
    """
    Build stdpopsim SampleSet list with ploidy fixed at 2.

    samples_json: None or JSON string like {"YRI": 5, "CHB": 5, "CEU": 0}
    """
    pop_ids = [pop.name for pop in model.populations]
    samples_list = defaultdict(str)

    if samples_json is None:
        # Default: one diploid individual per model population
        for pop in model.populations:
            samples_list[pop] = 1
        return samples_list
    try:
        requested = json.loads(samples_json)
        if not isinstance(requested, dict):
            raise ValueError
    except Exception as e:
        raise ValueError(
            f"--samples must be a JSON object mapping population IDs to DIPLOID counts; got: {samples_json}"
        ) from e

    for pop_label, n_indiv in requested.items():
        if pop_label not in pop_ids:
            raise ValueError(
                f"Population '{pop_label}' not found in model. "
                f"Available populations: {pop_ids}"
            )
        if not isinstance(n_indiv, int) or n_indiv < 0:
            raise ValueError(f"Diploid count must be a non-negative integer for '{pop_label}'.")
        if n_indiv == 0:
            continue
        samples_list[pop_label] = n_indiv

    if len(samples_list.keys()) == 0:
        raise ValueError("No samples requested (all zero). Provide at least one non-zero population count.")
    return samples_list

def main():
    args = parse_args()
    # Load species and model
    species = stdpopsim.get_species(args.species)
    try:
        model = species.get_demographic_model(args.model)
    except Exception as e:
        available = [m.id for m in species.demographic_models]
        raise SystemExit(
            f"Unknown model '{args.model}' for species '{args.species}'. "
            f"Available: {available}"
        ) from e

    # Build samples
    try:
        samples = build_samples_list(model, args.samples)
    except ValueError as e:
        print(f"Sample specification error: {e}", file=sys.stderr)
        sys.exit(2)
    
    engine = stdpopsim.get_engine("msprime")
    
    contig = species.get_contig(args.chrom, mutation_rate=model.mutation_rate, left = args.left, right = args.right)
    
    print("mean recombination rate:", f"{contig.recombination_map.mean_rate:.3}")
    print("contig mutation rate:", contig.mutation_rate)
    print("model mutation rate:", model.mutation_rate)

    ts = engine.simulate(
        demographic_model=model,
        contig=contig,
        samples=samples,
        seed=args.seed
    )

    with open(f"{args.out}.vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file)
    ts.dump(f"{args.out}.trees")

if __name__ == "__main__":
    main()
