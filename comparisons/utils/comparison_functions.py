# arg_utils/compare_functions.py

import numpy as np
import pandas as pd
import tskit
import dendropy as dp
from dendropy.calculate import treecompare


def getMedianTMRCA(tree, pairs, normalize_by_root = False):
    if normalize_by_root:
        root_time = getTimeOfRoot(tree)
        return np.median(np.array([tree.tmrca(x, y) / root_time for x, y in pairs]))
    else:
        return np.median(np.array([tree.tmrca(x, y) for x, y in pairs]))

def getTimeOfRoot(tree):
    root_node = tree.root
    return tree.time(root_node)
    
def getMedianTMRCAWindow(ts):
    median_tmrcas = []
    common_samples = sorted(set(tree.samples()).intersection(tree.samples()))
    pairs = list(combinations(common_samples, 2))
    
    for tree in ts.trees():
        if tree.num_roots != 1:
            continue
        median = getMedianTMRCA(ts, pairs)
        median_tmrcas.append(getMedianTMRCA(tree))
    return median_tmrcas

def getTimeOfRootWindow(ts):
    root_times = []
    for tree in ts.trees():
        if tree.num_roots != 1:
            continue
        root_times.append(getTimeOfRoot(tree))
    return root_times

def getTables(ts_table):
    '''
    input: e.g. ts.tables.nodes
    return pandas dataframe of tables
    '''
    return pd.DataFrame([row.asdict() for row in ts_table])

def getOffsetBreakpoints(ts, source):
    """
    Returns:
        offset (float): starting offset for coordinate system
        breakpoints (np.ndarray): array of breakpoint positions
    """
    all_breakpoints = np.unique(list(ts.breakpoints()))

    if source in ("simulation", "singer", "tsinfer"):
        return 0.0, all_breakpoints

    elif source in ("relate", "argweaver"):
        start_position = all_breakpoints[1]
        shifted = np.array([x - start_position for x in all_breakpoints if x - start_position >= 0])
        return start_position, np.unique(shifted)

    else:
        raise ValueError("Currently only coordinates of singer, tsinfer, ms_prime simulation, relate and argweaver are implemented")

def generateComparisonIntervals(ts1_intervals, ts2_intervals):
    """
    Generate aligned comparison intervals from two sets of breakpoints.

    Args:
        ts1_intervals (array-like): breakpoints from ts1 (shifted from 0 --> Seq length).
        ts2_intervals (array-like): breakpoints from ts2 (shifted from 0 --> Seq length).

    Returns:
        list of np.ndarray: each interval as shape (1, 2) array [[start, end]]
    """
    max_interval = min(ts1_intervals[-1], ts2_intervals[-1])
    all_intervals = np.sort(np.unique([x for x in ts1_intervals if x <= max_interval] + 
                                       [x for x in ts2_intervals if x <= max_interval]))
    starts = all_intervals[0:-1]
    ends = all_intervals[1:]
    intervals = [np.array([[start, end]]) for start, end in zip(starts, ends)]
    return intervals


def load_local_trees(ts, sampled_intervals, ts_offset):
    """
    Extract local trees and return as a pandas DataFrame.

    Parameters
    ----------
    ts : tskit.TreeSequence
        The tree sequence to pull trees from.
    sampled_intervals : list
        List of intervals [(start, end), ...].
    ts_offset : float
        Offset for ts intervals.
    """
    rows = []
    offset_intervals = [(interval - ts_offset, interval[0][0], interval[0][1])
                        for interval in sampled_intervals]

    for interval, interval_start, interval_end in offset_intervals:
        subset_ts = ts.keep_intervals(interval, simplify=False, record_provenance=False)
        for tree in subset_ts.trees():
            if tree.num_roots == 1:
                rows.append({
                    "start": interval_start,
                    "end": interval_end,
                    "newick": tree.as_newick()
                })

    return pd.DataFrame(rows, columns=["start", "end", "newick"])

def calculate_unweighted_rf(t1_newick, t2_newick):
    '''
    Use dendropy to calulate the Robinsons Fould distance between two newick strings
    '''
    shared_taxa = dp.TaxonNamespace()
    t1 = dp.Tree.get(data=t1_newick, schema="newick", taxon_namespace=shared_taxa)
    t2 = dp.Tree.get(data=t2_newick, schema="newick", taxon_namespace=shared_taxa)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return treecompare.unweighted_robinson_foulds_distance(t1, t2)


def getMedianCoalescence(tree, normalize = True):
    '''
    Given a local tree get median coalescence time (based on the time of each internal node)
    if normalize true the time is divided by time to root (so a value close to 1 means deep coalescences and 0 means shallow)
    '''
    if tree.num_roots != 1:
        return []
    root_time = getTimeOfRoot(tree)
    norm = 1
    if normalize:
        norm = root_time
    node_times = [ tree.time(node) for node in tree.nodes() ]
    normalized_node_times = [time / norm for time in node_times if time < root_time and time > 0]
    assert len(normalized_node_times) > 0, "node_times array in getMedianCoalescence is empty"
    return np.median(normalized_node_times)
    

def compute_tree_metrics(sub_tree1, sub_tree2, interval_start, interval_end, pairs):
    """
    Compute multiple topological and branch length aware similarity metrics between two local trees.

    This function compares two local trees extracted from an ARG (Ancestral Recombination Graph)
    over a specific genomic interval, quantifying how similar they are in topology, branch lengths,
    and pairwise coalescent properties.

    Parameters
    ----------
    sub_tree1 : tskit.Tree
        First local tree for the given interval.
    sub_tree2 : tskit.Tree
        Second local tree for the same interval.
    interval_start : float
        Start position of the genomic interval.
    interval_end : float
        End position of the genomic interval.
    pairs : list of tuple(int, int)
        List of sample ID pairs to compute pairwise metrics (e.g., TMRCA, path length).

    Returns
    -------
    dict or None
        Dictionary of tree-level and pairwise metrics for this interval, or None if
        either tree has multiple roots (i.e., invalid for direct comparison).

    Notes
    -----
    The computed metrics include:
      • **RF Distance (rf_distance):** Topological dissimilarity between two trees 
        based on shared bipartitions (unweighted Robinson–Foulds distance).
      • **KC Distances (kc_distance_branch_length / kc_distance_topology):**
        Kendall–Colijn distances that capture differences in branch lengths (λ=1) 
        or purely topology (λ=0).
      • **Euclidean Path and TMRCA distances:** Quantify how much pairwise distances
        and coalescent times differ between the two trees across sample pairs.
      • **R² metrics (path_distance_r2, pairwise_tmrca_r2):** Squared Pearson correlation
        coefficients measuring concordance of path lengths and TMRCA values.
    """
    
    # Skip if invalid roots
    if sub_tree1.num_roots != 1 or sub_tree2.num_roots != 1:
        return None

    # Tree-level metrics
    tree1_newick = sub_tree1.as_newick()
    tree2_newick = sub_tree2.as_newick()
    rf_distance = calculate_unweighted_rf(tree1_newick, tree2_newick)

    N = sub_tree1.num_samples()
    kc_bl = sub_tree1.kc_distance(sub_tree2, lambda_=1)
    kc_top = sub_tree1.kc_distance(sub_tree2, lambda_=0)

    # Pairwise metrics
    pairs_arr = np.asarray(pairs)
    # vectorized path length & TMRCA lookups aren’t available; list comps are fine here
    dist_vec_1 = np.array([sub_tree1.path_length(x, y) for x, y in pairs_arr])
    dist_vec_2 = np.array([sub_tree2.path_length(x, y) for x, y in pairs_arr])
    tmrca_vec_1 = np.array([sub_tree1.tmrca(x, y) for x, y in pairs_arr])
    tmrca_vec_2 = np.array([sub_tree2.tmrca(x, y) for x, y in pairs_arr])

    euclid_path = float(np.sqrt(np.sum((dist_vec_1 - dist_vec_2) ** 2)))
    euclid_tmrca = float(np.sqrt(np.sum((tmrca_vec_1 - tmrca_vec_2) ** 2)))

    # Pearson R^2 (guard small vectors)
    r2_path = (
        float(np.corrcoef(dist_vec_1, dist_vec_2)[0, 1] ** 2)
        if dist_vec_1.size > 1
        else np.nan
    )
    r2_tmrca = (
        float(np.corrcoef(tmrca_vec_1, tmrca_vec_2)[0, 1] ** 2)
        if tmrca_vec_1.size > 1
        else np.nan
    )

    return {
        "interval_start": interval_start,
        "interval_end": interval_end,
        "rf_distance": rf_distance,
        "kc_distance_branch_length": kc_bl,
        "kc_distance_topology": kc_top,
        "num_samples": N,
        "path_distance": euclid_path,
        "pairwise_tmrca": euclid_tmrca,
        "path_distance_r2": r2_path,
        "pairwise_tmrca_r2": r2_tmrca,
        "tree1_newick": tree1_newick,
        "tree2_newick": tree2_newick,
    }

##Specify exposed functions:
__all__ = [
    
    "generateComparisonIntervals",
    "load_local_trees",
    "calculate_unweighted_rf",
    "compute_tree_metrics",
    "getOffsetBreakpoints",
    "getTables",
    "getMedianTMRCA",
    "getTimeOfRoot",
    "getMedianTMRCAWindow",
    "getTimeOfRootWindow",
    "getMedianCoalescence"
]
