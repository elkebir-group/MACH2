def process_colormap_file(colormap_file):
    """Takes as input a colormap file and returns a colormap (dictionary)."""
    colormap = {}
    with open(colormap_file, 'r') as f:
        for line in f:
            site, color = line.split()
            colormap[site] = color
    return colormap

def get_colormap(locations):
    """Takes as input a collection of locations and returns a colormap (dictionary)."""
    return {s: (i+1) for i, s in enumerate(locations)}

def bootstrap_based_location_labeling_from_ccfs(df_mut, tree, n_boot=1000, boot_confidence=0.9, threshold=0.05, boot_func=None):
    '''
    Initializes the MACH2 solver with a clone tree, optional primary location, and criteria ordering.
        
    :param df_mut: pandas dataframe with index as mutation IDs, one column corresponding to the clone, and other columns are for CCFs in each location.
    :param tree: mach2.tree.Tree object.
    :param n_boot: Number of bootstrap replicates (default: 1000).
    :param boot_confidence: Fraction of bootstraps confirming a clone to be present in a location (default: 0.9).
    :param threshold: If parent clone CCF > sum of children clones' CCFs + threshold, then parent is present in the location with children (default: 0.05).
    :param boot_func: Function to summarize the bootstrapped CCFs (default: np.mean).
    '''

    import numpy as np
    import pandas as pd
    from collections import defaultdict

    if boot_func is None:
        boot_func = np.mean

    def bootstrap_values(series, n_boot=1000, func=np.mean):
        bootstrapped = [
            func(np.random.choice(series, size=len(series), replace=True))
            for _ in range(n_boot)
        ]
        return bootstrapped


    def bootstrap_df(df_group, n_boot=1000, func=np.mean):
        bootstrapped = {}
        for c in df_group.columns:
            if c != 'clone':
                bootstrapped[c] = [bootstrap_values(df_group[c], n_boot, func)]
        return pd.DataFrame(bootstrapped)

    def get_sum_children(df_mut, tree, u, l, i):
        d = {}
        for v in tree.get_descendants(u):
            if v not in tree.leaves:
                d[v] = max(sum(df_mut.loc[int(w)][l][i] for w in tree.get_children(v)), df_mut.loc[int(v)][l][i])
            else:
                d[v] = df_mut.loc[int(v)][l][i]
        return sum(max(d[w] for w in [v, *tree.get_descendants(v)]) for v in tree.get_children(u))

    
    df_mut = df_mut.groupby('clone').apply(lambda s: bootstrap_df(s, n_boot=n_boot, func=boot_func))
    df_mut = df_mut.reset_index(level=0).set_index('clone')

    counter = {}
    for u in tree.nodes:
        if u not in tree.leaves:
            counter[u] = {}
            for l in df_mut.columns:
                counter[u][l] = []
                for i in range(n_boot):
                    if df_mut.loc[int(u),l][i] > get_sum_children(df_mut, tree, u, l, i) + threshold:
                        counter[u][l].append(i)
        else:
            counter[u] = {}
            for l in df_mut.columns:
                counter[u][l] = []
                for i in range(n_boot):
                    if df_mut.loc[int(u),l][i] > 0:
                        counter[u][l].append(i)

    labelmap = defaultdict(list)
    for u in counter:
        for l in counter[u]:
            if len(counter[u][l]) > n_boot * boot_confidence:
                labelmap[u].append(l[:-9])

    return labelmap