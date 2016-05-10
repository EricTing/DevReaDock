#!/usr/bin/env python
"""Exploratory Data Analysis
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def distancesProfile():
    """Explore the distances profile
    """
    ifn = "/work/jaydy/working/PDBbind_refined07-core07.distances.json"
    with open(ifn) as ifs:
        dat = json.loads(ifs.read())

    # minimum distances
    max_dists = []
    min_dists = []
    mean_dists = []
    all_dists = []
    ligand_sizes_4_res = []
    for key, val in dat.iteritems():
        df = pd.DataFrame(val)
        min_dists.extend(df.dists.map(min))
        max_dists.extend(df.dists.map(max))
        mean_dists.extend(df.dists.map(np.mean))
        ligand_sizes_4_res.extend(df.dists.map(
            lambda dists: max(dists) - min(dists)))
        all_dists.extend(df.dists.values.tolist())

    all_dists = [d for l in all_dists for d in l]
    all_dists = pd.Series(all_dists)
    max_dists = pd.Series(max_dists)
    min_dists = pd.Series(min_dists)
    mean_dists = pd.Series(mean_dists)
    lig_sizes = pd.Series(ligand_sizes_4_res)

    plt.figure()
    all_dists.hist(bins=30)
    plt.xlabel("Distances [$\mathrm{\AA}$]")
    plt.savefig("/ddnB/work/jaydy/working/pdbbind/all_dists_hist.tiff")

    print("10% quatile of all distances: {}".format(all_dists.quantile(0.1)))

    plt.figure()
    lig_sizes.hist(bins=20)
    plt.xlabel("Ligand sizes [$\mathrm{\AA}$]")
    plt.savefig("/ddnB/work/jaydy/working/pdbbind/lig_sizes_hist.tiff")

    plt.figure()
    min_dists.hist(bins=30)
    plt.xlabel("Distances [$\mathrm{\AA}$]")
    plt.savefig("/ddnB/work/jaydy/working/pdbbind/min_dists.tiff")

    plt.figure()
    max_dists.hist(bins=30)
    plt.xlabel("Distances [$\mathrm{\AA}$]")
    plt.savefig("/ddnB/work/jaydy/working/pdbbind/max_dists.tiff")

    plt.figure()
    mean_dists.hist(bins=30)
    plt.xlabel("Distances [$\mathrm{\AA}$]")
    plt.savefig("/ddnB/work/jaydy/working/pdbbind/mean_dists.tiff")


def main():
    pass


if __name__ == '__main__':
    main()
