#!/usr/bin/env python
"""Exploratory Data Analysis
"""

from lists import REFINED_DIC, CORE_DIC
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.externals import joblib
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from myreduce import Dists15, Dists15Randomized, Dists15ShuffleLig

import aff_2015
import json
import luigi
import numpy as np
import pandas as pd
import numpy as np
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


class Tokens(luigi.Task):
    binning_size = luigi.Parameter(default=5.0)

    def output(self):
        ofns = [
            "/ddnB/work/jaydy/working/pdbbind/refined.{}.csv".format(
                self.binning_size),
            "/ddnB/work/jaydy/working/pdbbind/core.{}.csv".format(
                self.binning_size)
        ]

        return [luigi.LocalTarget(ofn) for ofn in ofns]

    def getTokens(self, profiles):
        df = pd.DataFrame(profiles)

        lig_span = df.dists.map(lambda dists: max(dists) - min(dists)).mean()

        tokens = []
        for profile in profiles:
            # ignore water molecule or residues too far from the ligand
            if profile["residue"] == "HOH":
                pass
            elif min(profile["dists"]) > 13.08:
                tokens.append(profile["residue"])
            else:
                ends = np.arange(13.08, 13.08 + lig_span, self.binning_size)
                type_dists = sorted(
                    zip(profile["atom_types"], profile["dists"]),
                    key=lambda x: x[1])

                token = profile["residue"]
                for end in ends:
                    my_type_dists = [
                        t
                        for t in type_dists
                        if t[1] > end and t[1] < (end + self.binning_size)
                    ]
                    if len(my_type_dists) > 0:
                        token = token + "-" + my_type_dists[0][0]
                tokens.append(token)

        return ' '.join(tokens)

    def run(self):
        ifn = "/work/jaydy/working/PDBbind_refined07-core07.distances.json"
        with open(ifn) as ifs:
            dat = json.loads(ifs.read())

        tokens = [self.getTokens(p) for p in dat.values()]
        myid_tokens = dict(zip(dat.keys(), tokens))

        refined_dat = [(myid, myid_tokens[myid], REFINED_DIC[myid])
                       for myid in dat.keys() if myid in REFINED_DIC]
        core_dat = [(myid, myid_tokens[myid], CORE_DIC[myid])
                    for myid in dat.keys() if myid in CORE_DIC]

        ofns = [output.path for output in self.output()]
        cols = ['myid', 'tokens', 'ki']
        pd.DataFrame(refined_dat, columns=cols).to_csv(ofns[0])
        pd.DataFrame(core_dat, columns=cols).to_csv(ofns[1])


class Tokens15(Tokens):
    def output(self):
        ofns = [
            "/ddnB/work/jaydy/working/pdbbind/refined.15.{}.csv".format(
                self.binning_size),
            "/ddnB/work/jaydy/working/pdbbind/core.15.{}.csv".format(
                self.binning_size)
        ]

        return [luigi.LocalTarget(ofn) for ofn in ofns]

    def requires(self):
        return Dists15()

    def run(self):
        ifn = self.requires().output().path

        with open(ifn) as ifs:
            dat = json.loads(ifs.read())

        tokens = [self.getTokens(p) for p in dat.values()]
        myid_tokens = dict(zip(dat.keys(), tokens))

        refined_dat = [(myid, myid_tokens[myid], aff_2015.refined_dat[myid])
                       for myid in dat.keys() if myid in aff_2015.refined_dat]
        core_dat = [(myid, myid_tokens[myid], aff_2015.core_dat[myid])
                    for myid in dat.keys() if myid in aff_2015.core_dat]

        ofns = [output.path for output in self.output()]
        cols = ['myid', 'tokens', 'ki']
        pd.DataFrame(refined_dat, columns=cols).to_csv(ofns[0])
        pd.DataFrame(core_dat, columns=cols).to_csv(ofns[1])


class Tokens15Randomized(Tokens):
    def output(self):
        ofns = [
            "/ddnB/work/jaydy/working/pdbbind/refined.15.rnd.{}.ki.csv".format(
                self.binning_size),
            "/ddnB/work/jaydy/working/pdbbind/refined.15.rnd.{}.kd.csv".format(
                self.binning_size),
        ]
        return [luigi.LocalTarget(ofn) for ofn in ofns]

    def requires(self):
        return Dists15Randomized()

    def run(self):
        ifn = self.requires().output().path

        with open(ifn) as ifs:
            dat = json.loads(ifs.read())

        tokens = [self.getTokens(p) for p in dat.values()]
        myid_tokens = dict(zip(dat.keys(), tokens))

        kds = set(aff_2015.Kds)
        kis = set(aff_2015.Kis)

        kds_dat = [(myid, myid_tokens[myid], aff_2015.refined_dat[myid])
                   for myid in dat.keys() if myid in kds]
        kis_dat = [(myid, myid_tokens[myid], aff_2015.refined_dat[myid])
                   for myid in dat.keys() if myid in kis]

        ofns = [output.path for output in self.output()]
        cols = ['myid', 'tokens', 'ki']
        pd.DataFrame(kis_dat, columns=cols).to_csv(ofns[0])
        pd.DataFrame(kds_dat, columns=cols).to_csv(ofns[1])


class Tokens15ShuffleLig(Tokens):
    def output(self):
        ofns = [
            "/ddnB/work/jaydy/working/pdbbind/refined.15.shuffled.{}.ki.csv".format(
                self.binning_size),
            "/ddnB/work/jaydy/working/pdbbind/refined.15.shuffled.{}.kd.csv".format(
                self.binning_size),
        ]
        return [luigi.LocalTarget(ofn) for ofn in ofns]

    def requires(self):
        return Dists15ShuffleLig()


class RF(luigi.Task):
    """try with the Random Forest algorithm
    """
    binning_size = luigi.Parameter(default=5.0)

    def requires(self):
        return Tokens(binning_size=self.binning_size)

    def read(self):
        task = self.requires()
        if not task.complete():
            raise Exception("{} not completed".format(task))

        refined_ifn, core_ifn = [_.path for _ in task.output()]

        refined_df = pd.read_csv(refined_ifn, index_col=0)
        core_df = pd.read_csv(core_ifn, index_col=0)

        return refined_df, core_df

    def split(self):
        refined_df, core_df = self.read()
        refined_df = refined_df[~refined_df['myid'].isin(core_df['myid'])]
        return refined_df, core_df

    def run(self):
        refined_df, core_df = self.split()
        tokens = refined_df.tokens.map(lambda x: x.split()).values
        unique_tokens = set([t for l in tokens for t in l])
        print("{} unique tokens".format(len(unique_tokens)))

        print("BinningSize MaxDf MinDf RMSE Corr")
        for max_df in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            for min_df in [0.0, 0.1, 0.2]:
                pipe_line = Pipeline([
                    ('tfidf', TfidfVectorizer(max_df=max_df,
                                              min_df=min_df,
                                              lowercase=False,
                                              token_pattern=r'(?u)\b\S+\b',
                                              analyzer='word')
                     ), ('model', RandomForestRegressor(n_estimators=50,
                                                        n_jobs=16))
                ])

                pipe_line.fit(refined_df['tokens'], refined_df['ki'])
                prediction = pipe_line.predict(core_df['tokens'])
                score = mean_squared_error(core_df['ki'], prediction)
                corr = pearsonr(core_df['ki'], prediction)[0]
                print("{} {} {} {} {}".format(self.binning_size, max_df,
                                              min_df, score, corr))

    def output(self):
        pass


def readDfireScores():
    """read the scores of Dfire on the PDBBind 2014 refined dataset
    """
    ifn = "../dat/pdbbind2014-dfire"
    df = pd.read_csv(ifn, index_col=0)
    return df


class RF15(RF):
    def requires(self):
        return Tokens15(binning_size=self.binning_size)


class RF15AgainstDfire(RF15):
    def output(self):
        ki_ofn = "/work/jaydy/working/pdbbind/PDBBind_refined_15/rf.ki.pkl"
        kd_ofn = "/work/jaydy/working/pdbbind/PDBBind_refined_15/rf.kd.pkl"
        return [luigi.LocalTarget(ki_ofn), luigi.LocalTarget(kd_ofn)]

    def run(self):
        refined_df, core_df = self.split()
        dfire_df = readDfireScores()
        merged = pd.merge(refined_df,
                          dfire_df,
                          left_on='myid',
                          right_on='pdbid')

        def train(df, target='ki'):
            y = df['ki']
            columns = df.columns.tolist()
            columns.remove('ki')
            X = df[columns]

            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.5)

            pipe_line = Pipeline([
                ('tfidf', TfidfVectorizer(max_df=1.0,
                                          min_df=0.0,
                                          lowercase=False,
                                          token_pattern=r'(?u)\b\S+\b',
                                          analyzer='word')
                 ), ('model', RandomForestRegressor(n_estimators=50,
                                                    n_jobs=1))
            ])

            pipe_line.fit(X_train['tokens'], y_train)
            ofn = "/work/jaydy/working/pdbbind/PDBBind_refined_15/rf.{}.pkl".format(
                target)
            joblib.dump(pipe_line, ofn)

            prediction = pipe_line.predict(X_test['tokens'])
            corr = pearsonr(y_test, prediction)

            dfire_corr = pearsonr(y_test, X_test['dfire'])
            uncorr_dfire_corr = pearsonr(y_test, X_test['uncor_dfire'])

            print("rf correlation: {}".format(corr))
            print("dfire correlation: {}".format(dfire_corr))
            print("uncorrelated-dfire correlation: {}".format(
                uncorr_dfire_corr))

        print("Ki")
        ki_df = merged[merged['myid'].isin(aff_2015.Kis)]
        train(ki_df, target='ki')

        print("Kd")
        kd_df = merged[merged['myid'].isin(aff_2015.Kds)]
        train(kd_df, target='kd')


def main():
    luigi.build(
        [
            # RF(binning_size=8.0),
            # RF(binning_size=7.0),
            # RF(binning_size=6.0),
            # RF(binning_size=5.0),

            # RF15(binning_size=8.0),
            # RF15(binning_size=7.0),
            # RF15(binning_size=6.0),
            # RF15(binning_size=5.0),
            RF15AgainstDfire(binning_size=7.0),
            # Tokens15Randomized(binning_size=7.0),
            # Tokens15ShuffleLig(binning_size=7.0),
        ],
        local_scheduler=True)


if __name__ == '__main__':
    main()
