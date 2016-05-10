#!/usr/bin/env python

import luigi
import json
import pandas as pd
from collections import defaultdict
from translate import BuildTokens, Distances
from lists import REFINED_DIC


class Dists(luigi.Task):
    def run(self):
        all_dat = {}
        for line in file("../dat/PDBbind_07.txt"):
            myid = line.rstrip()
            task = Distances(myid)
            if task.complete():
                with task.output().open('r') as ifs:
                    dat = json.loads(ifs.read())
                    all_dat[myid] = dat

        with self.output().open('w') as ofs:
            to_write = json.dumps(all_dat)
            ofs.write(to_write)

    def output(self):
        ofn = "/work/jaydy/working/PDBbind_refined07-core07.distances.json"
        return luigi.LocalTarget(ofn)


def refined_tokens():
    refined_ids = REFINED_DIC.keys()

    refined_dat = defaultdict(list)
    for myid in refined_ids:
        task = BuildTokens(myid)
        if task.complete():
            with task.output().open('r') as ifs:
                tokens = ifs.read()
                refined_dat['myid'].append(myid)
                refined_dat['tokens'].append(tokens)
                refined_dat['pKd/pKi'].append(REFINED_DIC[myid])

    refined_df = pd.DataFrame(refined_dat)
    refined_df.to_csv("../dat/PDBbind_refined07-core07.tokens.csv")


def main():
    luigi.build([Dists()], local_scheduler=True)


if __name__ == '__main__':
    main()
