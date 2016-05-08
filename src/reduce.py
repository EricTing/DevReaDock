#!/usr/bin/env python

import pandas as pd
from collections import defaultdict
from translate import BuildTokens
from lists import REFINED_DIC


def main():
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


if __name__ == '__main__':
    main()
