#!/usr/bin/env python

import pprint

core_ifn = "/work/jaydy/dat/pdbbind_v2015_docs/INDEX_core_data.2013"

core_dat = {}
for line in file(core_ifn):
    if not line.startswith('#'):
        tokens = line.split()
        myid = tokens[0]
        aff = tokens[3]
        core_dat[myid] = aff

print("core_dat =")
pprint.pprint(core_dat)

refined_ifn = "/work/jaydy/dat/pdbbind_v2015_docs/INDEX_refined_data.2015"

refined_dat = {}
for line in file(refined_ifn):
    if not line.startswith('#'):
        tokens = line.split()
        myid = tokens[0]
        aff = tokens[3]
        refined_dat[myid] = aff

print("refined_dat =")
pprint.pprint(refined_dat)

Kds = []
Kis = []
for line in file(refined_ifn):
    if not line.startswith('#'):
        tokens = line.split()
        myid = tokens[0]
        if 'Kd' in line:
            Kds.append(myid)
        elif 'Ki' in line:
            Kis.append(myid)

print("Kds = ")
pprint.pprint(Kds)

print("Kis = ")
pprint.pprint(Kis)
