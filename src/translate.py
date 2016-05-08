#!/usr/bin/env python

from Bio.PDB import PDBParser
from scipy.spatial.distance import euclidean
from openbabel import OBTypeTable

import numpy as np
import os
import pybel
import luigi
import paths


def Structure2Seq(lig, prt):
    """
    Keyword Arguments:
    lig -- ligand
    prt -- protein
    """
    tokens = []
    for residue in prt.get_residues():
        token = tokenize(residue, lig)
        tokens.append(token)

    return ' '.join(tokens)


def residueCenter(residue):
    """
    Keyword Arguments:
    residue -- 
    """
    coords = [a.coord for a in residue.get_atom()]
    return np.mean(coords, axis=0)


def tokenize(residue, lig):
    """
    Keyword Arguments:
    residue -- BioPython protein residue
    lig     -- pybel lig with hydrogen removed
    """
    typetable = OBTypeTable()
    typetable.SetFromType('INT')
    typetable.SetToType('SYB')

    residue_center = residueCenter(residue)
    l = []
    for idx, atom in enumerate(lig.atoms):
        atom_type = typetable.Translate(atom.type)
        dist = euclidean(atom.coords, residue_center)
        l.append((idx, dist, atom_type))

    l = [t for t in l if t[1] < 4.5]
    l = sorted(l, key=lambda x: x[1])
    if len(l) == 0:
        return residue.get_resname()
    elif len(l) == 1:
        return "{}_{}".format(residue.get_resname(), l[0][2])
    elif len(l) > 1:
        return "{}_{}_{}".format(residue.get_resname(), l[0][2], l[1][2])


class BuildTokens(luigi.Task):
    myid = luigi.Parameter()

    def run(self):
        mypath = paths.Paths07(self.myid)
        lig_ifn = mypath.sdf
        prt_ifn = mypath.pdb

        lig = pybel.readfile("sdf", lig_ifn).next()
        lig.removeh()
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('prt', prt_ifn)

        seq = Structure2Seq(lig, structure)
        with open(self.output().path, 'w') as ofs:
            ofs.write(seq)

    def output(self):
        mypath = paths.Paths07(self.myid)
        ofn = os.path.join(mypath.working,
                           "{}.07.tokens.txt".format(self.myid))
        return luigi.LocalTarget(ofn)


def test():
    luigi.build([BuildTokens("1ajx")], local_scheduler=True)


def main(myid):
    luigi.build([BuildTokens(myid)], local_scheduler=True)


if __name__ == '__main__':
    import sys
    main(sys.argv[1])