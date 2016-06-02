import os

PDBBind_core_14 = "/ddnB/work/jaydy/dat/pdbbind_v2014_core_set"
PDBBind_07 = "/ddnB/work/jaydy/dat/v2007"
PDBBind_refined_15 = "/work/jaydy/dat/pdbbinddb-refined-2015"

WORKING = "/work/jaydy/working/pdbbind"


class Paths07:
    def __init__(self, myid):
        self.myid = myid
        self.sdf = os.path.join(PDBBind_07, self.myid,
                                "{}_ligand.sdf".format(self.myid))
        self.pdb = os.path.join(PDBBind_07, self.myid,
                                "{}_protein.pdb".format(self.myid))

        self.working = os.path.join(WORKING, self.myid)

        try:
            os.makedirs(self.working)
        except:
            pass


class Paths15:
    def __init__(self, myid):
        "paths for dat of PDBBind 2015"
        self.myid = myid
        self.sdf = os.path.join(PDBBind_refined_15, self.myid,
                                "{}_ligand.sdf".format(self.myid))
        self.pdb = os.path.join(PDBBind_refined_15, self.myid,
                                "{}_protein.pdb".format(self.myid))

        self.working = os.path.join(WORKING, "PDBBind_refined_15", self.myid)

        try:
            os.makedirs(self.working)
        except:
            pass


class Paths15Rnd:
    def __init__(self, myid):
        "paths for the data of PDBBind 2015 with ligands randomized by AutoDock Vina"
        self.myid = myid
        self.sdf = os.path.join(PDBBind_refined_15, self.myid,
                                "{}_ligand-vina.mol2".format(self.myid))
        self.pdb = os.path.join(PDBBind_refined_15, self.myid,
                                "{}_protein.pdb".format(self.myid))

        self.working = os.path.join(WORKING, "PDBBind_refined_15", self.myid)

        try:
            os.makedirs(self.working)
        except:
            pass

