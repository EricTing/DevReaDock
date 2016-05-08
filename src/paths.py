import os

PDBBind_core_14 = "/ddnB/work/jaydy/dat/pdbbind_v2014_core_set"
PDBBind_07 = "/ddnB/work/jaydy/dat/v2007"

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
