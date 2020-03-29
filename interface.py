from pyrosetta import *
from utils import *
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.protocols.scoring import Interface
from pyrosetta.rosetta.protocols.fold_from_loops.movers import AlignByResidueSelectorMover
import numpy as np
import sys, os
# Requires python3.7. Change below to your pymol path.
#sys.path.insert(1,'/home/krivacic/software/pymol/lib/python3.7/site-packages')
#sys.path.insert(1,'/home/krivacic/software/pymol/lib/python3.7')
import pymol

class PyInterface(object):
    def __init__(self, pose):
        self.pose = pose
        self.scorefxn = create_score_function('ref2015')
        self.scorefxn(self.pose)
        self.cutoff_ = 6
        self.interfaces = []
        self.pdb_interfaces = []

    @property
    def cutoff(self):
        return self.cutoff_
    
    @cutoff.setter
    def cutoff(self, cutoff):
        self.cutoff_ = cutoff
    
    def find_interface(self):
        for jump in range(1, self.pose.fold_tree().num_jump() + 1):
            interface = Interface(jump)
            interface.calculate(self.pose)
            self.interfaces.append(interface)
        self.set_pdb_interface()

    def set_pdb_interface(self):
        self.pdb_interfaces = []
        for interface_obj in self.interfaces:
            for side in interface_obj.pair_list():
                self.pdb_interfaces.append(reslist_to_pdb_numbers(side, pose))

def cealign(reference, query, reference_reslist, query_reslist, window=8):
    query_selstr = "query and ({})".format(reslist_to_selstr(query_reslist))
    pymol.cmd.show('lines',query_selstr)
    reference_selstr = "reference and ({})".format(reslist_to_selstr(reference_reslist))
    pymol.cmd.show('lines',reference_selstr)
    alignment_str = "cealign {}, {}".format(reference_selstr,
            query_selstr)
    return pymol.cmd.cealign(reference_selstr, query_selstr, window=window)

'''
For testing
'''
init()
pymol.finish_launching(['pymol','-qc'])

pdbid = '3zpz'
pose = pose_from_rcsb(pdbid, 'test_inputs')
interface = PyInterface(pose)
interface.find_interface()

#query_resselectors = []

#reference_pose = pose_from_file(os.path.join('test_inputs','query.pdb'))
reference_interface = [126,127,125,19,124,80,79,128,78,18,76,129,11,16,14,81,13,116,12,115,15,9,8,7]
#reference_interface = [2,4,5,6,7]
#reference_selector = list_to_res_selector(reference_interface)
query_pymol = pymol.cmd.load(os.path.join('test_inputs', pdbid +
    '.clean.pdb'), 'query')
reference_pymol = pymol.cmd.load(os.path.join('test_inputs','reference.pdb'), 'reference')

i = 0
best_i = 0
best_rmsd = 999
for query_interface in interface.pdb_interfaces:
    try:
        alignment = cealign(reference_pymol, query_pymol, reference_interface,
                query_interface, window=3)
        pymol.cmd.save('test_{}.pse'.format(i))
        print(alignment['RMSD'])
        print(alignment)
        if alignment['RMSD'] < best_rmsd:
            best_rmsd = alignment['RMSD']
            best_i = i
    except:
        print('could not align')

    i += 1
print('Best interface: test_{}.pse'.format(best_i))
