'''
Usage: python3 interface.py <pdbid>

Finds interfaces in pdbid, then tries to cealign those interface
residues to the reference interface as defined in
reference_interface_definition.py.
'''
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
from reference_interface_definition import *

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

def cealign(reference, query, reference_reslist, query_reslist, pdbid, window=8):
    pymol.cmd.hide('lines','all')
    query_selstr = "{} and ({})".format(pdbid, reslist_to_selstr(query_reslist))
    pymol.cmd.show('lines',query_selstr)
    reference_selstr = "reference and ({})".format(reslist_to_selstr(reference_reslist))
    pymol.cmd.show('lines',reference_selstr)
    alignment_str = "cealign {}, {}".format(reference_selstr,
            query_selstr)
    pymol.util.cbc(selection=pdbid)
    pymol.cmd.color('white','reference and name c*')
    return pymol.cmd.cealign(reference_selstr, query_selstr, window=window)


def align_interfaces(query_interface, reference_interfaces, pdbid,
        reference_pdb, input_dir='test_inputs'):
    query_pymol = pymol.cmd.load(os.path.join(input_dir, pdbid +
        '.clean.pdb'), pdbid)
    #print(pymol.cmd.get_title(query_pymol, 0))
    reference_pymol = pymol.cmd.load(reference_pdb, 'reference')
    i = 0
    best_i = 0
    best_rmsd = 999
    for query_interface in query_interface.pdb_interfaces:
        for reference_interface in reference_interfaces:
            try:
                alignment = cealign(reference_pymol, query_pymol, reference_interface,
                        query_interface, pdbid, window=3)
                print(alignment)
                if alignment['RMSD'] < 3.0:
                    pymol.cmd.save('{}_{}_rmsd_{}.pse'.format(pdbid,i,alignment['RMSD']))
                if alignment['RMSD'] < best_rmsd:
                    best_rmsd = alignment['RMSD']
                    best_i = i
            except:
                print('could not align {} with {}'.format(query_interface, reference_interface))

            i += 1
    print('Best interface: test_{}.pse'.format(best_i))

'''
For testing
'''
if __name__=='__main__':
    init()
    pymol.finish_launching(['pymol','-qc'])
    input_dir = 'test_inputs'

    pdbid = sys.argv[1]
    pose = pose_from_rcsb(pdbid, 'test_inputs')
    interface = PyInterface(pose)
    interface.find_interface()

    #query_resselectors = []

    #reference_selector = list_to_res_selector(reference_interface)
    
    align_interfaces(interface, reference_interfaces, pdbid,
        reference_pdb)
