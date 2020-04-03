from pyrosetta.rosetta.core.select import residue_selector
import sys
import pandas as pd
from pyrosetta import *
from utils import *
from numeric import *

class Patches(object):
    '''
    Class which holds information about protein interface patches
    '''
    def __init__(self, pose):
        if pose.num_chains() > 1:
            print("WARNING: Pose has more than one chain. Interface "\
                    "residues between chains may not be counted as exposed,"\
                    " and therefore won't be used to determine surface "\
                    "patches. For better results, pass a single chain "\
                    "using pose.split_by_chain(chain_num)")
        self.pose = pose
        self.reslist = None

    def set_reslist(self, reslist):
        '''
        Give the object your own reslist. Useful if you want to only
        align to interfaces.
        Make sure to run determine_surface_residues if you want to
        narrow down this reslist to solvent-exposed residues.
        '''
        if type(reslist[0])==type(True):
            self.reslist = res_selector_to_size_list(reslist)
        else:
            self.reslist = intlist_to_vector1_size(reslist) 
    def determine_surface_residues(self):
        surface_selector = residue_selector.LayerSelector()
        surface_selector.set_layers(False, False, True)
        reslist =\
                res_selector_to_size_list(surface_selector.apply(self.pose))
        if self.reslist:
            # If a reslist is already defined, only  take surface
            # residues that are in that reslist
            reslist = [res for res in reslist if res in self.reslist]
        self.reslist = reslist

    def map_residues(self):
        resmap = {}
        for res1 in self.reslist:
            if res1 not in resmap:
                resmap[res1] = {}
            for res2 in self.reslist:
                if res2 in resmap: # Don't calculate twice
                    if res1 in resmap[res2]:
                        resmap[res1][res2] = resmap[res2][res1]
                else:
                    xyz1 = self.pose.residue(res1).xyz('CA')
                    xyz2 = self.pose.residue(res2).xyz('CA')
                    resmap[res1][res2] = euclidean_distance(xyz1, xyz2)
        resmap = pd.DataFrame(resmap).fillna(0).unstack().reset_index()
        resmap.columns = ['res1', 'res2', 'dist']
        self.resmap = resmap

    def nearest_n_residues(self, resnum, n, cutoff=15.0):
        neighbors = self.resmap[(self.resmap['res1']==resnum) &
                (self.resmap['dist'] <
                    cutoff)].sort_values(by='dist')['res2']
        return neighbors[0:n].tolist()


if __name__=='__main__':
    init()
    pdbid = sys.argv[1]
    pose = pose_from_rcsb(pdbid, 'test_inputs')

    patches = Patches(pose)
    #patches.set_reslist([7,11,14, 8])
    patches.determine_surface_residues()
    print('reslist:')
    print(patches.reslist)
    
    patches.map_residues()
    print(patches.resmap)
    #print(patches.resmap)
    print(patches.nearest_n_residues(11, 10))