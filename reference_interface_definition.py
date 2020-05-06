'''
Define reference interface. 
Required: reference_interfaces  (list of list), 
reference_pdb (path to reference pdb file)
'''
import os, sys
from patches import Patches
from utils import reslist_to_pdb_numbers
from pyrosetta import *


def get_reference_definition(reference_pdb, patches=True,
        return_surface=False):

    reference_interfaces = []
    if not patches:
        reference_interface_all = [126,127,125,19,124,80,79,128,78,18,76,129,11,16,
                14, 81, 13, 116, 12, 115, 15, 9, 8, 7]

        reference_interfaces = [[],[]]
        for resnum in reference_interface_all:
            if resnum < 90:
                reference_interfaces[0].append(resnum)
            else:
                reference_interfaces[1].append(resnum)

        reference_interfaces.append(reference_interface_all)

    else:
        #init()
        reference_pose = pose_from_file(reference_pdb)
        reference_patches = Patches(reference_pose)
        reference_patches.determine_surface_residues()
        reference_patches.map_residues()
        if return_surface:
            return reslist_to_pdb_numbers(reference_patches.reslist,
                    reference_pose)
        for res in reference_patches.reslist:
            reference_interfaces.append(reference_patches.nearest_n_residues(res, 8))

    return reference_interfaces


if __name__=="__main__":
    from utils import *
    for i in get_reference_definition(sys.argv[1]):
        print(reslist_to_selstr(i))
