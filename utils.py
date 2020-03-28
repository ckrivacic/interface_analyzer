from pyrosetta import *
from pyrosetta.rosetta.core.select import residue_selector
import os, wget


def pose_from_rcsb(pdbid, prefix=None):
    if prefix:
        path = os.path.join(prefix,pdbid)
    else:
        path = pdbid
    if not os.path.exists(path + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
    pyrosetta.toolbox.cleanATOM(path + '.pdb')
    pose = rosetta.core.import_pose.get_pdb_and_cleanup(path + '.clean.pdb')

    return pose


def list_to_str(l):
    return ','.join(list(str(i) for i in l))


def list_to_res_selector(l):
    return residue_selector.ResidueIndexSelector(list_to_str(l))


def res_selector_to_size_list(resselector):
    size_list = []
    for i, boolean in enumerate(resselector):
        if boolean == True:
            size_list.append(int(i + 1))

    return intlist_to_vector1_size(size_list)


def reslist_to_pdb_numbers(reslist, pose):
    poselist = []
    for res in reslist:
        poselist.append(pose.pdb_info().pose2pdb(res))

    return poselist


def reslist_to_selstr(reslist):
    selection = []
    for res in reslist:
        try:
            splt = res.split(' ')
            resnum = splt[0]
            chain = splt[1]
        except:
            # Default to chain A if given a list of integers instead of
            # resnum and chain. Can add an option to specify chain
            # later.
            resnum = res
            chain = 'A'
        selstr = ' (resi {} and chain {}) '.format(resnum, chain)
        selection.append(selstr)

    return 'or'.join(selection)
