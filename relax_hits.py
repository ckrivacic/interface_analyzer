from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta import *


def setup_movemap_from_resselectors(designable_selector, repackable_selector):
    """
    Temporary function. Ultimately we want a more sophisticated movemap,
    probably a neighbor movemap or something using clash-based shell selector.
    """
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in range(1, len(repackable_selector) + 1):
        if designable_selector[i] or repackable_selector[i]:
            mm.set_bb(i, True)
            mm.set_chi(i, True)

    # for i in residues_sc_movable:
    #    mm.set_chi(i, True)

    return mm


def setup_movemap(residues_bb_movable, residues_sc_movable):
    mm = rosetta.core.kinematics.MoveMap()

    if (residues_bb_movable is None) and (residues_sc_movable is None):
        # Default to everything movable
        mm.set_chi(True)
        mm.set_bb(True)
    else:
        mm.set_chi(False)
        mm.set_bb(False)

        for i in residues_bb_movable:
            mm.set_bb(i, True)
            mm.set_chi(i, True)

        for i in residues_sc_movable:
            mm.set_chi(i, True)

    return mm


def fast_relax(pose, residues_bb_movable=None, residues_sc_movable=None,
        selectors=False):
    '''Fast relax the pose'''
    if selectors==False:
        mm = setup_movemap(residues_bb_movable, residues_sc_movable)
    else:
        mm = setup_movemap_from_resselectors(residues_bb_movable,
                residues_sc_movable)
    #sfxn = setup_restrained_sfxn(['coordinate_constraint'],[2.0])
    sfxn = create_score_function('ref2015')

    fast_relax_rounds = 5
    fast_relax = rosetta.protocols.relax.FastRelax(sfxn, fast_relax_rounds)
    fast_relax.set_movemap(mm) 
    
    #fast_relax.apply(pose)
    return fast_relax


def relax_pose(pose):
    relaxer = fast_relax(pose, None, None)
    relaxer.apply(pose)


def score_dg(pose):
    #chain = pyrosetta.rosetta.std.set_int_t()
    chain = pose.pdb_info().pdb2pose('Z', 1)
    chain_poses = pose.split_by_chain()
    for chain in chain_poses:
        ch = chain.pdb_info().pose2pdb(1).split(' ')[1]
        print(ch)
    it = chainmap.begin()
    for ch in chainmap:
        print(ch.first)
    #interface = InterfaceAnalyzerMover(chainset.insert(chain))


if __name__=='__main__':
    init()
    testfile = \
            '/home/cody/sars/pymol_sessions/Nsp7/P61026/2eqb_cealign_18957_alignscore_n14.pdb'
    pose = pose_from_file(testfile)
    #fast_relax(pose)
    score_dg(pose)
