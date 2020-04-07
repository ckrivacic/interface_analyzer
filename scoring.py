from pyrosetta import *
import os, fnmatch, sys


def motif_scorefxn():
    '''
    Returns sforefunction with just motif_dock enabled
    '''
    sfxn = ScoreFunction()
    score_manager = rosetta.core.scoring.ScoreTypeManager()
    motif_term = score_manager.score_type_from_name('motif_dock')
    sfxn.set_weight(motif_term, 1)
    vdw_term = score_manager.score_type_from_name('interchain_vdw')
    sfxn.set_weight(vdw_term, 1)
    print(sfxn.get_nonzero_weighted_scoretypes())
    return sfxn


def find_pdbs(directory, pattern='*/combined/*.pdb'):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            filename = os.path.join(root, basename)
            if fnmatch.fnmatch(filename, pattern):
                yield filename


def score_pdbs(directory):
    sfxn = motif_scorefxn()
    out_dict = {}
    for filename in find_pdbs(directory):
        pose = pose_from_file(filename)
        out_dict[filename] = sfxn(pose)
    print(out_dict)

if __name__=='__main__':
    init('-docking_low_res_score motif_dock_score \
-mh:path:scores_BB_BB \
/home/krivacic/rosetta/database/additional_protocol_data/motif_dock/xh_16_ \
-mh:score:use_ss1 false \
-mh:score:use_ss2 false \
-mh:score:use_aa1 true \
-mh:score:use_aa2 true')
    score_pdbs(sys.argv[1])
