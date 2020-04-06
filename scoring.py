from pyrosetta import *
import os, fnmatch, sys


def motif_scorefxn():
    '''
    Returns sforefunction with just motif_dock enabled
    '''
    sfxn = ScoreFunction()
    score_manager = rosetta.core.scoring.ScoreTypeManager()
    score_term = score_manager.score_type_from_name('motif_dock')
    sfxn.set_weight(score_term, 1)
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
    init()
    score_pdbs(sys.argv[1])
