#! /wynton/home/kortemme/krivacic/software/anaconda36/bin/python3
#$ -l mem_free=4G
#$ -cwd
"""
Usage: cluster_scoring.py [options]

Options:
    --tasknum=[NUMBER]  Just run a specific task (for testing)
    --dimerize  Just score virus chain + target chain
"""
from scoring import *
import glob, docopt, os, sys
import pandas as pd


def parse_df_path(df_path):
    basename = os.path.basename(df_path)
    aligner = basename.split('.')[0].split('_')[1]
    return aligner


if __name__=="__main__":
    initialize()
    args = docopt.docopt(__doc__)
    dataframes = sorted(glob.glob("outputs/*/*/*.pkl"))
    print(dataframes)
    if args['--tasknum']:
        tasknum = int(args['--tasknum']) - 1
    else:
        tasknum = int(os.environ['SGE_TASK_ID']) - 1

    df_path = dataframes[tasknum]
    dimerize_pose = args['--dimerize']
    if dimerize_pose:
        df_path = df_path[:-4] + '_dimerized.pkl'
    print('OPENING DATAFRAME {}'.format(df_path))
    reference_pdbid = dataframes[tasknum].split('/')[1]
    reference_pdb = os.path.join('virus_pdbs', reference_pdbid.lower() +
            '.clean.pdb')
    aligner = parse_df_path(df_path)
    print(aligner)

    # Shouldn't need this anymore since I functionalized this
    # initialization
    '''
    init("-docking_low_res_score motif_dock_score \
-mh:path:scores_BB_BB  " + \
os.environ['HOME'] +
"/rosetta/database/additional_protocol_data/motif_dock/xh_16_ \
-mh:score:use_ss1 false \
-mh:score:use_ss2 false \
-mh:score:use_aa1 true \
-mh:score:use_aa2 true")
    '''

    df = pd.read_pickle(df_path)
    # Determine reference surface residues
    reference_surface = get_reference_definition(reference_pdb,
            return_surface=True)
    print(reference_surface)
    csvpath = df_path[:-4] + '_scored.csv'
    if not os.path.exists(csvpath):
        score_pdbs(df, aligner=aligner, reference_surface=reference_surface)
    else:
        print('{} already found - skipping scoring'.format(csvpath))
        sys.exit()

    df.to_pickle(df_path)
    df.to_csv(df_path[:-4] + '_scored.csv')
