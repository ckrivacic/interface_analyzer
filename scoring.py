from pyrosetta import *
import os, fnmatch, sys
from reference_interface_definition import *
import pandas as pd


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
    return sfxn


def find_pdbs(directory, pattern='*.pdb'):
    """Don't think I use this anywhere. Can probably delete."""
    for root, dirs, files in os.walk(directory):
        for basename in files:
            filename = os.path.join(root, basename)
            if fnmatch.fnmatch(filename, pattern):
                yield filename


def score_pdbs(dataframe, aligner='align'):
    """Score all pdbs in a 'patches' dataframe and update the dataframe."""
    sfxn = motif_scorefxn()
    for idx, row in dataframe.iterrows():
        pdb_path = row['{}_combined_pdb_path'.format(aligner)]
        if pd.isnull(pdb_path):
            continue
        else:
            pose = pose_from_file(pdb_path)
        dataframe.at[idx, 'pose_score'] = sfxn(pose)
        interface_total = 0.0
        residue_dict = {}
        for residue in row['rosetta_full_interface']:
            # Use Rosetta numbering
            residue_energy =\
                    pose.energies().residue_total_energy(residue)
            interface_total += residue_energy
            residue_dict[residue] = residue_energy
        
        patch_total = 0
        for residue in row['rosetta_patch']:
            residue_energy =\
                    pose.energies().residue_total_energy(residue)
            patch_total += residue_energy
            residue_dict[residue] = residue_energy

        reference_total = 0
        for residue in row['rosetta_reference_interface']:
            residue_energy =\
                    pose.energies().residue_total_energy(residue)
            reference_total += residue_energy
            residue_dict[residue] = residue_energy

        # Wrap dict in list so that pandas will accept it
        dataframe.at[idx, '{}_residue_scores'.format(aligner)] = [residue_dict]
        # Save total scores
        dataframe.at[idx, '{}_patch_score'.format(aligner)] = patch_total
        dataframe.at[idx, '{}_interface_score'.format(aligner)] = interface_total

    print(dataframe)


def score_pdb(filename, interface=None, sfxn = None):
    """Score a single pdb. Not used atm."""
    if not sfxn:
        sfxn = motif_scorefxn()
    row_dict = {}
    pose = pose_from_file(filename)
    row_dict['filename'] = filename
    row_dict['pose_total'] = sfxn(pose)
    interface_total = 0.0
    if interface:
        # Interface should have PDB numbering
        for residue in interface:
            rosettanum = pose.pdb_info().pdb2pose('Z', residue)
            residue_energy = pose.energies().residue_total_energy(rosettanum)
            interface_total += residue_energy
            print(pose.energies().residue_total_energy(rosettanum))
            row_dict[residue] = residue_energy
    row_dict['interface_total'] = interface_total

    return row_dict

if __name__=='__main__':
    init('-docking_low_res_score motif_dock_score \
-mh:path:scores_BB_BB \
/home/krivacic/rosetta/database/additional_protocol_data/motif_dock/xh_16_ \
-mh:score:use_ss1 false \
-mh:score:use_ss2 false \
-mh:score:use_aa1 true \
-mh:score:use_aa2 true')
    df = pd.read_pickle('outputs/Q969X5/patches_1.pkl')
    score_pdbs(df)
    df.to_pickle('test_pickle.pkl')
