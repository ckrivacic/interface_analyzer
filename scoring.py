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
    for root, dirs, files in os.walk(directory):
        for basename in files:
            filename = os.path.join(root, basename)
            if fnmatch.fnmatch(filename, pattern):
                yield filename


def score_pdbs(directory, interface=None):
    sfxn = motif_scorefxn()
    df_list = []
    for filename in find_pdbs(directory):
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

        df_list.append(row_dict)

    print(pd.DataFrame(df_list))


def score_pdb(filename, interface=None, sfxn = None):
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
    score_pdbs(sys.argv[1], reference_interface_all)
