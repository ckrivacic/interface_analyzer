'''
Usage: scoring.py <dataframe> [options]

Options:
    --aligner=STR  Which aligner to look at  [default: align]
    --skip=PDBS  Skip these PDBs
    --csv_out=STR  Save a csv file of the dataframe
'''
from pyrosetta import *
import os, fnmatch, sys
from reference_interface_definition import *
import pandas as pd
import docopt
import tqdm


def initialize():
        init('-docking_low_res_score motif_dock_score \
-mh:path:scores_BB_BB ' + \
os.environ['HOME'] + '/rosetta/database/additional_protocol_data/motif_dock/xh_16_ \
-mh:score:use_ss1 false \
-mh:score:use_ss2 false \
-mh:score:use_aa1 true \
-mh:score:use_aa2 true')


def fix_path(df_path, prefix=''):
    pathlist = df_path.split('/')[-5:]
    pathlist[0] = prefix + pathlist[0]
    pathlist.insert(0, 'outputs')
    return os.path.join(*pathlist)


def motif_scorefxn(vdw_weight=1.0):
    '''
    Returns sforefunction with just motif_dock enabled
    '''
    sfxn = ScoreFunction()
    #sfxn = create_score_function('motif_dock_score')
    
    score_manager = rosetta.core.scoring.ScoreTypeManager()
    motif_term = score_manager.score_type_from_name('motif_dock')
    sfxn.set_weight(motif_term, 1)
    #vdw_term = score_manager.score_type_from_name('interchain_vdw')
    vdw_term = score_manager.score_type_from_name('vdw')
    #vdw_term = score_manager.score_type_from_name('fa_rep')
    sfxn.set_weight(vdw_term, vdw_weight)
    
    return sfxn


def find_pdbs(directory, pattern='*.pdb'):
    """Don't think I use this anywhere. Can probably delete."""
    for root, dirs, files in os.walk(directory):
        for basename in files:
            filename = os.path.join(root, basename)
            if fnmatch.fnmatch(filename, pattern):
                yield filename


def score_pdbs(dataframe, aligner='align', skip=[], out='test_out.pkl',
        reference_surface=None, fix_paths=True):
    """Score all pdbs in a 'patches' dataframe and update the dataframe."""
    sfxn = motif_scorefxn()
    save_check = 0
    for idx, row in tqdm.tqdm(dataframe.iterrows(), total=dataframe.shape[0]):
        pdb_path = row['{}_combined_pdb_path'.format(aligner)]
        if pd.isnull(pdb_path):
            continue
        else:
            breakcheck = False
        print('USING PATH {}'.format(pdb_path))
        if fix_paths:
            pdb_path = fix_path(pdb_path)
        if pdb_path.split('/')[-1].split('_')[0] in skip:
            breakcheck = True
        if breakcheck:
            continue
        pose = pose_from_file(pdb_path)
        switch = SwitchResidueTypeSetMover("centroid")
        switch.apply(pose)
        save_check += 1
        dataframe.at[idx, 'pose_score'] = sfxn(pose)
        interface_total = 0.0
        residue_dict = {}
        for residue in row['pymol_target_interface']:
            resnum = int(residue.split(' ')[0])
            chain = residue.split(' ')[1]
            rosettanum = pose.pdb_info().pdb2pose(chain, resnum)
            #print(rosettanum)
            if rosettanum != 0:
                residue_energy =\
                        pose.energies().residue_total_energy(rosettanum)
                interface_total += residue_energy
                residue_dict[rosettanum] = residue_energy
        
        '''
        Patches are on the REPLACED chain so commenting this out
        as that chain is not in this pdb file.
        Should be looking at reference_total instead, which will be more
        useful if reference interfaces are defined as patches instead of
        whole interfaces.
        patch_total = 0
        for residue in row['rosetta_patch']:
            residue_energy =\
                    pose.energies().residue_total_energy(residue)
            patch_total += residue_energy
            residue_dict[residue] = residue_energy
        '''
        # Score entire surface of reference
        if reference_surface:
            ref_surface_total = 0
            for residue in reference_surface:
                resnum = int(residue.split(' ')[0])
                #chain = residue.split(' ')[1]
                chain = 'Z'
                rosettanum = pose.pdb_info().pdb2pose(chain, resnum)
                if rosettanum != 0:
                    residue_energy =\
                            pose.energies().residue_total_energy(rosettanum)
                    ref_surface_total += residue_energy
                    residue_dict[rosettanum] = residue_energy

        # Get total for residues in the interacting patch
        interacting_residues_total = 0
        for residue in row['pymol_target_patch']:
            resnum = int(residue.split(' ')[0])
            chain = residue.split(' ')[1]
            rosettanum = pose.pdb_info().pdb2pose(chain, resnum)
            if rosettanum != 0:
                residue_energy =\
                        pose.energies().residue_total_energy(rosettanum)
                interacting_residues_total += residue_energy
                residue_dict[rosettanum] = residue_energy


        # Total energy for residues in reference chain
        reference_total = 0
        for residue in row['pymol_reference_interface']:
            resnum = int(residue.split(' ')[0])
            rosettanum = pose.pdb_info().pdb2pose('Z', resnum)
            if rosettanum != 0:
                residue_energy =\
                        pose.energies().residue_total_energy(rosettanum)
                reference_total += residue_energy
                residue_dict[rosettanum] = residue_energy
        

        # Wrap dict in list so that pandas will accept it
        dataframe.at[idx, '{}_residue_scores'.format(aligner)] = [residue_dict]
        # Save total scores
        #dataframe.at[idx, '{}_patch_score'.format(aligner)] = patch_total
        dataframe.at[idx, '{}_interface_score'.format(aligner)] = interface_total
        dataframe.at[idx, '{}_reference_score'.format(aligner)] = reference_total
        dataframe.at[idx, '{}_target_patch_score'.format(aligner)] =\
                interacting_residues_total
        dataframe.at[idx, '{}_reference_surface_score'.format(aligner)]\
                = ref_surface_total
        if save_check == 100:
            dataframe.to_pickle(out)
            save_check = 0

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


def score_row(dataframe, idx, aligner='align'):
    """Score a single pdb. Not used atm."""
    sfxn = motif_scorefxn(vdw_weight=2.0)
    row = dataframe.iloc[idx]
    interface_total = 0.0
    pdb_path = row['{}_combined_pdb_path'.format(aligner)]
    pose = pose_from_file(pdb_path)
    switch = SwitchResidueTypeSetMover("centroid")
    switch.apply(pose)
    dataframe.at[idx, 'pose_score'] = sfxn(pose)
    interface_total = 0.0
    residue_dict = {}
    for residue in row['pymol_target_interface']:
        resnum = int(residue.split(' ')[0])
        chain = residue.split(' ')[1]
        rosettanum = pose.pdb_info().pdb2pose(chain, resnum)
        #print(rosettanum)
        if rosettanum != 0:
            residue_energy =\
                    pose.energies().residue_total_energy(rosettanum)
            interface_total += residue_energy
            residue_dict[rosettanum] = residue_energy
    
    '''
    Patches are on the REPLACED chain so commenting this out
    as that chain is not in this pdb file.
    Should be looking at reference_total instead, which will be more
    useful if reference interfaces are defined as patches instead of
    whole interfaces.
    patch_total = 0
    for residue in row['rosetta_patch']:
        residue_energy =\
                pose.energies().residue_total_energy(residue)
        patch_total += residue_energy
        residue_dict[residue] = residue_energy
    '''
    # Get total for residues in the interacting patch
    interacting_residues_total = 0
    for residue in row['pymol_target_patch']:
        resnum = int(residue.split(' ')[0])
        chain = residue.split(' ')[1]
        rosettanum = pose.pdb_info().pdb2pose(chain, resnum)
        if rosettanum != 0:
            print(pose.energies().residue_total_energies(rosettanum))
            residue_energy =\
                    pose.energies().residue_total_energy(rosettanum)
            interacting_residues_total += residue_energy
            residue_dict[rosettanum] = residue_energy


    # Total energy for residues in reference chain
    reference_total = 0
    for residue in row['pymol_reference_interface']:
        resnum = int(residue.split(' ')[0])
        rosettanum = pose.pdb_info().pdb2pose('Z', resnum)
        if rosettanum != 0:
            residue_energy =\
                    pose.energies().residue_total_energy(rosettanum)
            reference_total += residue_energy
            residue_dict[rosettanum] = residue_energy
    

    # Wrap dict in list so that pandas will accept it
    dataframe.at[idx, '{}_residue_scores'.format(aligner)] = [residue_dict]
    # Save total scores
    #dataframe.at[idx, '{}_patch_score'.format(aligner)] = patch_total
    dataframe.at[idx, '{}_interface_score'.format(aligner)] = interface_total
    dataframe.at[idx, '{}_reference_score'.format(aligner)] = reference_total
    dataframe.at[idx, '{}_target_patch_score'.format(aligner)] =\
            interacting_residues_total

    print(dataframe)

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    print(args)
    aligner = args['--aligner']
    if args['--skip']:
        skip = args['--skip'].split(',')
    else:
        skip = []
    input_dataframe = args['<dataframe>']
    init('-docking_low_res_score motif_dock_score \
-mh:path:scores_BB_BB ' + \
os.environ['HOME'] + '/rosetta/database/additional_protocol_data/motif_dock/xh_16_ \
-mh:score:use_ss1 false \
-mh:score:use_ss2 false \
-mh:score:use_aa1 true \
-mh:score:use_aa2 true')
    '''
    if args['--temp']:
        row = df.loc[156]
        row['align_combined_pdb_path'] =
                row['align_combined_pdb_path'[:-4] + '_0001.pdb']
    '''
    df = pd.read_pickle(input_dataframe)
    # Determine reference surface residues
    reference_surface = reference_patches.reslist
    score_pdbs(df, aligner=aligner, skip=skip, reference_surface=reference_surface)

    #score_pdbs(df)
    #score_row(df, 6962)
    df.to_pickle(input_dataframe)
    if args['--csv_out']:
        df.to_csv(args['--csv_out'])
