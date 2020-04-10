import pandas as pd
import pymol, sys
import matplotlib.cm
from pyrosetta import *

def make_pymol_session(df, idx, out='temp.pse', aligner='align'):
    """Make a pymol session that colors interface residues by
    motif_score + inter-chain vdw. Patch should be in sticks, interface
    in lines."""
    row = df.iloc[idx]

    pdb_path = row['{}_combined_pdb_path'.format(aligner)]
    pdbinfo = pose_from_file(pdb_path).pdb_info()
    #pymol.finish_launching(['pymol', '-qc'])
    pymol.cmd.load(pdb_path, 'combined')
    pymol.util.cbc()
    pymol.util.cbaw('chain Z')

    residue_scores = row['{}_residue_scores'.format(aligner)][0]
    minval = min(residue_scores.values())
    maxval = max(residue_scores.values())

    # Color all scored residues
    for residue in residue_scores:
        pdb_res = pdbinfo.pose2pdb(int(residue)).split(' ')
        resnum = pdb_res[0]
        chain = pdb_res[1]
        color = get_color(residue_scores[residue], minval, maxval)
        pymol.cmd.show('lines', 'resi {} and chain {}'.format(resnum,
            chain))
        pymol.cmd.color(color, 'chain {} and resi {} and name c*'.format(chain, resnum))

    # Show sticks for target patch
    for residue in row['pymol_target_patch']:
        splt = residue.split(' ')
        resnum = splt[0]
        chain = splt[1]
        pymol.cmd.show('sticks', 'resi {} and chain {}'.format(resnum,
            chain))

    for residue in row['pymol_reference_interface']:
        splt = residue.split(' ')
        resnum = splt[0]
        chain = 'Z'
        pymol.cmd.show('sticks', 'resi {} and chain {}'.format(resnum,
            chain))

    pymol.cmd.save(out)



def get_color(value, minval, maxval):
    cmap = matplotlib.cm.get_cmap('coolwarm')
    fraction = (value - minval) / (maxval - minval)
    rgba = cmap(fraction)
    return '0x' + matplotlib.colors.to_hex(rgba)[1:]


if __name__=='__main__':
    init()
    df = pd.read_pickle('test_pickle.pkl')
    make_pymol_session(df, 94)