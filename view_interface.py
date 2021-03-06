'''
Usage: view_interface.py <dataframe> [options]

Options:
    --aligner=STR  Which aligner to look at  [default: align]
    --path=STR  Use a different pdb path
    --row=INT, -r  Which row to plot
    --start=INT  When going through dataframe, start at this point
'''

import pandas as pd
from scoring import dimerize
import pymol, sys
import matplotlib.cm
from pyrosetta import *

def make_pymol_session(dfpath, idx, out='temp.pse', aligner='align',
        pdb_path=None):
    """Make a pymol session that colors interface residues by
    motif_score + inter-chain vdw. Patch should be in sticks, interface
    in lines."""
    pymol.cmd.reinitialize()
    df = pd.read_pickle(dfpath)
    row = df.iloc[idx]

    if not pdb_path:
        pdb_path = row['{}_combined_pdb_path'.format(aligner)]

    def parse_pdb_path(path):
        split = path.split('/')[-5:]
        split.insert(0, 'outputs')
        return os.path.join(*split)

    pdb_path = parse_pdb_path(pdb_path)
    #pdb_path = 'outputs/P55789/3r7c_align/combined/3r7c_48_length_13_rmsd_1.95_0001.pdb'
    pose = pose_from_file(pdb_path)
    if 'dimerized' in dfpath:
        print('DIMERIZING POSE')
        pose = dimerize(pose, row['pymol_target_patch'])
    pdbinfo = pose.pdb_info()
    #pymol.finish_launching(['pymol', '-qc'])
    basename = os.path.basename(pdb_path)
    pymol.cmd.load(pdb_path, basename)
    #pymol.util.cbc()
    pymol.util.cbao('chain Z')

    residue_scores = row['{}_residue_scores'.format(aligner)][0]
    #minval = min(residue_scores.values())
    minval = -5.0
    #maxval = max(residue_scores.values())
    maxval = 5.0

    # Color all scored residues
    make_colorbar = False
    for residue in residue_scores:
        pdb_res = pdbinfo.pose2pdb(int(residue)).split(' ')
        resnum = pdb_res[0]
        chain = pdb_res[1]
        color = get_color(residue_scores[residue], minval, maxval,
                make_colorbar=make_colorbar)
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

    pymol.cmd.remove('hydro')

    pymol.cmd.save(out)
    os.system('pymol temp.pse')



def get_color(value, minval, maxval, make_colorbar=False):
    cmap = matplotlib.cm.get_cmap('coolwarm')
    if make_colorbar:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 1))
        fig.subplots_adjust(bottom=0.5)
        norm = matplotlib.colors.Normalize(vmin=minval, vmax=maxval)
        cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                orientation='horizontal')
        cb1.set_label('Per-residue REU')
        fig.show()
        plt.show()
    fraction = (value - minval) / (maxval - minval)
    rgba = cmap(fraction)
    return '0x' + matplotlib.colors.to_hex(rgba)[1:]


if __name__=='__main__':
    import docopt

    args = docopt.docopt(__doc__)
    init()
    pdb_path = args['--path']
    dfpath = args['<dataframe>']
    if args['--row']:
        make_pymol_session(dfpath, int(args['--row']),
                aligner=args['--aligner'], pdb_path=pdb_path)
    else:
        if args['--start']:
            start = int(args['--start'])
        else:
            start = 0
        df = pd.read_pickle(dfpath)
        j = 0
        for idx, row in df.sort_values('pose_score').iterrows():
            if j < start:
                j += 1
                continue
            print('INDEX IS', idx)
            print('PDB PATH:',
                    row['{}_combined_pdb_path'.format(args['--aligner'])])
            make_pymol_session(dfpath, idx, aligner=args['--aligner'],
                    pdb_path=pdb_path)
            j += 1
    #os.system('pymol temp.pse')
