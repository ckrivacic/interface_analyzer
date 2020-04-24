"""
Usage: run [options]

Options:
    --aligner=STR  Run cealign instead of default align in PyMOL.
    [default: align]
    --percent_id=NUM  What percentage identity for query proteins  [default: 70]
"""
from blast import *
from interface import *
import sys
from reference_interface_definition import *
import pickle as pkl
import docopt


def load_blast(prey):
    print('loading blast for {}'.format(prey))
    seq_path = os.path.join('seqs', '{}.pkl'.format(prey))
    if not os.path.exists(seq_path):
        return None
    
    with open(seq_path, 'rb') as f:
        seq = pkl.load(f)
    
    s = seq[0]
    blast_path = os.path.join('blasts', '{}_blast.pkl'.format(s.id))
    if not os.path.exists(blast_path):
        return None

    with open(blast_path, 'rb') as f:
        blast = pkl.load(f)
        return blast

def load_clean_pose(pdbid, prefix='prey_pdbs'):
    path = os.path.join(prefix, pdbid + '.clean.pdb')
    if not os.path.exists(path):
        return None
    else:
        return rosetta.core.import_pose.get_pdb_and_cleanup(path)


def parse_bait(row, prefix='SARS-CoV2 ', folder='virus_pdbs'):
    baitname = row['Bait']
    baitname = baitname[len(prefix):].lower()
    base = os.path.join(folder, baitname)
    if not os.path.exists(base + '.clean.pdb'):
        pyrosetta.toolbox.cleanATOM(base + '.pdb')
    return base + '.clean.pdb'


if __name__=='__main__':
    args = docopt.docopt(__doc__)
    datafile = '2020-03-18_Krogan_SARSCoV2_27baits.txt'
    df = pd.read_csv(datafile, sep='\t')
    rownum = int(os.environ['SGE_TASK_ID'] - 1) // 2
    row = df.iloc[rownum]
    aligner = 'cealign' if (rownum%2 == 0) else 'align'
    #row = df.iloc[rownum]
    percent_id = int(args['--percent_id'])
    aligner = args['--aligner']
    uniprot_id = row['Preys']
    init()
    pymol.finish_launching(['pymol','-qc'])

    blast = load_blast(uniprot_id)
    reference_pdb = parse_bait(row)
    reference_pose = pose_from_file(reference_pdb)
    df = empty_interface_dataframe()

    outdir = os.path.join('outputs', row['Bait'], uniprot_id)
    if not os.path.exists(outdir):
        print('making outdirs')
        os.makedirs(outdir, exist_ok=True)

    pdbs = blast.getHits(percent_identity=percent_id)
    print(pdbs)
    for pdbid in pdbs:
        cached_df = empty_interface_dataframe()
        pickle_path = os.path.join(
                'outputs', row['Bait'], uniprot_id, '{}_{}'.format(pdbid,
                    aligner), 'patches.pkl'
                )
        if os.path.exists(pickle_path):
            cached_df = pd.read_pickle(pickle_path)
            done = '{}_rmsd'.format(aligner) in cached_df.columns
            #df = df.append(cached_df, ignore_index=True)
        if not os.path.exists(pickle_path) or not done:
            print('Running alignments on {}'.format(pdbid))
            pymol.cmd.reinitialize()
            pose = load_clean_pose(pdbid)
            if pose is None:
                continue
            interfaces = PyInterface(pose, reference_pose)
            interfaces.find_interface()
            interfaces.set_reference_interfaces(reference_interfaces)
            interfaces.set_dataframe(cached_df)
            interfaces.find_patches()
            print(interfaces.dataframe)

            interface_aligner = PyMOLAligner(aligner, interfaces, pdbid,
                    reference_pdb,
                    output_dir=os.path.join('outputs', row['Bait'], uniprot_id))
            interface_aligner.align_patches()
            df = df.append(interface_aligner.interface.dataframe,
                    ignore_index=True)
            del interface_aligner

    df.to_pickle(os.path.join('outputs', row['Bait'], uniprot_id,
        'dataframe_{}.pkl'.format(aligner)))
    df.to_csv(os.path.join('outputs', row['Bait'], uniprot_id,
        'dataframe_{}.csv'.format(aligner)))
