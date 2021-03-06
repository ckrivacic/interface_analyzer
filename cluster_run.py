#! /wynton/home/kortemme/krivacic/software/anaconda36/bin/python3
#$ -l mem_free=4G
#$ -cwd
"""
Usage: cluster_run.py [options]

Options:
    --aligner=STR  Run cealign instead of default align in PyMOL.
    [default: align]
    --percent_id=NUM  What percentage identity for query proteins  [default: 70]
    --tasknum=[NUMBER]  Just run a specific task
    --test  Only do one pdb
"""
from blast import *
from interface import *
import sys
from reference_interface_definition import *
import pickle as pkl
import docopt
import glob


def finish_io(temp, final, prefix=''):
    from shutil import copyfile
    import subprocess
    print('Finishing IO')
    if not os.path.exists(final):
        os.makedirs(final, exist_ok=True)
    outfile = '{}_outputs.tar.gz'.format(prefix)
    cmd = 'tar -czvf ' + os.path.join(temp, outfile) + ' --directory=' + temp + ' .'
    print('initial command...')
    print(cmd)
    cmd = cmd.split(' ')
    print('RUNNING COMMAND {}'.format(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
    result.stdout.decode('utf-8')
    print(result)
    #os.system(cmd)
    copyfile(os.path.join(temp, outfile), 
            os.path.join(final, outfile))

def load_blast(prey):
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
        otherpath = os.path.join(prefix, pdbid + '.pdb')
        if os.path.exists(otherpath):
            pyrosetta.toolbox.cleanATOM(otherpath)
        else:
            return None
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
    #print('ROWNUM: {}'.format(rownum))
    if args['--tasknum']:
        tasknum = int(args['--tasknum'])
    else:
        tasknum = (int(os.environ['SGE_TASK_ID']) - 1)
    rownum = tasknum // 2
    row = df.iloc[rownum]
    if tasknum%2 == 0:
        aligner = 'cealign'
    else:
        aligner = 'align'
    print('Using aligner {}'.format(aligner))
    percent_id = 30
    #aligner = args['--aligner']
    uniprot_id = row['Preys']
    init("-total_threads 1")
    #pymol.finish_launching(['pymol','-qc'])

    blast = load_blast(uniprot_id)
    print('Blast record:')
    print(blast)
    reference_pdb = parse_bait(row)
    reference_interfaces = get_reference_definition(reference_pdb)
    print('REFERENCE INTERFACES MAIN')
    print(reference_interfaces)
    print(reference_pdb)
    print('Loading reference pose from {}'.format(reference_pdb))
    reference_pose = pose_from_file(reference_pdb)
    df = empty_interface_dataframe()

    if 'TMPDIR' in os.environ:
        os_tmp = os.environ['TMPDIR']
    else:
        os_tmp = os.path.join('temp', os.environ['USER'])
    
    bait = row['Bait'][len('SARS-CoV2 '):]
    tempdir = os.path.join(os_tmp, bait, uniprot_id)
    if not os.path.exists(tempdir):
        print('making temp outdirs')
        os.makedirs(tempdir, exist_ok=True)
    outdir = os.path.join('outputs', bait, uniprot_id)
    if not os.path.exists(outdir):
        print('making outdirs')
        os.makedirs(outdir, exist_ok=True)

    pdbs = blast.getHits(percent_identity=percent_id)
    print(pdbs.keys())
    more_pdbs = glob.glob('prey_pdbs/{}*.pdb'.format(uniprot_id))
    print('EXTENDING PDBLIST WITH MODBASE. FINAL PDB LIST:')
    for more_pdb in more_pdbs:
        pdbs[more_pdb.split('/')[-1].split('.')[0]] = []
    print(pdbs.keys())

    for pdbid in pdbs:
        print('Running alignments for {}'.format(pdbid))
        cached_df = empty_interface_dataframe()
        pickle_path = os.path.join(
                'outputs', bait, uniprot_id, '{}_{}'.format(pdbid,
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
            print('POSE HERE ')
            print(pose)
            if pose is None:
                continue
            interfaces = PyInterface(pose, reference_pose)
            interfaces.find_interface()
            interfaces.set_reference_interfaces(reference_interfaces)
            interfaces.set_dataframe(cached_df)
            interfaces.find_patches()
            print('INTERFACE DATAFRAME')
            print(interfaces.dataframe)

            interface_aligner = PyMOLAligner(aligner, interfaces, pdbid,
                    reference_pdb,
                    output_dir=tempdir)
            print('ALIGNING!!')
            interface_aligner.align_patches()
            df = df.append(interface_aligner.interface.dataframe,
                    ignore_index=True)
            del interface_aligner
        if args['--test']:
            break

    finish_io(tempdir, outdir, prefix='{}_{}_{}'.format(
        row['Bait'][len('SARS-CoV2 '):].lower(), aligner, uniprot_id))

    df.to_pickle(os.path.join('outputs', bait, uniprot_id,
        'dataframe_{}.pkl'.format(aligner)))
    df.to_csv(os.path.join('outputs', bait, uniprot_id,
        'dataframe_{}.csv'.format(aligner)))
