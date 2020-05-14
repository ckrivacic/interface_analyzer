#! /wynton/home/kortemme/krivacic/software/anaconda36/bin/python3
#$ -l mem_free=4G
#$ -cwd
"""
Usage: relax_preys.py <xml> [options]

Options:
    --tasknum=NUMBER  Just run a specific task
    --parent_dir=STR  Where are the pdb files?
    [default: /wynton/home/kortemme/krivacic/sars/figures/pymol_sessions]
    --outdir=STR  Where to put the output files?  [default: relaxed_outputs]
    """
import glob, os, docopt, subprocess, sys


def run_command(command):
    print('Working directory:', os.getcwd())
    print('Command:', ' '.join(command))
    sys.stdout.flush()

    process = subprocess.Popen(command)

    print('Process ID:', process.pid)
    print()
    sys.stdout.flush()
    
    process.wait()


if __name__=='__main__':
    rosetta_scripts_path = \
            '/wynton/home/kortemme/krivacic/rosetta/source/bin/rosetta_scripts.linuxgccrelease'
            #'/kortemmelab/home/ckrivacic/rosetta/source/bin/rosetta_scripts.linuxgccrelease'
            #'/home/cody/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease'
    rosetta_database_path = \
            '/wynton/home/kortemme/krivacic/rosetta/database'
            #'/kortemmelab/home/ckrivacic/rosetta/database'
            #'/home/cody/rosetta/main/database'
    args = docopt.docopt(__doc__)
    parent_dir = args['--parent_dir']
    pdbs = sorted(glob.glob(parent_dir + '/*/*/*.pdb'))
    outdir = args['--outdir']
    print(pdbs)

    if args['--tasknum']:
        tasknum = int(args['--tasknum'])
    else:
        tasknum = int(os.environ['SGE_TASK_ID']) - 1

    xml = args['<xml>']

    pdb = pdbs[tasknum]
    pdbid = os.path.basename(pdb).split('.')[0].split('_')[0]
    print(pdbid)
    infile = os.path.join('prey_pdbs', pdbid + '.clean.pdb')

    rosetta_cmd = [
            rosetta_scripts_path, 
            '-database', rosetta_database_path, 
            '-in:file:s', infile,
            '-out:prefix', outdir + '/',
            '-out:suffix', '_' + str(tasknum),
            '-out:no_nstruct_label',
            '-parser:protocol', xml,
            '-ex1', '-ex2',
            '-ignore_unrecognized_res'
            ]
    print(rosetta_cmd)
    
    run_command(rosetta_cmd)
