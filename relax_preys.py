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


chain_dict = {
        '2jkr_cealign.pdb': 'U',
        '2jkr.pdb': 'A',
        '2xa7_cealign.pdb': 'M',
        '3kdp_align_1.pdb': 'D',
        '2a1t_cealign_1.pdb': 'D',
        '1t9g_cealign_1.pdb': 'D',
        '1oqc_cealign_11166_alignscore_n8.pdb': 'D',
        '1oqc_cealign_8933_alignscore_n8.pdb': 'C',
        '3r7c_cealign_1.pdb': 'C',
        '3r7c_cealign_2743_alignscore_n10.pdb': 'A',
        '3r7c_cealign_5059_alignscore_n9.pdb': 'C',
        '3r7c_cealign_5060_alignscore_n9.pdb': 'C',
        '3r7c_cealign_6343_alignscore_n10.pdb': 'D',
        '3r7c_cealign_6408_alignscore_n5.pdb': 'D',
        '6c4d_align_8494_alignscore_p5.pdb': 'D',
        '6c4d_cealign_8351_alignscore_n12.pdb': 'C',
        'b4y_align_12196_alignscore_n1.pdb': 'A',
        '4g1c_align_10121_alignscore_3.pdb': 'B',
        '4gf1c_cealign_8465_alignscore_n6.pdb': 'A',
        '4utn_cealign_14627_alignscore_n18.pdb': 'B',
        '6eo0_cealign_25566_alignscore_n10.pdb': 'D',
        '4lhx_align_1.pdb': 'D',
        '4lhx_align_2.pdb': 'D',
        '4lhx_cealign_1.pdb': 'D',
        '2ocy_align_1.pdb': 'B',
        '4lhx_align_2_alignscore_n5.pdb': 'D',
        '2eqb_cealign_18957_alignscore_n14.pdb': 'C',
        '2ocy_cealign_47658_alignscore_n4.pdb': 'A',
        '2ocy_cealign_48198_alignscore_n6.pdb': 'A',
        '2ocy_cealign_52112_alignscore_n16.pdb': 'B',
        '3l0i_cealign_33202_alignscore_n11.pdb': 'B',
        '4fmb_cealign_21585_alignscore_n10.pdb': 'A',
        '4lhx_cealign_5670_alignscore_n14.pdb': 'D',
        '4lhx_cealign_5671_alignscore_n14.pdb': 'D',
        '5o74_cealign_53863_alignscore_n11.pdb': 'A',
        '6ekk_cealign_68348_alignscore_n10.pdb': 'C'
        }


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

    pdb = pdbs[tasknum//10]
    outdir = outdir + '/' + pdb + '_origin'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    chain = chain_dict[os.path.basename(pdb)]
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
            '-parser:script_vars', 'chain={}'.format(chain),
            '-ex1', '-ex2',
            '-ignore_unrecognized_res'
            ]
    print(rosetta_cmd)
    
    run_command(rosetta_cmd)
