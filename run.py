"""
Usage: run <UniProt> [options]

Options:
    --aligner=STR  Run cealign instead of default align in PyMOL.
    [default: 'align']
"""
from blast import *
from interface import *
import sys
from reference_interface_definition import *
import pickle as pkl
import docopt

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    aligner = args['--aligner']
    uniprot_id = args['<UniProt>']
    init()
    pymol.finish_launching(['pymol','-qc'])
    #pdbid = sys.argv[1]
    #run_blast(pdbid)
    seq = get_sequence(uniprot_id)
    blasts = []
    for s in seq: # Generally should only be 1 seq
        if not os.path.exists('blasts'):
            os.mkdir('blasts')
        pickle_path = os.path.join('blasts','{}_blast.pkl'.format(s.id))
        if not os.path.exists(pickle_path):
            blast = run_blast(str(s.seq))
            if blast:
                with open(pickle_path,'wb') as f:
                    pkl.dump(blast, f)
        else:
            with open(pickle_path,'rb') as f:
                blast = pkl.load(f)

        blasts.append(blast)

    for blast in blasts: # Generally should only be 1 blast record
        pdbs = blast.getHits(percent_identity=70)
        f = open('{}_pdbs.txt'.format(sys.argv[1]), 'w')
        for pdbid in pdbs:
            f.write(pdbid + '\n')
            print('Running alignments on {}'.format(pdbid))
            pymol.cmd.reinitialize()
            pose = pose_from_rcsb(pdbid, 'test_inputs')
            interfaces = PyInterface(pose)
            interfaces.find_interface()

            aligner = PyMOLAligner(aligner, interfaces, reference_interfaces, pdbid,
                    reference_pdb, output_dir=uniprot_id)
            aligner.align_interfaces()
        f.close()
