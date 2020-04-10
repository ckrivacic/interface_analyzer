"""
Usage: run <UniProt> [options]

Options:
    --aligner=STR  Run cealign instead of default align in PyMOL.
    [default: align]
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
            print(blast)
            if blast:
                with open(pickle_path,'wb') as f:
                    pkl.dump(blast, f)
        else:
            with open(pickle_path,'rb') as f:
                blast = pkl.load(f)

        blasts.append(blast)

    reference_pose = pose_from_file(reference_pdb)
    df = empty_interface_dataframe()
    for blast in blasts: # Generally should only be 1 blast record
        pdbs = blast.getHits(percent_identity=69)
        f = open('{}_pdbs.txt'.format(sys.argv[1]), 'a')
        for pdbid in pdbs:
            pickle_path = os.path.join(
                    'outputs', uniprot_id, '{}_{}'.format(pdbid,
                        aligner), 'patches.pkl'
                    )
            if not os.path.exists(pickle_path):
                f.write(pdbid + '\n')
                print('Running alignments on {}'.format(pdbid))
                pymol.cmd.reinitialize()
                pose = pose_from_rcsb(pdbid, 'test_inputs')
                interfaces = PyInterface(pose, reference_pose)
                interfaces.find_interface()
                interfaces.set_reference_interfaces(reference_interfaces)
                interfaces.set_dataframe(df)
                interfaces.find_patches()

                interface_aligner = PyMOLAligner(aligner, interfaces, pdbid,
                        reference_pdb,
                        output_dir=os.path.join('outputs',uniprot_id))
                interface_aligner.align_patches()
                df = interface_aligner.interface.dataframe
                del interface_aligner
            
            else:
                cached_df = pd.read_pickle(pickle_path)
                df.append(cached_df, ignore_index=True)

        f.close()
    df.to_pickle(os.path.join('outputs', uniprot_id, 'dataframe.pkl'))
    df.to_csv(os.path.join('outputs', uniprot_id, 'dataframe.csv'))
