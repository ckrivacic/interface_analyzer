"""
Usage: run <UniProt> [options]

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

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    percent_id = int(args['--percent_id'])
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
        pdbs = blast.getHits(percent_identity=percent_id)
        f = open('{}_pdbs.txt'.format(sys.argv[1]), 'a')
        for pdbid in pdbs:
            cached_df = empty_interface_dataframe()
            pickle_path = os.path.join(
                    'outputs', uniprot_id, '{}_{}'.format(pdbid,
                        aligner), 'patches.pkl'
                    )
            if os.path.exists(pickle_path):
                cached_df = pd.read_pickle(pickle_path)
                done = '{}_rmsd' in cached_df.columns
                #df = df.append(cached_df, ignore_index=True)
            if not os.path.exists(pickle_path) or not done:
                f.write(pdbid + '\n')
                print('Running alignments on {}'.format(pdbid))
                pymol.cmd.reinitialize()
                pose = pose_from_rcsb(pdbid, 'test_inputs')
                interfaces = PyInterface(pose, reference_pose)
                interfaces.find_interface()
                interfaces.set_reference_interfaces(reference_interfaces)
                interfaces.set_dataframe(cached_df)
                interfaces.find_patches()
                print(interfaces.dataframe)

                interface_aligner = PyMOLAligner(aligner, interfaces, pdbid,
                        reference_pdb,
                        output_dir=os.path.join('outputs',uniprot_id))
                interface_aligner.align_patches()
                df = df.append(interface_aligner.interface.dataframe,
                        ignore_index=True)
                del interface_aligner

        f.close()
    df.to_pickle(os.path.join('outputs', uniprot_id,
        'dataframe_{}.pkl'.format(aligner)))
    df.to_csv(os.path.join('outputs', uniprot_id,
        'dataframe_{}.csv'.format(aligner)))
