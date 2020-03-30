from blast import *
from interface import *
import sys
from reference_interface_definition import *

if __name__=='__main__':
    print(reference_interfaces)
    init()
    pymol.finish_launching(['pymol','-qc'])
    #pdbid = sys.argv[1]
    #run_blast(pdbid)
    seq = get_sequence(sys.argv[1])
    blasts = []
    for s in seq: # Generally should only be 1 seq
        blast = run_blast(str(s.seq))
        blasts.append(blast)

    for blast in blasts: # Generally should only be 1 blast record
        print(blast)
        pdbs = blast.getHits(percent_identity=70)
        f = open('{}_pdbs.txt'.format(sys.argv[1]))
        for pdbid in pdbs:
            f.write(pdbid + '\n')
            print('Running alignments on {}'.format(pdbid))
            pymol.cmd.reinitialize()
            pose = pose_from_rcsb(pdbid, 'test_inputs')
            interfaces = PyInterface(pose)
            interfaces.find_interface()

            align_interfaces(interfaces, reference_interfaces, pdbid,
                    reference_pdb)
        f.close()
