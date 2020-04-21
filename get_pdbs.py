"""
Usage: get_pdbs.py <uniprot_list>
"""
from blast import *
import pickle as pkl
#import docopt
import sys, os
import pandas as pd
from pyrosetta import *

def download_and_clean_pdb(pdbid, prefix=None):
    if prefix:
        path = os.path.join(prefix, pdbid)
    else:
        path = pdbid

    if not os.path.exists(path + '.clean.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
    pyrosetta.toolbox.cleanATOM(path + '.pdb')


if __name__=='__main__':
    datafile = '2020-03-18_Krogan_SARSCoV2_27baits.txt'
    df = pd.read_csv(datafile, sep='\t')
    preys = df['Preys']
    if not os.path.exists('blasts'):
        os.mkdir('blasts')
    blasts = []
    for prey in preys:
        seq = get_sequence(prey)
        for s in seq:
            pickle_path = os.path.join('blasts',
                    '{}_blast.pkl'.format(s.id))
            if not os.path.exists(pickle_path):
                blast = run_blast(str(s.seq))
                if blast:
                    with open(pickle_path, 'wb') as f:
                        pkl.dump(blast, f)
            else:
                with open(pickle_path, 'rb') as f:
                    blast = pkl.load(f)
            
            blasts.append(blast)

    for blast in blasts:
        pdbs = blast.getHits(percent_identity=50)
        for pdb in pdbs:
            download_and_clean_pdb(pdb, prefix='prey_pdbs')
