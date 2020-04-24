"""
Usage: get_pdbs.py <uniprot_list>
"""
from blast import *
import pickle as pkl
#import docopt
import sys, os, time, wget
import pandas as pd
from pyrosetta import *

def download_and_clean_pdb(pdbid, prefix=None):
    if prefix:
        if not os.path.exists(prefix):
            os.mkdir(prefix)
        path = os.path.join(prefix, pdbid)
    else:
        path = pdbid

    if not os.path.exists(path + '.clean.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
        pyrosetta.toolbox.cleanATOM(path + '.pdb')
        os.remove(path + '.pdb')

def get_blasts(prey, blast_errors=[], retry_errors=False,
        blast_errors_path='blast_errors.pkl', force=False,
        load_only=False):
    if (prey not in blast_errors) or (retry_errors):
        blasts = []
        print('Collecting info on {}'.format(prey))
        seq_pickle = os.path.join('seqs', '{}.pkl'.format(prey))
        if not os.path.exists(seq_pickle):
            seq = get_sequence(prey)
            with open(seq_pickle, 'wb') as f:
                pkl.dump(seq, f)
        else:
            print('opening existing seq file')
            with open(seq_pickle, 'rb') as f:
                seq = pkl.load(f)
        for s in seq:
            pickle_path = os.path.join('blasts',
                    '{}_blast.pkl'.format(s.id))
            if not os.path.exists(pickle_path) or force:
                if load_only:
                    continue
                try:
                    print('Running BLAST on {}'.format(prey))
                    blast = run_blast(str(s.seq))
                    print('')
                    time.sleep(1)
                    if blast is not None:
                        print('Saving BLAST file')
                        with open(pickle_path, 'wb') as f:
                            pkl.dump(blast, f)
                except Exception as e:
                    print('The following error has occurred with prey {}:'.format(prey))
                    print(e)
                    blast_errors.append(prey)
                    with open(blast_errors_path, 'wb') as f:
                        pkl.dump(blast_errors, f)
                    continue
            else:
                print('Opening existing BLAST file')
                with open(pickle_path, 'rb') as f:
                    blast = pkl.load(f)

            blasts.append(blast)
    else:
        print('An error has previously occurred with prey {}'.format(prey))
        blasts = []

    return blasts

if __name__=='__main__':
    datafile = '2020-03-18_Krogan_SARSCoV2_27baits.txt'
    df = pd.read_csv(datafile, sep='\t')
    preys = df['Preys']
    if not os.path.exists('seqs'):
        os.mkdir('seqs')
    if not os.path.exists('blasts'):
        os.mkdir('blasts')
    blasts = []
    blast_errors_path = 'blast_errors.pkl'
    if os.path.exists(blast_errors_path):
        with open(blast_errors_path, 'rb') as f:
            blast_errors = pkl.load(f)
    else:
        blast_errors = []
    retry_errors = False
    for prey in preys:
        new_blasts = get_blasts(prey, blast_errors=blast_errors,
                load_only=True)
            
        blasts.extend(new_blasts)

    for blast in blasts:
        if blast is not None:
            pdbs = blast.getHits(percent_identity=50)
            for pdb in pdbs:
                print('Downloading {}...'.format(pdb))
                try:
                    download_and_clean_pdb(pdb, prefix='prey_pdbs')
                except:
                    print('error downloading {}'.format(pdb))
