from prody import *
from Bio import SwissProt
import sys
import requests as r
from Bio import SeqIO
from io import StringIO

def run_blast(sequence):
    #return blastPDB(sequence)
    return blastPDB(sequence, timeout=240)


def get_sequence(cID):

    #cID='P04637'

    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)

    Seq=StringIO(cData)
    return list(SeqIO.parse(Seq,'fasta'))


if __name__=='__main__':
    #pdbid = sys.argv[1]
    #run_blast(pdbid)
    seq = get_sequence(sys.argv[1])
    print(seq)
    blasts = []
    for s in seq:
        print(str(s.seq))
        blast = run_blast(str(s.seq))
        blasts.append(blast)

    for blast in blasts:
        print(blast)
        pdbs = blast.getHits(percent_identity=70)

