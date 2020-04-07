'''
Usage: python3 interface.py <pdbid>

Finds interfaces in pdbid, then tries to cealign those interface
residues to the reference interface as defined in
reference_interface_definition.py.
'''
from pyrosetta import *
from utils import *
from numeric import *
from patches import Patches
import pandas as pd
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.protocols.scoring import Interface
from pyrosetta.rosetta.protocols.fold_from_loops.movers import AlignByResidueSelectorMover
import numpy as np
import sys, os
# Requires python3.7. Change below to your pymol path.
#sys.path.insert(1,'/home/krivacic/software/pymol/lib/python3.7/site-packages')
#sys.path.insert(1,'/home/krivacic/software/pymol/lib/python3.7')
import pymol
from reference_interface_definition import *


class PyInterface():
    '''
    Class to handle interface definitions.
    '''
    def __init__(self, pose):
        self.pose = pose
        scorefxn = create_score_function('ref2015')
        scorefxn(self.pose)
        self.cutoff_ = 8.0
        self.patch_residues = 5
        self.patch_cutoff = 20
        self.interfaces = []
        self.pdb_interfaces = []
        self.dataframe = None

    @property
    def cutoff(self):
        return self.cutoff_
    
    @cutoff.setter
    def cutoff(self, cutoff):
        self.cutoff_ = cutoff
    
    def find_interface(self):
        for jump in range(1, self.pose.fold_tree().num_jump() + 1):
            interface = Interface(jump)
            interface.distance(self.cutoff_)
            interface.calculate(self.pose)
            self.interfaces.append(interface)
        self.set_pdb_interface()

    def set_pdb_interface(self):
        self.pdb_interfaces = []
        for interface_obj in self.interfaces:
            for side in interface_obj.pair_list():
                self.pdb_interfaces.append(reslist_to_pdb_numbers(side,
                    self.pose))

    def find_patches(self):
        '''
        Build dataframe of patches for all interfaces
        '''
        df_list = []
        for interface in self.interfaces:
            # Get PDB info from interface
            pair_list = vector1_to_python_list(interface.pair_list())
            sideA = vector1_to_python_list(pair_list[0])
            sideA_pymol = reslist_to_pdb_numbers(sideA, self.pose)
            sideB = vector1_to_python_list(pair_list[1])
            sideB_pymol = reslist_to_pdb_numbers(sideB, self.pose)

            chainA = pose.pdb_info().pose2pdb(sideA[0]).split(' ')[1]
            rosetta_chainA = pose.chain(sideA[0])
            chainB = pose.pdb_info().pose2pdb(sideB[0]).split(' ')[1]
            rosetta_chainB = pose.chain(sideB[0])

            # Get patches for side A
            patches = Patches(self.pose)
            # Limit residues to those found by Interface obj
            patches.set_reslist(sideA)
            # Further limit residues to surface residues
            patches.determine_surface_residues()
            patches.map_residues()

            # Build row for side A
            for resi in sideA:
                patch = patches.nearest_n_residues(resi,
                        self.patch_residues, cutoff=self.patch_cutoff)
                # No need to add a patch twice
                if not any(d['rosetta_patch'] == patch for d in
                        df_list):
                    row = {
                            'pymol_chain': chainA,
                            'rosetta_chain': rosetta_chainA,
                            'pymol_tarchain': chainB,
                            'rosetta_tarchain': rosetta_chainB,
                            'rosetta_full_interface': sideA,
                            'pymol_full_interface': sideA_pymol,
                            'rosetta_patch': patch,
                            'pymol_patch': reslist_to_pdb_numbers(patch,
                                self.pose)
                            }
                    df_list.append(row)

            # Get patches for side B
            del patches
            patches = Patches(self.pose)
            # Limit residues to those found by Interface obj
            patches.set_reslist(sideB)
            # Further limit residues to surface residues
            patches.determine_surface_residues()
            patches.map_residues()

            # Build row for side B
            for resi in sideB:
                patch = patches.nearest_n_residues(resi,
                        self.patch_residues, cutoff=self.patch_cutoff)
                # No need to add a patch twice
                if not any(d['rosetta_patch'] == patch for d in
                        df_list):
                    row = {
                            'pymol_chain': chainB,
                            'rosetta_chain': rosetta_chainB,
                            'pymol_tarchain': chainA,
                            'rosetta_tarchain': rosetta_chainA,
                            'rosetta_full_interface': sideB,
                            'pymol_full_interface': sideB_pymol,
                            'rosetta_patch': patch,
                            'pymol_patch': reslist_to_pdb_numbers(patch,
                                self.pose)
                            }
                    df_list.append(row)
        self.dataframe = pd.DataFrame(df_list)


class PyMOLAligner(object):
    def __init__(self, aligner, query_interface, reference_interfaces,
            pdbid, reference_pdb, output_dir='.',
            input_dir='test_inputs', window=3, cycles=5):
        # String, either 'cealign' or 'align'
        self.aligner=aligner
        # PyInterface object
        self.query_interface = query_interface
        # Reference interfaces list
        self.reference_interfaces = reference_interfaces
        # PDBID for naming
        self.pdbid = pdbid
        # Path to reference PDB
        self.reference_pdb = reference_pdb
        # IO
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.input_dir = input_dir
        # Parameter for align
        self.cycles=cycles
        # Parameter for cealign
        self.window=3
        self.query_pymol = pymol.cmd.load(os.path.join(self.input_dir, self.pdbid +
            '.clean.pdb'), self.pdbid)
        #print(pymol.cmd.get_title(query_pymol, 0))
        self.reference_pymol = pymol.cmd.load(self.reference_pdb, 'reference')
        pymol.cmd.alter('reference', "chain='Z'")

    def align(self, reference_reslist, query_reslist):
        '''
        Align a reference patch with a query patch using the defined
        aligner.
        '''
        pymol.cmd.hide('lines','all')
        query_selstr = "{} and ({})".format(self.pdbid, reslist_to_selstr(query_reslist))
        pymol.cmd.show('lines', query_selstr)
        reference_selstr = "reference and ({})".format(reslist_to_selstr(reference_reslist, chain='Z'))
        pymol.cmd.show('lines', reference_selstr)
        alignment_str = "cealign {}, {}".format(reference_selstr,
                query_selstr)
        pymol.util.cbc(selection=self.pdbid)
        pymol.cmd.color('white','reference and name c*')
        # return pymol.cmd.cealign(reference_selstr, query_selstr, window=window)
        if self.aligner == 'cealign':
            return pymol.cmd.cealign(reference_selstr, query_selstr,
                                   window=self.window)
        elif self.aligner == 'align':
            return pymol.cmd.align(reference_selstr, query_selstr,
                                   cycles=self.cycles)
        else:
            # Default to align
            return pymol.cmd.align(reference_selstr, query_selstr,
                                   cycles=self.cycles)

    def align_interfaces(self, align_patches=True):
        '''
        Align all interfaces defined in reference and query proteins.
        '''
        i = 0
        best_i = 0
        best_rmsd = 999
        formatted_outdir = os.path.join(self.output_dir,
                '{}_{}'.format(self.pdbid, self.aligner))
        if align_patches:
            query_iterator = self.query_interface.dataframe['pymol_patch']
        else:
            query_iterator = self.query_interface.dataframe['pymol_full_interface'].unique()
        for query_interface in query_iterator:
            if len(query_interface) == 0:
                continue
            query_chain = query_interface[0].split(' ')[1]
            for reference_interface in self.reference_interfaces:
                try:
                    alignment = self.align(reference_interface,
                            query_interface)
                    print(i)
                    print(alignment)
                    if self.aligner=='cealign':
                        if alignment['RMSD'] < 3.0:
                            if not os.path.exists(formatted_outdir):
                                os.mkdir(formatted_outdir)
                            if not os.path.exists(os.path.join(formatted_outdir,
                                'combined')):
                                os.mkdir(os.path.join(formatted_outdir,
                                    'combined'))
                            pymol.cmd.center('reference')
                            name = '{}_{}_length_{}_rmsd_{:.2f}'.format(self.pdbid,
                                    i, alignment['alignment_length'],
                                    alignment['RMSD'])
                            pymol.cmd.save(os.path.join(formatted_outdir,
                                name + '.pse'))
                            pymol.cmd.create('combined', 'reference or ('
                                    + self.pdbid + ' and not chain ' +
                                    query_chain + ')')
                            pymol.cmd.save(os.path.join(formatted_outdir,
                                        'combined', name + '.pdb'),
                                        'combined')
                            pymol.cmd.delete('combined')
                        if alignment['RMSD'] < best_rmsd:
                            best_rmsd = alignment['RMSD']
                            best_i = i
                    elif self.aligner=='align': 
                        if alignment[0] < 3.0:
                            if not os.path.exists(formatted_outdir):
                                os.mkdir(formatted_outdir)
                            if not os.path.exists(os.path.join(formatted_outdir,
                                'combined')):
                                os.mkdir(os.path.join(formatted_outdir,
                                    'combined'))
                            pymol.cmd.center('reference')
                            name = '{}_{}_length_{}_rmsd_{:.2f}'.format(self.pdbid,
                                    i, alignment[1],
                                    alignment[0])
                            pymol.cmd.save(os.path.join(formatted_outdir,
                                name + '.pse'))
                            pymol.cmd.create('combined', 'reference or ('
                                    + self.pdbid + ' and not chain ' +
                                    query_chain + ')')
                            pymol.cmd.save(os.path.join(formatted_outdir,
                                        'combined', name + '.pdb'),
                                        'combined')
                            pymol.cmd.delete('combined')
                        if alignment[0] < best_rmsd:
                            best_rmsd = alignment[0]
                            best_i = i
                except:
                    print('could not align {} with {}'.format(query_interface, reference_interface))

                i += 1
        print('Best alignment by RMSD: alignment {}'.format(best_i))

'''
For testing
'''
if __name__=='__main__':
    init()
    pymol.finish_launching(['pymol','-qc'])
    input_dir = 'test_inputs'

    pdbid = sys.argv[1]
    pose = pose_from_rcsb(pdbid, 'test_inputs')
    interface = PyInterface(pose)
    interface.find_interface()
    interface.find_patches()

    #query_resselectors = []

    #reference_selector = list_to_res_selector(reference_interface)
    
    #align_interfaces(interface, reference_interfaces, pdbid,
    #    reference_pdb)
    
    aligner='align'
    interface_aligner = PyMOLAligner(aligner, interface, reference_interfaces, pdbid,
            reference_pdb,
            output_dir=os.path.join('outputs','Q969X5'))
    interface_aligner.align_interfaces()
