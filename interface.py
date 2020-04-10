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


def empty_interface_dataframe():
    dataframe = pd.DataFrame(
            columns=[
                'pymol_chain',
                'rosetta_chain',
                'pymol_tarchain',
                'rosetta_tarchain',
                'rosetta_full_interface',
                'rosetta_target_interface',
                'pymol_full_interface',
                'pymol_target_interface',
                'rosetta_patch',
                'pymol_patch',
                'rosetta_target_patch',
                'pymol_target_patch',
                'pymol_reference_interface',
                'rosetta_reference_interface'
                ]
            )
    return dataframe


class PyInterface():
    '''
    Class to handle interface definitions.
    '''
    def __init__(self, pose, reference_pose):
        self.pose = pose
        self.reference_pose = reference_pose
        scorefxn = create_score_function('ref2015')
        scorefxn(self.pose)
        self.cutoff_ = 8.0
        self.patch_residues = 8
        self.patch_cutoff = 30
        self.interfaces = []
        self.pdb_interfaces = []
        self.reference_interfaces = []
        self.dataframe = empty_interface_dataframe()

    @property
    def cutoff(self):
        return self.cutoff_
    
    @cutoff.setter
    def cutoff(self, cutoff):
        self.cutoff_ = cutoff

    def set_dataframe(self, df):
        '''Load a previously generated dataframe'''
        self.dataframe = df
    
    def find_interface(self):
        for jump in range(1, self.pose.fold_tree().num_jump() + 1):
            interface = Interface(jump)
            interface.distance(self.cutoff_)
            interface.calculate(self.pose)
            self.interfaces.append(interface)
        self.set_pdb_interface()

    def set_pdb_interface(self, residue_list, chain='Z'):
        '''
        Don't use this at the moment; the idea was to use the same class
        to define the reference interface, but that may not work since
        this is a monomer interface. Keping this code for now just in
        case.
        '''
        pdbinfo = self.pose.pdb_info()
        self.pdb_interfaces = []
        self.interfaces = []
        for resi in residue_list:
            self.pdb_interfaces.append(' '.join([str(resi), chain]))
            self.interfaces.append(pdbinfo.pdb2pose(chain, resi))

    def add_pdb_interface(self, residue_list, chain='Z'):
        '''
        Don't use this at the moment; the idea was to use the same class
        to define the reference interface, but that may not work since
        this is a monomer interface. Keping this code for now just in
        case.
        '''
        pdbinfo = self.pose.pdb_info()
        temp_interface = []
        temp_pdb_interface = []
        for resi in residue_list:
            temp_pdb_interface.append(' '.join([str(resi), chain]))
            temp_interface.append(pdbinfo.pdb2pose(chain, resi))
        self.interfaces.append(temp_interface)
        self.pdb_interfaces.append(temp_pdb_interface)

    def set_pdb_interface(self):
        self.pdb_interfaces = []
        for interface_obj in self.interfaces:
            for side in interface_obj.pair_list():
                self.pdb_interfaces.append(reslist_to_pdb_numbers(side,
                    self.pose))

    def set_reference_interfaces(self, list_of_lists):
        '''
        Define reference interface
        '''
        self.reference_interfaces = list_of_lists

    def add_reference_interface(self, list_of_pdb_residues):
        '''
        Add to reference interface definition
        '''
        self.reference_interfaces.append(list_of_pdb_residus)

    def find_patches(self):
        '''
        Build dataframe of patches for all interfaces
        '''
        #df_list = []
        for interface in self.interfaces:
            # Get PDB info from interface
            pair_list = vector1_to_python_list(interface.pair_list())
            sideA = vector1_to_python_list(pair_list[0])
            sideA_pymol = reslist_to_pdb_numbers(sideA, self.pose)
            sideB = vector1_to_python_list(pair_list[1])
            sideB_pymol = reslist_to_pdb_numbers(sideB, self.pose)
            print('SIDEA',sideA)
            print('SIDEB',sideB)

            chainA = self.pose.pdb_info().pose2pdb(sideA[0]).split(' ')[1]
            rosetta_chainA = self.pose.chain(sideA[0])
            chainB = self.pose.pdb_info().pose2pdb(sideB[0]).split(' ')[1]
            rosetta_chainB = self.pose.chain(sideB[0])

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
                if patch:

                    # Get "partner patch"
                    interacting_residues = []
                    for resi in patch:
                        partners = vector1_to_python_list(vector1_to_python_list(
                                interface.contact_list())[resi - 1])
                        interacting_residues += partners
                    # No need to add a patch twice
                    if not patch in self.dataframe['rosetta_patch'].values:
                        for reslist in self.reference_interfaces:
                            row = {
                                    'pymol_chain': chainA,
                                    'rosetta_chain': rosetta_chainA,
                                    'pymol_tarchain': chainB,
                                    'rosetta_tarchain': rosetta_chainB,
                                    'rosetta_full_interface': sideA,
                                    'rosetta_target_interface': sideB,
                                    'pymol_full_interface': sideA_pymol,
                                    'pymol_target_interface': sideB_pymol,
                                    'rosetta_patch': patch,
                                    'pymol_patch': reslist_to_pdb_numbers(patch,
                                        self.pose),
                                    'rosetta_target_patch': interacting_residues,
                                    'pymol_target_patch': set(reslist_to_pdb_numbers(interacting_residues,
                                        self.pose)),
                                    'pymol_reference_interface':
                                    int_list_to_pdb_numbers(reslist),
                                    'rosetta_reference_interface':
                                    rosetta_numbers_from_pdb(reslist,
                                        self.reference_pose),
                                    }
                            self.dataframe = self.dataframe.append(row, ignore_index=True)

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

                if patch:

                    # Get "partner patch"
                    interacting_residues = []
                    for resi in patch:
                        partners = vector1_to_python_list(vector1_to_python_list(
                                interface.contact_list())[resi - 1])
                        interacting_residues += partners
                    # No need to add a patch twice
                    if not patch in self.dataframe['rosetta_patch'].values:
                        for reslist in self.reference_interfaces:
                            row = {
                                    'pymol_chain': chainB,
                                    'rosetta_chain': rosetta_chainB,
                                    'pymol_tarchain': chainA,
                                    'rosetta_tarchain': rosetta_chainA,
                                    'rosetta_full_interface': sideB,
                                    'rosetta_target_interface': sideA,
                                    'pymol_full_interface': sideB_pymol,
                                    'pymol_target_interface': sideA_pymol,
                                    'rosetta_patch': patch,
                                    'pymol_patch': reslist_to_pdb_numbers(patch,
                                        self.pose),
                                    'rosetta_target_patch': interacting_residues,
                                    'pymol_target_patch': set(reslist_to_pdb_numbers(interacting_residues,
                                        self.pose)),
                                    'pymol_reference_interface':
                                    int_list_to_pdb_numbers(reslist),
                                    'rosetta_reference_interface':
                                    rosetta_numbers_from_pdb(reslist,
                                        self.reference_pose),
                                    }
                            self.dataframe = self.dataframe.append(row, ignore_index=True)
        #self.dataframe = pd.DataFrame(df_list)


class PyMOLAligner(object):
    def __init__(self, aligner, interface,
            pdbid, reference_pdb, output_dir='.',
            input_dir='test_inputs', window=3, cycles=0):
        # String, either 'cealign' or 'align'
        self.aligner=aligner
        # PyInterface object
        self.interface = interface
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

        # Make reference chain Z so that we can combine it with the
        # aligned structure later (THIS MEANS IT ONLY WORKS WITH
        # MONOMER REFERENCE STRUCTURES)
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

    def align_patches(self):
        '''
        Align all interface patches defined in reference and query proteins.
        '''
        i = 0
        best_i = 0
        best_rmsd = 999
        formatted_outdir = os.path.join(self.output_dir,
                '{}_{}'.format(self.pdbid, self.aligner))
        # Leave option to align entire interface instead of patches
        # To do: add option to align reference interface by patches
        for idx, row in self.interface.dataframe.iterrows():
            query_interface = row['pymol_patch']
            if len(query_interface) == 0:
                # Nothing to be done, no interface
                continue
            query_chain = row['pymol_chain']
            reference_interface = row['pymol_reference_interface']
            self.interface.dataframe.at[idx,'alignment'] = \
                    self.aligner
            try:
                alignment = self.align(reference_interface,
                        query_interface)
                print(i)
                print(alignment)
                if self.aligner=='cealign':
                    self.interface.dataframe.at[idx,'cealign_rmsd'] =\
                            alignment['RMSD']
                    self.interface.dataframe.at[idx,'cealign_len'] =\
                            alignment['alignment_length']
                    if 0.0 < alignment['RMSD'] < 3.0:
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
                        self.interface.dataframe.at[idx,'cealign_combined_pdb_path']=\
                                os.path.join(formatted_outdir, 'combined',
                                name + '.pdb')
                        pymol.cmd.delete('combined')
                    if alignment['RMSD'] < best_rmsd:
                        best_rmsd = alignment['RMSD']
                        best_i = i
                elif self.aligner=='align': 
                    self.interface.dataframe.at[idx,'align_rmsd'] =\
                            alignment[0]
                    self.interface.dataframe.at[idx,'align_atoms'] =\
                            alignment[1]
                    if 0.0 < alignment[0] < 3.0:
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
                        self.interface.dataframe.at[idx,'align_combined_pdb_path']=\
                                os.path.join(formatted_outdir, 'combined',
                                name + '.pdb')
                        pymol.cmd.delete('combined')
                    if alignment[0] < best_rmsd:
                        best_rmsd = alignment[0]
                        best_i = i
            except:
                print('could not align {} with {}'.format(query_interface, reference_interface))

            i += 1
        print('Best alignment by RMSD: alignment {}'.format(best_i))
        num = 0
        if not os.path.exists(formatted_outdir):
            os.mkdir(formatted_outdir)
        df_file = os.path.join(formatted_outdir, 'patches.pkl')
        while os.path.exists(df_file):
            num += 1
            df_file = os.path.join(formatted_outdir,
                    'patches_{}.pkl'.format(num))
        print(df_file)
        self.interface.dataframe.to_pickle(df_file)

    def align_interfaces(self):
        '''
        Align all interfaces defined in reference and query proteins.
        '''
        i = 0
        best_i = 0
        best_rmsd = 999
        formatted_outdir = os.path.join(self.output_dir,
                '{}_{}'.format(self.pdbid, self.aligner))
        # Leave option to align entire interface instead of patches
        # To do: add option to align reference interface by patches
        query_iterator = self.interface.dataframe['pymol_full_interface'].unique()
        reference_interfaces = self.interface.DataFrame['pymol_reference_interface'].unique()
        for query_interface in query_iterator:
            if len(query_interface) == 0:
                # Nothing to be done, no interface
                continue
            query_chain = query_interface[0].split(' ')[1]
            for reference_interface in reference_interfaces:
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
    reference_pose = pose_from_file(reference_pdb)
    interface = PyInterface(pose, reference_pose)
    interface.find_interface()
    interface.set_reference_interfaces(reference_interfaces)
    interface.find_patches()
    print(interface.dataframe)

    #query_resselectors = []

    #reference_selector = list_to_res_selector(reference_interface)
    
    #align_interfaces(interface, reference_interfaces, pdbid,
    #    reference_pdb)
    
    aligner='align'
    interface_aligner = PyMOLAligner(aligner, interface, pdbid,
            reference_pdb,
            output_dir=os.path.join('outputs','Q969X5'))
    interface_aligner.align_patches()
    print(interface_aligner.interface.dataframe)
