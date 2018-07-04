'''
Elliot Williams
June 1st, 2018
ContactNetwork

This file aims to construct a dataframe of points based on the PDB input files,
corresponding the the CA position in residues and C1' position in nucleotides.

The result of this code will be a csv file output to `data/positions/`
'''

import pymol
import time
import numpy as np
import sys
import pandas as pd
import glob
# import networkit as nk

# TODO: Enable option to download PDB structures as need be followed by deletion,
#       autodownload structures, or just use structures in directory only

# TODO: Figure out how to map PyMol chain names to PDB chain names
#       (they're not the same ...)

# TODO: Wrap this all within a class definition so it's both pretty and usable

# TODO: Add an option for changing what the nucleotide atom of choice is


'''
This class defines and contains the PyMol interface code with which we obtain
the dataframe containing the residue index, Euclidean location, and
corresponding chain name of the specific amino acid / nucleotide atoms we are
using to construct a contact network.
'''
class CoordConstruct:

    '''
    Class initialization function
    Input:
        struct_names : String list of structures to construct networks for
                       (these should be located within pdb_structs/ ),
        type         : Contact network type ("aa", "nt", or "both"),
        min_thresh   : Minimum distance threshold for 'contact' in angstroms,
        max_thresh   : Maximum distance threshold for 'contact' in angstroms,
        valid_chains : List of chains that we actually want to analyze
                       (ie one might want to exclude tRNA, mRNA chains)

    '''
    def __init__(self, struct_names, type="both", min_thresh=None,
                 max_thresh=8, valid_chains=None):

        # The cmd object is equivalent to the PyMol command line, and refers to the
        # object which has functions corresponding to PyMol's API (poorly documented)
        # See: https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha
        self.cmd = pymol.cmd
        cmd = pymol.cmd
        self.stored = pymol.stored
        stored = pymol.stored

        for struct_name in struct_names:

            actual_name = struct_name.split('/')[-1]
            print(actual_name)

            print("> Opening {}".format(struct_name))

            # Launches headless PyMol quietly
            pymol.finish_launching(['pymol', '-qc'])

            # Opens structure in headless PyMol
            cmd.load(struct_name)

            # Assigns default values to amino acid, nucleotide parameters, to
            # deal with cases in which we don't want to include whole network
            # data.

            aa_x, aa_y, aa_z, aa_chains, aa_index, aa_resi, aa_resn = [ [] ]*7
            nt_x, nt_y, nt_z, nt_chains, nt_index, nt_resi, nt_resn = [ [] ]*7

            if type in ["aa", "both"]:
                # Gets Euclidean position, chain name, PyMol index for amino acids
                aa_xyz, aa_chains, aa_index, aa_resi, aa_resn = self.get_aa_data()
                aa_x, aa_y, aa_z = zip(*aa_xyz)
                print("There are {} amino acids in the network".format(len(aa_index)))

            if type in ["nt", "both"]:
                # Gets Euclidean position, chain name, PyMol index for nucleotides
                nt_xyz, nt_chains, nt_index, nt_resi, nt_resn = self.get_nt_data()
                nt_x, nt_y, nt_z = zip(*nt_xyz)
                print("There are {} nucleotides in the network".format(len(nt_index)))

            # Constructs DataFrame from the above PyMol gleaned data
            df = pd.DataFrame({'x': nt_x + aa_x, 'y': nt_y + aa_y, 'z': nt_z + aa_y,
                                 'chain': nt_chains + aa_chains,
                                 'index': nt_index  + aa_index,
                                 'resi' : nt_resi   + aa_resi ,
                                 'resn' : nt_resn   + aa_resn ,
                                 'type': ['nucleotide']*len(nt_x)+['amino acid']*len(aa_x)})

            print("Writing positional CSV to `data/positions/`")
            df.to_csv("data/positions/{}.csv".format(actual_name[:-4]))

            # Deletes the molecule from PyMol, preventing memory issues
            cmd.delete("all")


    def get_aa_data(self):

        cmd = self.cmd
        stored = self.stored

        # Selects the alpha carbon associated with each amino acid
        cmd.select("aa_ca", "name ca")

        # Gets Euclidean position of each alpha carbon in aa_ca
        stored.aa_xyz = []
        cmd.iterate_state(1,"aa_ca","stored.aa_xyz.append([x,y,z])")

        # Gets the chain associated with each alpha carbon in aa_ca
        stored.aa_chains = []
        cmd.iterate_state(1, "aa_ca", "stored.aa_chains.append(chain)")

        # Gets the PyMol unique index associated with each alpha carbon in aa_ca
        stored.aa_index = []
        cmd.iterate_state(1, "aa_ca", "stored.aa_index.append(index)")

        # Gets the residue index associated with each alpha carbon in aa_ca
        stored.aa_resi = []
        cmd.iterate_state(1, "aa_ca", "stored.aa_resi.append(resi)")

        # Gets the amino acid type associated w each alpha carbon in aa_ca
        stored.aa_resn = []
        cmd.iterate_state(1, "aa_ca", "stored.aa_resn.append(resn)")

        return (stored.aa_xyz, stored.aa_chains, stored.aa_index,
               stored.aa_resi, stored.aa_resn)



    def get_nt_data(self):

        cmd = self.cmd
        stored = self.stored

        # Selects the C1' associated with each nucleotide
        cmd.select("nucleo_C1", "name C1' and resn a+c+g+u")

        # Gets Euclidean position of each C1' within nucleo_C1
        stored.nt_xyz = []
        cmd.iterate_state(1,"nucleo_C1","stored.nt_xyz.append([x,y,z])")
        nt_x, nt_y, nt_z = zip(*stored.nt_xyz)

        # Gets the chain associated with each C1' in nucleo_C1
        stored.nt_chains = []
        cmd.iterate_state(1, "nucleo_C1", "stored.nt_chains.append(chain)")

        # Gets the PyMol unique index associated with each C1' in nucleo_C1
        stored.nt_index = []
        cmd.iterate_state(1, "nucleo_C1", "stored.nt_index.append(index)")

        # Gets the residue index associated with each C1' in nucleo_C1
        stored.nt_resi = []
        cmd.iterate_state(1, "nucleo_C1", "stored.nt_resi.append(resi)")

        # Gets the base type associated with each C1' in nucleo_C1
        stored.nt_resn = []
        cmd.iterate_state(1, "nucleo_C1", "stored.nt_resn.append(resn)")

        return (stored.nt_xyz, stored.nt_chains, stored.nt_index,
               stored.nt_resi, stored.nt_resn)

if __name__ == "__main__":

    struct_names = glob.glob("data/pdb_structs/*.pdb")
    print(">> Reading in {} structures".format(len(struct_names)))
    CoordConstruct(struct_names, type="both", min_thresh=None,
                 max_thresh=8, valid_chains=None)


# TODO: See if it's faster to calculate distances yourself, or to use pymol
#       distance commands (I think the latter will be orders of magnitude faster)

# TODO: Add edge construction code in either networkit or graph-tool
#       (Requires benchmarking them first)
