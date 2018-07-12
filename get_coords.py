'''
Elliot Williams
July 3rd, 2018
ContactNetwork

This file aims to construct a dataframe of points based on the PDB input files,
corresponding the the CA position in residues and C1' position in nucleotides.

The result of this code will be a csv file output to

'''

import pymol
import time
import numpy as np
import sys
import pandas as pd
import glob

# TODO: Enable option to download PDB structures as need be followed by deletion,
#       autodownload structures, or just use structures in directory only

# TODO: Figure out how to map PyMol chain names to PDB chain names
#       (they're not the same ...)


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
        type         : Contact network type ("aa", "nt", "li", or "all"),
        min_thresh   : Minimum distance threshold for 'contact' in angstroms,
        max_thresh   : Maximum distance threshold for 'contact' in angstroms,
        grain        : Whether or not the contact network should be based on
                       alpha carbon distances, or all-atom distances
                       ("allatom", or "ca")

    '''
    def __init__(self, struct_names, type="all", grain="allatom"):

        # The cmd object is equivalent to the PyMol command line, and refers to the
        # object which has functions corresponding to PyMol's API (poorly documented)
        # See: https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha
        self.cmd = pymol.cmd
        cmd = pymol.cmd
        self.stored = pymol.stored
        stored = pymol.stored

        # Catches common argument errors
        assert type in ["aa", "nt", "li", "all"]
        assert grain in ["allatom","ca"]

        # Generates coordinate-generating PyMol Selection
        if grain is "allatom":
            self.nt_sel = "resn a+t+c+g+u and !(hetatm)"
            self.aa_sel = "byres name ca and !(hetatm)"
        else:
            self.nt_sel = "name C1' and resn a+t+c+g+u and !(hetatm)'"
            self.aa_sel = "name ca and !(hetatm)"

        # Note: This is guaranteed to contain all ligand atoms
        # TODO: Ask KMT whether or not this is reasonable behaviour
        self.li_sel = "not((resn a+t+c+g+u) | (byres name ca)) or hetatm"

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

            aa_x, aa_y, aa_z, aa_chains, aa_index, aa_resi, aa_resn = [tuple([])]*7
            nt_x, nt_y, nt_z, nt_chains, nt_index, nt_resi, nt_resn = [tuple([])]*7
            li_x, li_y, li_z, li_chains, li_index, li_resi, li_resn = [tuple([])]*7

            if type in ["aa", "all"]:
                # Gets Euclidean position, chain name, PyMol index for amino acids
                aa_xyz, aa_chains, aa_index, aa_resi, aa_resn = self.get_aa_data()
                aa_x, aa_y, aa_z = zip(*aa_xyz)
                chain_resi = list(map(lambda x: x[0] + str(x[1]), zip(aa_chains, aa_resi)))
                num_aa = len(set(chain_resi))
                print("There are {} amino acids in the network".format(num_aa))
                print("There are {} amino acid atoms".format(len(set(aa_index))))

            if type in ["nt", "all"]:
                # Gets Euclidean position, chain name, PyMol index for nucleotides
                nt_xyz, nt_chains, nt_index, nt_resi, nt_resn = self.get_nt_data()
                nt_x, nt_y, nt_z = zip(*nt_xyz)
                chain_resi = list(map(lambda x: x[0] + str(x[1]), zip(nt_chains, nt_resi)))
                num_nt = len(set(chain_resi))
                print("There are {} nucleotides in the network".format(num_nt))
                print("There are {} nucleotide atoms".format(len(set(nt_index))))

            if type in ["li", "all"]:
                # Gets Euclidean position, chain name, PyMol index for nucleotides
                li_xyz, li_chains, li_index, li_resi, li_resn = self.get_li_data()
                li_x, li_y, li_z = zip(*li_xyz)
                chain_resi = list(map(lambda x: x[0] + str(x[1]), zip(li_chains, li_resi)))
                num_li = len(set(chain_resi))
                print("There are {} ligands in the network".format(num_li))
                print("There are {} ligand atoms".format(len(set(li_index))))

            # Constructs DataFrame from the above PyMol gleaned data
            df = pd.DataFrame({'x': nt_x + aa_x + li_x,
                               'y': nt_y + aa_y + li_y,
                               'z': nt_z + aa_z + li_z,
                               'chain': list(nt_chains) + list(aa_chains) + list(li_chains),
                               'index': list(nt_index)  + list(aa_index)  + list(li_index),
                               'resi' : list(nt_resi)   + list(aa_resi)   + list(li_resi) ,
                               'resn' : list(nt_resn)   + list(aa_resn)   + list(li_resn),
                               'type': ['nucleotide']*len(nt_x) +
                                       ['amino acid']*len(aa_x) +
                                       ['ligand']*len(li_x)
                               })

            print("Writing positional CSV to `data/positions/`")
            df.to_csv("data/positions/{}.csv".format(actual_name[:-4]))

            # Deletes the molecule from PyMol, preventing memory issues
            cmd.delete("all")


    def get_aa_data(self):

        cmd = self.cmd
        stored = self.stored

        # Selects all amino acid atoms (ie all residues containing alpha carbon)
        cmd.select("aa", self.aa_sel)

        # Gets Euclidean position of each atom in aa selection
        stored.aa_xyz = []
        cmd.iterate_state(1,"aa","stored.aa_xyz.append([x,y,z])")

        # Gets the chain associated with each atom in aa selection
        stored.aa_chains = []
        cmd.iterate_state(1, "aa", "stored.aa_chains.append(segi)")

        # Gets the PyMol unique index associated with each atom in aa selection
        stored.aa_index = []
        cmd.iterate_state(1, "aa", "stored.aa_index.append(index)")

        # Gets the residue index associated with each atom in aa
        stored.aa_resi = []
        cmd.iterate_state(1, "aa", "stored.aa_resi.append(resi)")

        # Gets the amino acid type associated w each atom in aa_ca
        stored.aa_resn = []
        cmd.iterate_state(1, "aa", "stored.aa_resn.append(resn)")

        return (stored.aa_xyz, stored.aa_chains, stored.aa_index,
               stored.aa_resi, stored.aa_resn)



    def get_nt_data(self):

        cmd = self.cmd
        stored = self.stored

        # Selects all atoms associated with each nucleotide
        cmd.select("nucleo", self.nt_sel)

        # Gets Euclidean position of each atom
        stored.nt_xyz = []
        cmd.iterate_state(1,"nucleo","stored.nt_xyz.append([x,y,z])")
        nt_x, nt_y, nt_z = zip(*stored.nt_xyz)

        # Gets the chain associated with each atom
        stored.nt_chains = []
        cmd.iterate_state(1, "nucleo", "stored.nt_chains.append(segi)")

        # Gets the PyMol unique index associated with each atom
        stored.nt_index = []
        cmd.iterate_state(1, "nucleo", "stored.nt_index.append(index)")

        # Gets the residue index associated with each atom
        stored.nt_resi = []
        cmd.iterate_state(1, "nucleo", "stored.nt_resi.append(resi)")

        # Gets the base type associated with each atom
        stored.nt_resn = []
        cmd.iterate_state(1, "nucleo", "stored.nt_resn.append(resn)")

        return (stored.nt_xyz, stored.nt_chains, stored.nt_index,
               stored.nt_resi, stored.nt_resn)

    def get_li_data(self):

        cmd = self.cmd
        stored = self.stored

        # Selects all atoms associated with each nucleotide
        cmd.select("ligands", self.li_sel)

        # Gets Euclidean position of each atom
        stored.li_xyz = []
        cmd.iterate_state(1,"ligands","stored.li_xyz.append([x,y,z])")
        li_x, li_y, li_z = zip(*stored.li_xyz)

        # Gets the chain associated with each atom
        stored.li_chains = []
        cmd.iterate_state(1, "ligands", "stored.li_chains.append(segi)")

        # Gets the PyMol unique index associated with each atom
        stored.li_index = []
        cmd.iterate_state(1, "ligands", "stored.li_index.append(index)")

        # Gets the residue index associated with each atom
        stored.li_resi = []
        cmd.iterate_state(1, "ligands", "stored.li_resi.append(resi)")

        # Gets the base type associated with each atom
        stored.li_resn = []
        cmd.iterate_state(1, "ligands", "stored.li_resn.append(resn)")

        return (stored.li_xyz, stored.li_chains, stored.li_index,
               stored.li_resi, stored.li_resn)

if __name__ == "__main__":

    if len(sys.argv) == 2:
        type = sys.argv[1]
    else:
        type = "all"

    struct_names = glob.glob("data/pdb_structs/*.pdb")
    print(">> Reading in {} structures".format(len(struct_names)))
    CoordConstruct(struct_names, type=type, grain="allatom")
