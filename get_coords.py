'''
Elliot Williams
June 1st, 2018
ContactNetwork
'''

import pymol
import time
import numpy as np
import sys
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
                       (these should be located within data/ ),
        type         : Contact network type ("aa", "nt", or "both"),
        min_thresh   : Minimum distance threshold for 'contact' in angstroms,
        max_thresh   : Maximum distance threshold for 'contact' in angstroms,
        valid_chains : List of chains that we actually want to analyze
                       (ie one might want to exclude tRNA, mRNA chains)

    '''
    def __init__(self, struct_names, type="both", min_thresh=None,
                 max_thresh=8, valid_chains=None):
        pass # Unimplemented


# Easy way to launch with/without GUI (fix later)
if len(sys.argv) == 2:
    # Launches PyMol quietly with GUI
    pymol.finish_launching(['pymol', '-q'])
else:
    # Launches headless PyMol quietly
    pymol.finish_launching(['pymol', '-qc'])


# The cmd object is equivalent to the PyMol command line, and refers to the
# object which has functions corresponding to PyMol's API (poorly documented)
# See: https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha
cmd = pymol.cmd
stored = pymol.stored

struct_name = "5jup.pdb"

# Opens test structure (in this case, PDB code 5JUP Korostelev Ribosome struct)
cmd.load("pdb_structs/{}".format(struct_name))

chain_names = cmd.get_chains()

print("Structure {} has chains {}".format(struct_name, chain_names))

# Hides everything (for testing visualization -- TODO: Delete later)
cmd.hide("everything","all")

# Selects the alpha carbon associated with each amino acid
cmd.select("aa_ca", "name ca")
cmd.show("dots", "aa_ca") # TODO: Delete later, for testing visualization
cmd.color("red", "aa_ca") # TODO: Delete later, for testing visualization

# Selects the C1' associated with each nucleotide
cmd.select("nucleo_C1", "name C1' and resn a+c+g+u")
cmd.show("dots", "nucleo_C1") # TODO: Delete later, for testing visualization
cmd.color("blue", "nucleo_C1") # TODO: Delete later, for testing visualization

cmd.deselect() # TODO: Delete later, for testing visualization

# Gets Euclidean position of each alpha carbon in aa_ca
stored.aa_xyz = []
cmd.iterate_state(1,"aa_ca","stored.aa_xyz.append([x,y,z])")

# Gets the chain associated with each alpha carbon in aa_ca
stored.aa_chains = []
cmd.iterate_state(1, "aa_ca", "stored.aa_chains.append(chain)")

# Gets the PyMol unique index associated with each alpha carbon in aa_ca
stored.aa_index = []
cmd.iterate_state(1, "aa_ca", "stored.aa_index.append(index)")
# print(len(stored.aa_index))

# Gets Euclidean position of each C1' within nucleo_C1
stored.nt_xyz = []
cmd.iterate_state(1,"nucleo_C1","stored.nt_xyz.append([x,y,z])")

# Gets the chain associated with each alpha carbon in aa_ca
stored.nt_chains = []
cmd.iterate_state(1, "nucleo_C1", "stored.nt_chains.append(chain)")

# Gets the PyMol unique index associated with each alpha carbon in aa_ca
stored.nt_index = []
cmd.iterate_state(1, "nucleo_C1", "stored.nt_index.append(index)")

print("There are {} nucleotides in the network".format(len(stored.nt_index)))
print("There are {} amino acids in the network".format(len(stored.aa_index)))

# TODO: Add code which combines the data (preferably into a dataframe of some
#       sort which then can be manipulated/returned in additional code)

# TODO: See if it's faster to calculate distances yourself, or to use pymol
#       distance commands (I think the latter will be orders of magnitude faster)



# TODO: Add edge construction code in either networkit or graph-tool
#       (Requires benchmarking them first)

# exit()
