'''
Elliot Williams
July 4th, 2018
ContactNetwork -- Network Construction

This file aims to construct graphs readily available for analysis in the
networkit package -- allowing us to reason about the topological structures
associated with ribosome networks across different crystal structures.

The result of this code will be graphs saved in METIS within
`data/CSN_graphs`, or a networkit object return
(not recommended -- this takes a while for ribosome structures).
'''

import pandas as pd
import glob
import sys
import networkit as nk
from math import sqrt
import numpy as np
import os

pd.options.mode.chained_assignment = None

'''
This class will construct a network from the relevant positional file,
and either return it or output it to the `data/CSN_graphs` directory

Input:
    position_path - path to the node position file used to construct the graph
    csn_type      - type of network constructed (in ["thresh", "shadow"])
    min_thresh    - the minimum distance bw nodes for contacts to be considered
    max_thresh    - the maximum distance bw nodes for contacts to be considered
    out_path      - Output graph save path (None if returning networkit object)
    verbose       - whether or not the constructor should be verbose
'''
class NetworkConstructor:
    def __init__(self, position_path, csn_type, min_thresh=0, max_thresh=8,
                 out_path=None, verbose = True):
        assert csn_type in ["thresh"] # shadow networks are unimplemented

        if verbose:
            print("Reading in {}".format(position_path.split("/")[-1]))

        # Reads in positional dataset
        df = pd.read_csv(position_path)
        df['node_id'] = range(df.shape[0])

        if csn_type is "thresh":
            print("Starting network connecting code")

            # Creates a new networkit graph with the correct number of nodes
            G = nk.graph.Graph(n=df.shape[0])

            i = 0
            resid = None
            while df.shape[0] > 0:

                # Gets the positional and index info for current node
                curr_node  = df['resi'][i]
                curr_chain = df['chain']
                cx, cy, cz = df['x'][i], df['y'][i], df['z'][i]

                if verbose:
                    if resid != curr_node:
                        print(">>Currently dealing with residue {}".format(curr_node))
                        resid = curr_node

                df = df.iloc[1:] # Deletes first row of dataframe

                # Filters out nodes that can't be within contact range
                df_filt = df.query('@cx - @max_thresh < x')
                df_filt = df_filt.query('x < @cx + @max_thresh')

                df_filt = df_filt.query('@cy - @max_thresh < y')
                df_filt = df_filt.query('y < @cy + @max_thresh')

                df_filt = df_filt.query("@cz - @max_thresh < z")
                df_filt = df_filt.query("z < @cz + @max_thresh")

                # And filters out self links
                df_filt = df_filt.query("!(resi == @curr_node and chain == @curr_chain)")

                # Does actual distance calculation, for the remaining nodes
                df_filt['dx2']  = (df_filt['x'] - cx) ** 2
                df_filt['dy2']  = (df_filt['y'] - cy) ** 2
                df_filt['dz2']  = (df_filt['z'] - cz) ** 2
                df_filt['dist'] = np.sqrt(df_filt['dx2'] + df_filt['dy2'] + df_filt['dz2'])

                # Selects only those nodes within contact range
                df_filt = df_filt.query("dist > @min_thresh and dist < @max_thresh")

                # Adds edges between the contacting nodes
                for neigh_node in df_filt['resi']:
                    if not G.hasEdge(curr_node, neigh_node):
                        G.addEdge(curr_node, neigh_node)

                i += 1

            if out_file is not None:
                if verbose:
                    print("Writing graph to {}".format(out_file))
                nk.graphio.writeGraph(G, out_file, nk.Format.METIS)
            else:
                if verbose:
                    print("Returning graph object")
                return G

if __name__ == "__main__":

    # Should take the distance cutoffs as arguments to the Python file
    assert len(sys.argv) == 3
    min_thresh = int(sys.argv[1])
    max_thresh = int(sys.argv[2])

    print("-- Welcome to the Contact Network Constructor --")
    print("(Constructing basic network with min={}, max={})".format(
        min_thresh, max_thresh))
    print("\nWe hope your stay with us is pleasant and enjoyable")

    # Iterates over every position CSV within the folder
    for pos_file in glob.glob('data/positions/*.csv'):
        out_file = "data/CSN_graphs/{}_{}A_{}.graph".format(
            pos_file.split("/")[-1][:-4], max_thresh, "thresh")

        if not os.path.isfile(out_file):
            NetworkConstructor(pos_file, "thresh", 0, max_thresh, out_file, True)
        else:
            print("Already Exists: {}".format(out_file))
