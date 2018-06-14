'''
Elliot Williams
June 1st, 2018
ContactNetwork

This file aims to construct graphs readily available for analysis in the
networkit package -- allowing us to reason about the topological structures
associated with ribosome networks across different crystal structures.

The result of this code will be graphs saved in [some format] within
`data/CSN_graphs`, or a networkit object return.
'''

import pandas as pd
import glob
import networkit as nk
from math import sqrt
import numpy as np

pd.options.mode.chained_assignment = None

'''
This class will construct a network from the relevant positional file,
and either return it or output it to the `data/CSN_graphs` directory

Input:
    position_path - path to the node position file used to construct the graph
    csn_type  - the type of network constructed (in ["thresh", "shadow"])
    min_thresh    - the minimum distance bw nodes for contacts to be considered
    max_thresh    - the maximum distance bw nodes for contacts to be considered
    out_path      - Output graph save path (None if returning networkit object)
'''
class NetworkConstructor:
    def __init__(self, position_path, csn_type, min_thresh=0, max_thresh=8, out_path=None):
        assert csn_type in ["thresh"] # shadow networks are unimplemented

        print("Starting reading in the csv")

        # Reads in positional dataset
        df = pd.read_csv(position_path)
        df['node_id'] = range(df.shape[0])

        print("Read in the csv")

        if csn_type is "thresh":
            print("Starting network connecting code")

            G = nk.graph.Graph(n=df.shape[0])

            i = 0
            while df.shape[0] > 0:
                print("Currently dealing with node {}".format(i))

                curr_node = df['node_id'][i]
                cx, cy, cz = df['x'][i], df['y'][i], df['z'][i]

                df = df.iloc[1:] # Deletes first row of dataframe

                df_filt = df.query('@cx - @max_thresh < x')
                df_filt = df_filt.query('x < @cx + @max_thresh')

                df_filt = df_filt.query('@cy - @max_thresh < y')
                df_filt = df_filt.query('y < @cy + @max_thresh')

                df_filt = df_filt.query("@cz - @max_thresh < z")
                df_filt = df_filt.query("z < @cz + @max_thresh")


                df_filt['dx2']  = (df_filt['x'] - cx) ** 2
                df_filt['dy2']  = (df_filt['y'] - cy) ** 2
                df_filt['dz2']  = (df_filt['z'] - cz) ** 2
                df_filt['dist'] = np.sqrt(df_filt['dx2'] + df_filt['dy2'] + df_filt['dz2'])

                df_filt = df_filt.query("dist > @min_thresh and dist < @max_thresh")

                for neigh_node in df_filt['node_id']:
                    G.addEdge(curr_node, neigh_node)

                i += 1

            print("Outputting graph")
            nk.graphio.writeGraph(G, "data/CSN_graphs/5jup.graphviz", Format.GraphViz)




if __name__ == "__main__":
    NetworkConstructor("data/positions/5jup.csv", "thresh")
