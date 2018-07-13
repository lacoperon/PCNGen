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
import csv
import ray
import random

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

@ray.remote
def writeEdgeList(df, NodeDict, min_thresh, max_thresh, out_file, csn_type="thresh", verbose=True, dependency=None):
    pd.options.mode.chained_assignment = None
    print("Writing edgelist")
    i = 0
    resid = None
    rand_name = out_file + str(random.getrandbits(100)) + ".tmp" #TODO: Make this invisible
    edge_set = set()

    while df.shape[0] > 0:
        # Gets the positional and index info for current node
        try:
            curr_node  = df['node_id'].iloc[i]
        except IndexError:
            break
        node_num   = NodeDict[curr_node]
        cx, cy, cz = df['x'].iloc[i], df['y'].iloc[i], df['z'].iloc[i]

        if verbose:
            if resid != curr_node:
                print(">>At residue {}".format(curr_node))
                resid = curr_node

        df = df.iloc[1:] # Deletes first row of dataframe

        # Filters out nodes that can't be within contact range
        df_filt = df.query('@cx - @max_thresh < x')
        df_filt = df_filt.query('x < @cx + @max_thresh')

        df_filt = df_filt.query('@cy - @max_thresh < y')
        df_filt = df_filt.query('y < @cy + @max_thresh')

        df_filt = df_filt.query("@cz - @max_thresh < z")
        df_filt = df_filt.query("z < @cz + @max_thresh")

        # # And filters out self links
        df_filt = df_filt.query("~(node_id == @resid)")

        # Does actual distance calculation, for the remaining nodes
        df_filt['dx2']  = (df_filt['x'] - cx) ** 2
        df_filt['dy2']  = (df_filt['y'] - cy) ** 2
        df_filt['dz2']  = (df_filt['z'] - cz) ** 2
        df_filt['dist'] = np.sqrt(df_filt['dx2'] + df_filt['dy2'] + df_filt['dz2'])

        # Selects only those nodes within contact range
        df_filt = df_filt.query("dist > @min_thresh and dist < @max_thresh")

        f= open(rand_name,"a")

        # Writes edges between the contacting nodes
        for neigh_node in df_filt['node_id']:
            neigh_node_num = NodeDict[neigh_node]
            edge_string = "{}\t{}\n".format(node_num, neigh_node_num)
            if edge_string not in edge_set:
                f.write(edge_string)
                edge_set.add(edge_string)
                if verbose:
                    print("Added edge bw {} and {}".format(node_num, neigh_node_num))
        i += 1
    return "x"

@ray.remote
def createGraph(out_file, NodeDict, verbose = True, dependency=None):
    # Creates a new networkit graph with the correct number of nodes
    G = nk.graph.Graph(n=len(NodeDict))

    # Reads in adjacency files
    for edge_file in glob.glob("{}*.tmp".format(out_file)):
        f = open(edge_file, "r")

        # Converts edge file to tuple list
        edges = map(x.split("/t") for x in f.readlines())
        edges = map(lambda x: (int(x[0]), int(x[1])))

        # Adds edges to graph if they're not there already
        for edge in edges:
            if not G.hasEdge(edge[0],edge[1]):
                G.addEdge(edge[0],edge[1])

    if out_file is not None:
        if verbose:
            print("Writing graph to {}".format(out_file))
        nk.graphio.writeGraph(G, out_file, nk.Format.GML)
    else:
        if verbose:
            print("Returning graph object")
        return G

class NetworkConstructor:

    def makeNodeDict(self, unique_nodes, out_file):
        # Creates a mapping between all residues and node number in networkit, and outputs that csv
        NodeDict = dict(zip(unique_nodes, range(0, len(unique_nodes))))
        with open('{}{}.csv'.format(out_file.split(".")[0], out_file.split(".")[1]), 'w') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in NodeDict.items():
                writer.writerow([key, value])
        print("There are {} nodes".format(len(NodeDict)))
        return(NodeDict)

    def segmentCoords(self, df, part_axis, boundaries, max_thresh):

        # Segments the df into `num_threads` dfs
        s_dfs = []
        for i in range(self.num_threads):
            minv = boundaries[i]
            maxv = boundaries[i+1]
            s_dfs.append(df.query("{0}<=@maxv & {0}>=@minv".format(part_axis)))

        # Creates `num_threads`-1 boundary region dfs
        b_dfs = []
        for i in range(1,self.num_threads):
            bound_val = boundaries[i]
            bmin = bound_val - max_thresh
            bmax = bound_val + max_thresh
            b_dfs.append(df.query("{0}>=@bmin & {0}<=@bmax".format(part_axis))) # TODO: Add threshold

        return(s_dfs, b_dfs)

    def __init__(self, position_path, csn_type, min_thresh=0, max_thresh=8,
                 out_path=None, verbose = True):


        self.num_threads = 8

        assert csn_type in ["thresh"] # shadow networks are unimplemented (should implement maybe, for completeness)

        if verbose:
            print("Reading in {}".format(position_path.split("/")[-1]))

        # Reads in positional dataset
        df = pd.read_csv(position_path)

        # Creates node for each unique residue in positional dataset
        df["node_id"]  = df['chain'] + df["resi"].map(str)
        unique_nodes = list(set((df["node_id"])))

        NodeDict = self.makeNodeDict(unique_nodes, out_file)

        mins = {
         "x": np.amin(df['x']),
         "y": np.amin(df['y']),
         "z": np.amin(df['z'])
        }

        maxs = {
         "x": np.amax(df['x']),
         "y": np.amax(df['y']),
         "z": np.amax(df['z'])
        }

        axes = ["x","y","z"]
        ranges = [maxs[c]-mins[c] for c in axes]
        range_index = ranges.index(max(ranges))
        part_axis = axes[range_index]

        # Calculates boundaries
        bound_dists = np.multiply(range(self.num_threads+1), ranges[range_index]/(self.num_threads))
        boundaries = [mins[part_axis] + y for y in bound_dists]

        for axis in axes:
            print("{} goes from {} to {}".format(axis, mins[axis], maxs[axis]))
        print("Segmenting on {}".format(part_axis))
        print("Boundaries for {} given by {}".format(part_axis, boundaries))

        # Gets segmented dataframes to send to each thread, alongside
        # boundary condition dataframes that need to be dealt with at the end
        s_dfs, b_dfs = self.segmentCoords(df, part_axis, boundaries, max_thresh)

        # TODO: Remove these print statements
        print("---")
        [print(df.size) for df in s_dfs]
        print("---")
        [print(df.size) for df in b_dfs]

        # Calls ray for all segments
        x = [writeEdgeList.remote(sdf, NodeDict, min_thresh,
            max_thresh, out_file) for sdf in s_dfs]
        xres = ray.get(x)

        # Calls ray for all boundaries (requires x to finish)
        y = [writeEdgeList.remote(bdf, NodeDict, min_thresh,
            max_thresh, out_file, dependency=xres) for bdf in b_dfs]
        yres = ray.get(y)

        # Outputs final graph (requires y to finish)
        G = ray.get(createGraph.remote(out_file, NodeDict, dependency=yres))

        # Removes temporary files
        for file in glob.glob("data/CSN_graphs/*.tmp"):
            os.remove(file)
        raise




if __name__ == "__main__":

    ray.init()

    # Cleans any leftover temporary files from previous runs
    for file in glob.glob("data/CSN_graphs/*.tmp"):
        os.remove(file)

    # Should take the distance cutoffs as arguments to the Python file
    assert len(sys.argv) == 3
    min_thresh = float(sys.argv[1])
    max_thresh = float(sys.argv[2])

    print("-- Welcome to the Contact Network Constructor --")
    print("(Constructing basic network with min={}, max={:.2f})".format(
        min_thresh, max_thresh))
    print("\nWe hope your stay with us is pleasant and enjoyable")


    try:
        # Iterates over every position CSV within the folder
        for pos_file in glob.glob('data/positions/*.csv'):
            out_file = "data/CSN_graphs/{}_{:.2f}A_{}.gml".format(
                pos_file.split("/")[-1][:-4], max_thresh, "thresh")

            if not os.path.isfile(out_file):
                NetworkConstructor(pos_file, "thresh", 0, max_thresh, out_file, True)
            else:
                print("Already Exists: {}".format(out_file))
    except:
        # Cleans up temporary files so they're not persistent
        for file in glob.glob("data/CSN_graphs/*.tmp"):
            os.remove(file)
        raise
