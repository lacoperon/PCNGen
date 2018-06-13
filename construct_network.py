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
    __init__(self, position_path, csn_type, min_thresh=0,
             max_thresh=8, out_path=None):
             assert csn_type in ["thresh"] # shadow networks are unimplemented


    # Reads in positional dataset
    df = pd.read_csv(position_path)

    if csn_type is "thresh":
        print(df['resi'])
