import pandas as pd
import torch

def load_edgelist():
    edgelist_df = pd.read_csv("data/edgelist.csv")
    el = []    
    for index, row in edgelist_df.iterrows():
        el.append([row['bee_node_id'], row['plant_node_id']])
    return torch.tensor(el, dtype=torch.long).transpose(0,1)
