import torch

import torch_geometric.transforms as T
from torch_geometric.datasets import Planetoid
from torch_geometric.nn import GAE, VGAE, GCNConv
from torch_geometric.data import Data
from torch_geometric.utils import negative_sampling
from torch_geometric.nn import GATv2Conv

import sys 
import pandas as pd
import os

def load_edgelist():
    edgelist_df = pd.read_csv("./artifacts/vgae/edgelist.csv")
    el = []    
    for index, row in edgelist_df.iterrows():
        el.append([row['bee']-1, row['plant']-1])
    return torch.tensor(el, dtype=torch.long).transpose(0,1)

def load_node_ids():
    return pd.read_csv("./artifacts/vgae/node_ids.csv")

def lookup_species_name(node_ids, id):
    return(node_ids[node_ids.node_id.values == (id+1)].species.values)

class GCNEncoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels)
        self.conv2 = GCNConv(2 * out_channels, out_channels)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv2(x, edge_index)


class VariationalGCNEncoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels)
        self.conv_mu = GCNConv(2 * out_channels, out_channels)
        self.conv_logstd = GCNConv(2 * out_channels, out_channels)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv_mu(x, edge_index), self.conv_logstd(x, edge_index)

def train(model, optimizer, train_data, variational=False):
    model.train()
    optimizer.zero_grad()
    z = model.encode(train_data.x, train_data.edge_index)
    loss = model.recon_loss(z, train_data.pos_edge_label_index)
    
    if variational:
        loss = loss + (1 / train_data.num_nodes) * model.kl_loss()
    loss.backward()
    optimizer.step()
    return float(loss)

@torch.no_grad()
def test(model, data):
    model.eval()
    z = model.encode(data.x, data.edge_index)
    return model.test(z, data.pos_edge_label_index, data.neg_edge_label_index)


def get_embedding_df(embedding, node_ids):
    df = pd.DataFrame()

    num_species, emb_dim = embedding.size()
    
    for i in range(emb_dim):
        df["dim_%d" % i] = [0.0 for i in range(num_species)]
    df["species"] = ["" for i in range(num_species)]
    
    for i,x in enumerate(embedding):
        sp = lookup_species_name(node_ids, i)
        df["species"][i] = sp[0]
        for j,y in enumerate(x):
            df["dim_%d" % j][i] = y.item()
        
        #for j in range(emb_dim):
            #=for k in range(num_species):
            #   df["dim_%d" % j][k] = embedding[i][j,k].item()
    return df


def fit_model(edgelist, node_ids, input_feature_dim, embed_dims):   
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    transform = T.Compose([
        T.NormalizeFeatures(),
        T.ToDevice(device),
        T.RandomLinkSplit(num_val=0., num_test=0.5, is_undirected=True,
                        split_labels=True, add_negative_train_samples=True,
                        neg_sampling_ratio=1.),
    ])
    
    noise_features = torch.randn_like(torch.zeros(NUM_SPECIES, input_feature_dim))
    data = Data(edge_index=edgelist, x=noise_features)
    train_data, val_data, test_data = transform(data)
    
    gae = GAE(GCNEncoder(input_feature_dim, embed_dims))
    vgae = VGAE(VariationalGCNEncoder(input_feature_dim, embed_dims))

    models = [gae, vgae]
    is_variational = [False, True]

    for i,model in enumerate(models):
        model = model.to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

        num_epochs = 50
        for epoch in range(1, num_epochs + 1):
            loss = train(model, optimizer, train_data, variational=is_variational[i])
            auc, ap = test(model, test_data)
            print(f'Model: {model}, Epoch: {epoch:03d}, AUC: {auc:.4f}, AP: {ap:.4f}')

    embeddings = [model.encode(noise_features, data.edge_index) for model in models]
    
    gae_df = get_embedding_df(embeddings[0], node_ids)
    vgae_df = get_embedding_df(embeddings[1], node_ids)
    
    os.makedirs("./artifacts/vgae/fits",exist_ok=True)
    
    vgae_out_path = "./artifacts/vgae/fits/VGAE_inputdim_%d_embeddim_%d.csv"  % (input_feature_dim, embed_dims)
    gae_out_path = "./artifacts/vgae/fits/GAE_inputdim_%d_embeddim_%d.csv" % (input_feature_dim, embed_dims)
    gae_df.to_csv(gae_out_path)
    vgae_df.to_csv(vgae_out_path)
    return  
    


edgelist = load_edgelist()
node_ids = load_node_ids()


NUM_SPECIES = 175

input_feature_dim = [1028, 256, 64, 16]
embed_dims = [64, 32, 16, 8, 4]


for inp in input_feature_dim:
    for emb in embed_dims:
        fit_model(edgelist, node_ids, inp, emb)



