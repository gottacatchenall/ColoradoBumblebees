import torch

import torch_geometric.transforms as T
from torch_geometric.datasets import Planetoid
from torch_geometric.nn import GAE, VGAE, GCNConv
from torch_geometric.data import Data
from torch_geometric.utils import negative_sampling
from torch_geometric.nn import GATv2Conv

import sys 
sys.path.append('./scripts/02_analysis/02_engineer_features/01_vgae')
from load_metaweb import load_edgelist

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
transform = T.Compose([
    T.NormalizeFeatures(),
    T.ToDevice(device),
    T.RandomLinkSplit(num_val=0.05, num_test=0.1, is_undirected=True,
                      split_labels=True, add_negative_train_samples=True,
                      neg_sampling_ratio=4.0),
])

edgelist = load_edgelist()

noise_features = torch.randn_like(feature_mat)

data = Data(edge_index=edgelist, x=noise_features)
train_data, val_data, test_data = transform(data)



#dataset = Planetoid("", "Cora", transform=transform)
#cora, _, _ = dataset[0]


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


class LinearEncoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv = GCNConv(in_channels, out_channels)

    def forward(self, x, edge_index):
        return self.conv(x, edge_index)


class VariationalLinearEncoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv_mu = GCNConv(in_channels, out_channels)
        self.conv_logstd = GCNConv(in_channels, out_channels)

    def forward(self, x, edge_index):
        return self.conv_mu(x, edge_index), self.conv_logstd(x, edge_index)

""" 
if not args.variational and not args.linear:
    model = GAE(GCNEncoder(in_channels, out_channels))
elif not args.variational and args.linear:
    model = GAE(LinearEncoder(in_channels, out_channels))
elif args.variational and not args.linear:
    model = VGAE(VariationalGCNEncoder(in_channels, out_channels))
elif args.variational and args.linear:
    model = VGAE(VariationalLinearEncoder(in_channels, out_channels))
"""

def train(variational=False):
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
def test(data):
    model.eval()
    z = model.encode(data.x, data.edge_index)
    return model.test(z, data.pos_edge_label_index, data.neg_edge_label_index)


in_channels, out_channels = train_data.num_features, 8


models = [
    GAE(LinearEncoder(in_channels, out_channels)),
    VGAE(VariationalLinearEncoder(in_channels, out_channels)),
    GAE(GCNEncoder(in_channels, out_channels)),
    VGAE(VariationalGCNEncoder(in_channels, out_channels)),
]

is_variational = [False, True, False, True]


for i,model in enumerate(models):
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

    num_epochs = 1500
    for epoch in range(1, num_epochs + 1):
        loss = train(variational=is_variational[i])
        auc, ap = test(test_data)
        print(f'Model: {model}, Epoch: {epoch:03d}, AUC: {auc:.4f}, AP: {ap:.4f}')



# this is the embedding for each node
embeddings = [model.encode(noise_features, data.edge_index) for model in models]

import pandas as pd
feature_df = pd.read_csv("data/features.csv")

def lookup_species_name(id):
    return feature_df[feature_df["node_id"] == id].species.values[0]

model_names = ["Linear_GAE", "Linear_VGAE", "GAE", "VGAE"]

dict = {}

for model_id,e in enumerate(embeddings):
    this_dict = {}
    for i,x in enumerate(e):
        this_dict[lookup_species_name(i)] = str(list(x.detach().numpy()))
    dict[model_names[model_id]] = this_dict



import json
json_object = json.dumps(dict, indent=4)

json_object

with open("emb.json", "w") as outfile:
    outfile.write(json_object)