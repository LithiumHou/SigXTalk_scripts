from torch.utils.data import DataLoader
import torch
import torch.nn.functional as F
import sys

work_dir = '/home/jiawen/myMLnet/benchmark/GENELink/Code'
if work_dir not in sys.path:
    sys.path.append(work_dir)

from torch.optim import Adam
from scGNN import GENELink
from torch.optim.lr_scheduler import StepLR
import scipy.sparse as sp
from utils import scRNADataset, load_data, adj2saprse_tensor, Evaluation,  Network_Statistic
import pandas as pd
from torch.utils.tensorboard import SummaryWriter
from PytorchTools import EarlyStopping
import numpy as np
import random
import csv

from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.model_selection import train_test_split

import itertools
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--lr', type=float, default=3e-3, help='Initial learning rate.')
parser.add_argument('--epochs', type=int, default = 100, help='Number of epoch.')
parser.add_argument('--num_head', type=list, default=[3,3], help='Number of head attentions.')
parser.add_argument('--alpha', type=float, default=0.2, help='Alpha for the leaky_relu.')
parser.add_argument('--hidden_dim', type=int, default=[128,64,32], help='The dimension of hidden layer')
parser.add_argument('--output_dim', type=int, default=16, help='The dimension of latent layer')
parser.add_argument('--batch_size', type=int, default=64, help='The size of each batch')
parser.add_argument('--loop', type=bool, default=False, help='whether to add self-loop in adjacent matrix')
parser.add_argument('--seed', type=int, default=128, help='Random seed')
parser.add_argument('--Type',type=str,default='dot', help='score metric')
parser.add_argument('--flag', type=bool, default=False, help='the identifier whether to conduct causal inference')
parser.add_argument('--reduction',type=str,default='concate', help='how to integrate multihead attention')
parser.add_argument('--round', type=int, default=1, help='Current round')


args = parser.parse_args()
args.flag = False
seed = args.seed
random.seed(args.seed)
torch.manual_seed(args.seed)
np.random.seed(args.seed)
device = torch.device("cuda:0")

########## 
## Construct the Rec-TF network

def Generate_negative_2(df):
    # Generate all possible combinations of X, Y values

    df2 = df.copy()
    df2.columns = ['From','To','label']
    N1 = set(df2.iloc[:,0])
    N2 = set(df2.iloc[:,1])

    all_combinations = list(itertools.product(N1, N2))
    df_all = pd.DataFrame(all_combinations, columns=['From','To'])

    # Filter out combinations that are already in the DataFrame
    merged_df = df_all.merge(df2, on=['From','To'], how='left', indicator=True)
    neg_samples = merged_df[merged_df['_merge'] == 'left_only']
    neg_samples = neg_samples.drop(columns='_merge')
    neg_samples['label'] = 0
    
    del all_combinations, merged_df
    return neg_samples

data_input = pd.read_csv('/home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_rectf/ExpressionData.csv',index_col=0) # row: genes, col: cells
loader = load_data(data_input)
feature_rectf = loader.exp_data()
feature_rectf = torch.from_numpy(feature_rectf)

RecTFDB = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_RecTF.csv',header=0, sep = ';')
rec = [int(item[4:]) if item.startswith('Gene') else item for item in list(set(RecTFDB['From']))]
rec = sorted(rec)
tf = [int(item[4:]) if item.startswith('Gene') else item for item in list(set(RecTFDB['To']))]
tf = sorted(tf)

allgenes = data_input.index.tolist()
string_to_index_rectf = {string: index for index, string in enumerate(allgenes)}
def map_strings_to_indexes_rectf(value):
    return string_to_index_rectf.get(value, value)

all_data = RecTFDB.copy()
all_data = all_data.applymap(map_strings_to_indexes_rectf)
all_data['label'] = 1

all_data.columns = ['From','To','label']
all_false = Generate_negative_2(all_data)
all_false = all_false.sample(n = len(all_data))
all_data = pd.concat([all_data,all_false],axis = 0)
training_data, val_data = train_test_split(all_data, test_size=0.2, stratify=all_data['label'])

data_feature_rectf = feature_rectf.to(device)
rec = torch.from_numpy(pd.DataFrame(rec).values.astype(np.int64))
rec = rec.to(device)

train_data = training_data.values
validation_data = val_data.values

train_load = scRNADataset(train_data, feature_rectf.shape[0], flag=args.flag)
adj_rectf = train_load.Adj_Generate(rec,loop=args.loop)
adj_rectf = adj2saprse_tensor(adj_rectf)
adj_rectf = adj_rectf.to(device)

train_data = torch.from_numpy(train_data).to(device)
validation_data = torch.from_numpy(validation_data).to(device)

model_rectf = GENELink(input_dim=feature_rectf.size()[1],
                hidden1_dim=args.hidden_dim[0],
                hidden2_dim=args.hidden_dim[1],
                hidden3_dim=args.hidden_dim[2],
                output_dim=args.output_dim,
                num_head1=args.num_head[0],
                num_head2=args.num_head[1],
                alpha=args.alpha,
                device=device,
                type=args.Type,
                reduction=args.reduction
                )
model_rectf = model_rectf.to(device)

optimizer = Adam(model_rectf.parameters(), lr=args.lr)
scheduler = StepLR(optimizer, step_size=1, gamma=0.99)
for epoch in range(args.epochs):
    running_loss = 0.0

    for train_x, train_y in DataLoader(train_load, batch_size=args.batch_size, shuffle=True):
        model_rectf.train()
        optimizer.zero_grad()

        if args.flag:
            train_y = train_y.to(device)
        else:
            train_y = train_y.to(device).view(-1, 1)


        # train_y = train_y.to(device).view(-1, 1)
        pred = model_rectf(data_feature_rectf, adj_rectf, train_x)

        #pred = torch.sigmoid(pred)
        if args.flag:
            pred = torch.softmax(pred, dim=1)
        else:
            pred = torch.sigmoid(pred)
        loss_BCE = F.binary_cross_entropy(pred, train_y)


        loss_BCE.backward()
        optimizer.step()
        scheduler.step()

        running_loss += loss_BCE.item()


    model_rectf.eval()
    score = model_rectf(data_feature_rectf, adj_rectf, validation_data)
    if args.flag:
        score = torch.softmax(score, dim=1)
    else:
        score = torch.sigmoid(score)
    
    val_y = validation_data[:, -1].to(device).view(-1, 1)
    loss_val = F.binary_cross_entropy(score, val_y.type(torch.float32))
    # score = torch.sigmoid(score)

    AUC, AUPR, AUPR_norm = Evaluation(y_pred=score, y_true=validation_data[:, -1],flag=args.flag)
        #
    # if epoch % 10 == 9:
    #     print('Epoch:{}'.format(epoch + 1),
    #         'train loss:{:.3F}'.format(running_loss),
    #         'val loss:{:.3F}'.format(loss_val),
    #         'AUC:{:.3F}'.format(AUC),
    #         'AUPR:{:.3F}'.format(AUPR))

#####
## Construct the TF-TG network

data_input = pd.read_csv('/home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_tftg/ExpressionData.csv',index_col=0) # row: genes, col: cells
loader = load_data(data_input)
feature_tftg = loader.exp_data()
feature_tftg = torch.from_numpy(feature_tftg)

TFTGDB = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_TFTG.csv',header=0, sep = ';')
tf = [int(item[4:]) if item.startswith('Gene') else item for item in list(set(TFTGDB['From']))]
tf = sorted(tf)
target = [int(item[4:]) if item.startswith('Gene') else item for item in list(set(TFTGDB['To']))]
target = sorted(target)

allgenes = data_input.index.tolist()
string_to_index_tftg = {string: index for index, string in enumerate(allgenes)}
def map_strings_to_indexes_tftg(value):
    return string_to_index_tftg.get(value, value)

all_data = TFTGDB.copy()
all_data = all_data.applymap(map_strings_to_indexes_tftg)
all_data['label'] = 1

all_data.columns = ['From','To','label']
all_false = Generate_negative_2(all_data)
all_false = all_false.sample(n = len(all_data))
all_data = pd.concat([all_data,all_false],axis = 0)
training_data, val_data = train_test_split(all_data, test_size=0.2, stratify=all_data['label'])

data_feature_tftg = feature_tftg.to(device)
tf = torch.from_numpy(pd.DataFrame(tf).values.astype(np.int64))
tf = tf.to(device)

train_data = training_data.values
validation_data = val_data.values

train_load = scRNADataset(train_data, feature_tftg.shape[0], flag=args.flag)
adj_tftg = train_load.Adj_Generate(rec,loop=args.loop)
adj_tftg = adj2saprse_tensor(adj_tftg)
adj_tftg = adj_tftg.to(device)

train_data = torch.from_numpy(train_data)
validation_data = torch.from_numpy(validation_data)
train_data = train_data.to(device)
validation_data = validation_data.to(device)


model_tftg = GENELink(input_dim=feature_tftg.size()[1],
                hidden1_dim=args.hidden_dim[0],
                hidden2_dim=args.hidden_dim[1],
                hidden3_dim=args.hidden_dim[2],
                output_dim=args.output_dim,
                num_head1=args.num_head[0],
                num_head2=args.num_head[1],
                alpha=args.alpha,
                device=device,
                type=args.Type,
                reduction=args.reduction
                )
model_tftg = model_tftg.to(device)

optimizer = Adam(model_tftg.parameters(), lr=args.lr)
scheduler = StepLR(optimizer, step_size=1, gamma=0.99)
for epoch in range(args.epochs):
    running_loss = 0.0

    for train_x, train_y in DataLoader(train_load, batch_size=args.batch_size, shuffle=True):
        model_tftg.train()
        optimizer.zero_grad()

        if args.flag:
            train_y = train_y.to(device)
        else:
            train_y = train_y.to(device).view(-1, 1)


        # train_y = train_y.to(device).view(-1, 1)
        pred = model_tftg(data_feature_tftg, adj_tftg, train_x)

        #pred = torch.sigmoid(pred)
        if args.flag:
            pred = torch.softmax(pred, dim=1)
        else:
            pred = torch.sigmoid(pred)

        loss_BCE = F.binary_cross_entropy(pred, train_y)
        loss_BCE.backward()
        optimizer.step()
        scheduler.step()

        running_loss += loss_BCE.item()


    model_tftg.eval()
    score = model_tftg(data_feature_tftg, adj_tftg, validation_data)
    if args.flag:
        score = torch.softmax(score, dim=1)
    else:
        score = torch.sigmoid(score)
    
    val_y = validation_data[:, -1].to(device).view(-1, 1)
    loss_val = F.binary_cross_entropy(score, val_y.type(torch.float32))
    # score = torch.sigmoid(score)

    AUC, AUPR, AUPR_norm = Evaluation(y_pred=score, y_true=validation_data[:, -1],flag=args.flag)
    
    # if epoch % 10 == 9:
    #     print('Epoch:{}'.format(epoch + 1),
    #         'train loss:{:.3F}'.format(running_loss),
    #         'val loss:{:.3F}'.format(loss_val),
    #         'AUC:{:.3F}'.format(AUC),
    #         'AUPR:{:.3F}'.format(AUPR))

    
# validate
model_rectf.eval()
model_tftg.eval()

val_set = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_val.csv',index_col = 0,header = 0)
val1 = val_set[['Receptor','TF','label']].values
val1 = torch.from_numpy(val1)
pred_rectf = model_rectf(data_feature_rectf,adj_rectf,val1)
if args.flag:
    pred_rectf = torch.softmax(pred_rectf, dim=1)
else:
    pred_rectf = torch.sigmoid(pred_rectf)
pred_rectf = pred_rectf.cpu().detach().numpy()

val2 = val_set[['TF','TG','label']].values
val2[:,0] = val2[:,0]-len(rec)
val2[:,1] = val2[:,1]-len(rec)

val2 = torch.from_numpy(val2)
pred_tftg = model_rectf(data_feature_tftg,adj_tftg,val2)
if args.flag:
    pred_tftg = torch.softmax(pred_tftg, dim=1)
else:
    pred_tftg = torch.sigmoid(pred_tftg)
pred_tftg = pred_tftg.cpu().detach().numpy()

preds = (pred_rectf)*(pred_tftg)
preds_new =  [0 if np.isnan(x) else x for x in preds]

labels = val_set['label'].astype(int)
AUC_final = roc_auc_score(y_true=labels, y_score=preds_new)
AUPR_final = average_precision_score(y_true=labels,y_score=preds_new)

print(' ------ The performance for GENELINK for round ' + str(args.round) + '------')
print('AUROC:{:.3F}'.format(AUC_final),
      'AUPR:{:.3F}'.format(AUPR_final))

# result = [args.round, 'GENELINK',AUC_final,AUPR_final]
# filename = '/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_comparison.csv'
# with open(filename, 'a', newline='') as file:
#     writer = csv.writer(file)
#     # Write the list as a new row in the CSV file
#     writer.writerow(result)
