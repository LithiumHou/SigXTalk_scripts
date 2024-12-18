import torch
import torch.nn.functional as F
import random
import numpy as np
import pandas as pd
import csv
import argparse
import dhg
from sklearn.model_selection import train_test_split
import sys

work_dir = '/home/jiawen/myMLnet/pythoncodes'
if work_dir not in sys.path:
    sys.path.append(work_dir)

from preprocessing import *
from predictor import *
from training import *

parser = argparse.ArgumentParser(description="Hyperlink prediction.")

parser.add_argument('--thres', type=float, default=[0.4], help='The thresholds for hypergraph construction')
parser.add_argument('--corr_type', type=str, default='Spearman', help='The method to calculate correlations')
parser.add_argument('--sample_size', type=int, default=600, help='Total number of training+validation samples.')
parser.add_argument('--train_size', type=float, default=0.75, help='Fraction of training samples.')
parser.add_argument('--hgnn_dims', type=int, default=[128,64], help='The dimension of hidden layer')
parser.add_argument('--linear_dims', type=int, default=[32, 16], help='The dimensions of the MLP layer')
parser.add_argument('--lr', type=float, default= 0.01, help='Initial learning rate.')
parser.add_argument('--epochs', type=int, default = 200, help='Number of epoch.')
parser.add_argument('--batch_size', type=int, default=64, help='The size of each batch')
parser.add_argument('--seed', type=int, default=2024, help='Random seed')
parser.add_argument('--epoch_output', type=bool, default=False, help='Whether output metrics for each epoch')
parser.add_argument('--metrics_output', type=bool, default=False, help='Whether output the auroc and aupr curves')
parser.add_argument('--round', type=int, default=1, help='Current round')


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")
    
args = parser.parse_args()
args.seed = 10*int(args.round)

random.seed(args.seed)
torch.manual_seed(args.seed)
np.random.seed(args.seed)

# Load the data and prior knowledge
input_exp = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_data.csv', sep = "\t") 
input_exp = input_exp.astype(np.float32) # row: genes, col: cells
RecTFDB = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_RecTF.csv',header=0, sep = ';')
TFTGDB = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_TFTG.csv',header=0, sep = ';')
RecTFDB.columns = ['Receptor','TF']
TFTGDB.columns = ['TF','TG']
RecTFTG = pd.merge(RecTFDB,TFTGDB,on='TF') # all the pathways

# Filter the input expression matrix
gene_all = input_exp.index.tolist()
input_exp = input_exp.loc[gene_all]
input_exp_nor = data_scale(input_exp,pca=False,n_comp=250)

input_used = input_exp_nor.copy()
input_tensor = torch.tensor(input_used.values, dtype=torch.float32)
input_tensor = input_tensor.to(device)


string_to_index = {string: index for index, string in enumerate(gene_all)}
def map_strings_to_indexes(value):
    return string_to_index.get(value, value)
    

# Generate the negative samples
positive_samples = RecTFTG.applymap(map_strings_to_indexes)
positive_samples['label'] = 1
negative_samples = Generate_negative_3(positive_samples,args)
negative_samples = negative_samples.sample(n = len(positive_samples))
negative_samples['label'] = 0


all_samples = pd.concat([positive_samples,negative_samples],axis = 0)
all_samples[['Receptor','TF','TG']] = all_samples[['Receptor','TF','TG']].applymap(map_strings_to_indexes)
training_data, val_data = train_test_split(all_samples, test_size=1-args.train_size-0.001, stratify=all_samples['label'],random_state=0)
val_data.to_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_val.csv')

# Construct the hypergraph
pos_training = training_data[training_data['label'] == 1]
pos_rectf = list(set([tuple(x) for x in pos_training[['Receptor','TF']].values]))
pos_tftg = list(set([tuple(x) for x in pos_training[['TF','TG']].values]))
pos_rectf = pd.DataFrame(pos_rectf, columns=['Receptor','TF'])
pos_tftg = pd.DataFrame(pos_tftg, columns=['TF','TG'])

num_genes = len(gene_all)

edge_list = [tuple(row) for row in pos_training[['Receptor','TF','TG']].values]
hg = dhg.Hypergraph(num_genes, edge_list)
hg = hg.to(device)

# Construct the model
mymodel = HGNNPredictor(
    in_channels = input_tensor.size()[1],
    hgnn_channels = args.hgnn_dims,
    linear_channels = args.linear_dims
).to(device)
    
mymodel, auroc, aupr = train(args, input_tensor, mymodel, hg, training_data, val_data, device)

# auroc = metrics['AUROC'].iloc[-1]
# aupr = metrics['AUPR'].iloc[-1]

# print(' ------ The performance for SigXTalk ------')
# print(' - AUROC:{:.3F}'.format(auroc),', AUPR: {:.3F}'.format(aupr))

### Predict the results
mymodel.eval()
pred = Predict(Exp = input_used, model = mymodel, input_samples = val_data,
                       hypergraph = hg, genes = gene_all, device = device)
pred['label'] = val_data['label'].astype(int)
AUC_final = roc_auc_score(y_true=pred['label'], y_score=pred['pred_label'])
AUPR_final = average_precision_score(y_true=pred['label'], y_score=pred['pred_label'])

print(' ------ The performance for SigXTalk for round ' + str(args.round) + '------')
print('AUROC:{:.3F}'.format(AUC_final),
      'AUPR:{:.3F}'.format(AUPR_final))


result = [args.round, args.seed,AUC_final,AUPR_final]
filename = '/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_randomseed.csv'
with open(filename, 'a', newline='') as file:
    writer = csv.writer(file)
    # Write the list as a new row in the CSV file
    writer.writerow(result)

# dfs.to_csv('/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_curve.csv')