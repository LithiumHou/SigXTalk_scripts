import csv
import numpy as np
import sys
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score

cur_round = sys.argv[1]
# cur_round = 1
# load the validation set
val_set = pd.read_csv("/home/jiawen/myMLnet/benchmark/datasets/sergio_val.csv", index_col=0)
val_set['Receptor'] = ['Gene'+str(i) for i in val_set['Receptor']]
val_set['TF'] = ['Gene'+str(i) for i in val_set['TF']]
val_set['TG'] = ['Gene'+str(i) for i in val_set['TG']]
RecTF = val_set[['Receptor','TF']]
TFTG = val_set[['TF','TG']]
# load all methods
all_method = ['PIDC', 'GRNVBEM', 'GRNBOOST2','PPCOR','SCODE','SINCERITIES','LEAP', 'GRISLI','SINGE','SCRIBE','SCSGL']
dff = pd.DataFrame(columns=['round','method','AUC','AUPR'])

for i in range(len(all_method)):
    method = all_method[i]

    file_rectf = '/home/jiawen/myMLnet/benchmark/Beeline/outputs/example/bench_rectf/'+method+'/rankedEdges.csv'
    res_rectf = pd.read_csv(file_rectf,sep="\t")
    res_rectf.columns = ['Receptor','TF','Weight']

    file_tftg = '/home/jiawen/myMLnet/benchmark/Beeline/outputs/example/bench_tftg/'+method+'/rankedEdges.csv'
    res_tftg = pd.read_csv(file_tftg,sep="\t")
    res_tftg.columns = ['TF','TG','Weight']

    if len(res_rectf)*len(res_tftg) == 0:
        result = [cur_round, method, 0.5, 0.5]
    else:
        res_rectf['Weight'] = res_rectf['Weight']/max(res_rectf['Weight'])
        res_tftg['Weight'] = res_tftg['Weight']/max(res_tftg['Weight'])
        filtered_rectf = pd.merge(res_rectf,RecTF,on=['Receptor','TF'],how='right')
        filtered_rectf.fillna(0,inplace=True)
        pred1 = filtered_rectf['Weight'].tolist()    
        filtered_tftg = pd.merge(res_tftg,TFTG,on=['TF','TG'],how='right')
        filtered_tftg.fillna(0,inplace=True)
        pred2 = filtered_tftg['Weight'].tolist()

        pred = [abs(x) * abs(y) for x, y in zip(pred1, pred2)]

        labels = val_set['label'].astype(int)
        AUC_final = roc_auc_score(y_true=labels, y_score=pred)
        AUPR_final = average_precision_score(y_true=labels,y_score=pred)
        result = [cur_round, method, AUC_final, AUPR_final]
        
    dff.loc[len(dff)] = result

filename = '/home/jiawen/myMLnet/benchmark/Beeline/outputs/comparison/comparison_'+str(cur_round)+'.csv'
dff.to_csv(filename, sep=',')
