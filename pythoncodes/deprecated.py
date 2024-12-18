import numpy as np
import pandas as pd
import xgboost as xgb
import sys
from tqdm import tqdm
work_dir = '/home/jiawen/myMLnet/pythoncodes'
sys.path.append(work_dir)

from preprocessing import *
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance

def XGBmodel_single(cur_row,exp):

    cur_TF = cur_row.iloc[0,1]
    cur_TG = cur_row.iloc[0,2]
    Recs_target = cur_row['Receptor'].tolist()
    exp_TF = np.array(exp[cur_TF])
    exp_TG = np.array(exp[cur_TG])

    exp_TFTG = exp_TF*exp_TG
    # exp_TFTG = np.log2(1+exp_TFTG)
    exp_recs = exp[Recs_target]
    # exp_recs = np.log2(1+exp_recs)

    xgb_model = xgb.XGBRegressor(n_estimators=1000, learning_rate=0.1, max_depth=5, random_state=2024,device = 'cuda')
    xgb_model.fit(exp_recs,exp_TFTG)

    importance = xgb_model.get_booster().get_score(importance_type='gain')
    df_imp = pd.DataFrame(list(importance.items()), columns=['Receptor', 'importance'])
    df_imp['TF'] = cur_TF
    df_imp['TG'] = cur_TG
    df_imp = df_imp[['Receptor','TF','TG','importance']]

    return df_imp

def RFmodel_single(cur_row,exp):

    cur_TF = cur_row.iloc[0,1]
    cur_TG = cur_row.iloc[0,2]
    Recs_target = cur_row['Receptor'].tolist()
    exp_TF = np.array(exp[cur_TF])
    exp_TG = np.array(exp[cur_TG])

    exp_TFTG = exp_TF*exp_TG
    # exp_TFTG = np.log2(1+exp_TFTG)
    exp_recs = exp[Recs_target]
    # exp_recs = np.log2(1+exp_recs)

    rf_model = RandomForestRegressor(n_estimators=200, learning_rate=0.1, max_depth=5, random_state=2024,device = 'cuda')
    rf_model.fit(exp_recs,exp_TFTG)

    importance = permutation_importance(rf_model, exp_recs,exp_TFTG, n_repeats=20,random_state=0)

    df_imp = pd.DataFrame(data = {'Receptor': Recs_target, 'importance': importance.importances_mean})
    df_imp['TF'] = cur_TF
    df_imp['TG'] = cur_TG
    df_imp = df_imp[['Receptor','TF','TG','importance']]

    return df_imp

def PCR_single(cur_row,exp,n_comp = 5):

    cur_TF = cur_row.iloc[0,1]
    cur_TG = cur_row.iloc[0,2]
    Recs_target = cur_row['Receptor'].tolist()
    exp_TF = np.array(exp[cur_TF])
    exp_TG = np.array(exp[cur_TG])

    exp_TFTG = exp_TF*exp_TG
    exp_TFTG = np.log2(1+exp_TFTG)
    exp_recs = exp[Recs_target].values
    exp_recs = np.log2(1+exp_recs)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(exp_recs)
    # Perform PCA
    pca = PCA()
    X_train_pca = pca.fit_transform(X_scaled)
    X_train_pca_reduced = X_train_pca[:, :n_comp]
    model = LinearRegression()
    model.fit(X_train_pca_reduced, exp_TFTG)

    # Compute the original variable coefficients
    pc_coefficients = model.coef_
    loadings = pca.components_[:n_comp].T
    importance = loadings @ pc_coefficients

    return tuple(importance)

def Regressor(df, exp, method = "xgboost"):
    # method: xgboost, RF, pcregression

    # Ensure the DataFrame has a unique index
    df = df.reset_index(drop=True)
    exp_nor = NormalizeData(exp)

    # Grouping DataFrame by 'y' and 'z'
    grouped = df.groupby(['TF', 'TG'])
    
    print("Start regression ...")
    
    # with mp.Pool(24) as pool:
    #     results = pool.starmap(RFmodel_single, tasks)
    count = 0
    res_list = [None for _ in range(len(grouped))]
    with tqdm(total=len(grouped)) as pbar:
        pbar.set_description('Processing:')
        for _,group in grouped:
            imp = XGBmodel_single(group,exp_nor)
            res_list[count] = imp
            count+=1
            pbar.update(1)

    # Flatten the list of results and sort by the original index
    flat_results = pd.concat(res_list)

    return flat_results

# filter the databases
def filter_DB(RecTFDB, TFTGDB, gene_all, recs = None):

    # RecTFDB = RecTFDB[RecTFDB.apply(lambda row: all(value in gene_all for value in row), axis=1)]
    # TFTGDB = TFTGDB[TFTGDB.apply(lambda row: all(value in gene_all for value in row), axis=1)]
    RecTFDB = RecTFDB.iloc[:,range(2)]
    TFTGDB = TFTGDB.iloc[:,range(2)]
    RecTFDB.columns = ['From','To']
    TFTGDB.columns = ['From','To']

    if recs is None:
        recs = list(set(RecTFDB['From']) - set(RecTFDB['To']))

    recs = set(recs).intersection(set(gene_all))
    recs = recs.intersection(set(RecTFDB['From']))
    recs = list(recs)
    tfs_1 = list(set(TFTGDB['From']))
    tfs_2 = list(set(RecTFDB['To']))
    tfs_ori = set(tfs_1).intersection(set(tfs_2))
    tfs_ori = list(set(tfs_ori).intersection(set(gene_all)))
    tgs_ori = list(set(TFTGDB['To']).intersection(set(gene_all)))

    tgs_new = list(set(tgs_ori) - set(recs+tfs_1))
    tfs_new = list(set(tfs_ori) - set(tgs_new+recs))
    recs_new = recs

    RecTFDB = RecTFDB[RecTFDB['From'].isin(recs_new)]
    RecTFDB = RecTFDB[RecTFDB['To'].isin(tfs_new)]
    TFTGDB = TFTGDB[TFTGDB['From'].isin(tfs_new)]
    TFTGDB = TFTGDB[TFTGDB['To'].isin(tgs_new)]

    return RecTFDB,TFTGDB

def construct_graph(cor_mat,g_thres):
    
    cor_mat = cor_mat - np.eye(len(cor_mat))

    top_percent_threshold = np.percentile(cor_mat.values, (1-g_thres)*100)
    filtered_mat = (cor_mat > top_percent_threshold).astype(int)
    filtered_mat = filtered_mat.values
    # Convert the adj matrix to edge lists
    rows, cols = np.triu_indices_from(filtered_mat, k=1)
    edges_indices = filtered_mat[rows, cols] > 0
    edges = list(zip(rows[edges_indices], cols[edges_indices]))
    
    genes = cor_mat.columns
    G = dhg.Graph(len(genes),edges)
    
    return G

# construct the hypergraph using the correlation/mi matrix
def construct_hypergraph(df1, df2, hg_thres):

    df1.columns = ['Receptor','TF','Weight1']
    df2.columns = ['TF','TG','Weight2']

    df_all = pd.merge(df1,df2,on='TF')
    df_all['Weight'] = df_all['Weight1']*df_all['Weight2']
    df_all['Weight'] = df_all['Weight'].abs()
    del df_all['Weight1']
    del df_all['Weight2']

    df_sorted = df_all.sort_values(by='Weight',ascending=False)
    df_sorted.reindex(columns=['Receptor','TF','TG','Weight'])

    all_len = len(df_sorted)
    hg_index = int(hg_thres*all_len)
    df_hg = df_sorted.iloc[:hg_index]

    df_hg_filtered = df_hg[['Receptor','TF','TG']]

    return df_hg_filtered, df_sorted

# generate positive training samples for the model
def Generate_positive(df,sample_thres):

    df_sorted = df.sort_values(by='Weight',ascending=False)
    df_sorted.reindex(columns=['Receptor','TF','TG','Weight'])

    all_len = len(df_sorted)
    sample_index = int(sample_thres*all_len)    
    df_samples = df_sorted.iloc[:sample_index]
    df_samples['label'] = 1
    del df_samples['Weight']

    return df_samples

def Generate_negative_2(df):
    # Generate all possible combinations of X, Y values

    df2 = df.copy()
    df2.columns = ['From','To','label']
    N1 = set(df2.iloc[:,0])
    N2 = set(df2.iloc[:,1])

    all_combinations = list(itertools.product(N1, N2))
    df_all = pd.DataFrame(all_combinations, columns=['From','To'])

    # Filter out combinations that are not already in the DataFram
    merged_df = df_all.merge(df2, on=['From','To'], how='left', indicator=True)
    neg_samples = merged_df[merged_df['_merge'] == 'left_only']
    neg_samples = neg_samples.drop(columns='_merge')
    neg_samples['label'] = 0
    
    del all_combinations, merged_df
    return neg_samples