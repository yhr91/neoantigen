import pandas as pd
import glob
import numpy as np

# Read in data
xls_files = glob.glob('../scripts/results/*xls')
xls_files = {int(name.split('_')[-1].split('.xls')[0]):name for name in xls_files}
ID_dict = np.load('../tim_data/ID_dict.npy').item()

# Helper functions
# Given XID get corresponding header
def get_header(key):
    key = key.replace('_',':')
    header = ID_dict[key]
    header_items = header.split('|')
    return [h.split(':')[-1] for h in header_items]

# Given a list of ID's get all corresponding headers in a dataframe format
def get_header_df(headers):
    header_df = pd.DataFrame(np.array([get_header(v) for v in headers]))
    header_df = header_df.rename(columns = {0:'gene_id',
                                1:'transcript_id',
                                2:'gene_name',
                                3:'type',
                                4:'sgRNA',
                                5:'strand'})
    return header_df

def add_header_df(df):
    header_df = get_header_df(df['ID'].values)
    base_df = df.reset_index()
    colnames = np.concatenate([base_df.columns,header_df.columns])
    concat_df = pd.concat([base_df,header_df], axis=1, ignore_index=True)
    concat_df = concat_df.rename(columns = 
                        {i:colnames[i] for i in range(len(colnames))})
    return concat_df

# Main script
all_dfs = []
for i,f in xls_files.items():
    print(i)
    df = pd.read_csv(f, delimiter = '\t', index_col=0, header=1)
    all_dfs.append(add_header_df(df))
all_dfs_concat = pd.concat(all_dfs, ignore_index=True)
all_dfs_concat.to_csv('combined.csv')
