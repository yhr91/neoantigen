import pandas as pd
import numpy as np
from glob import glob
import os

ID_dict = np.load('../tim_data/ID_dict.npy').item()
path = "./results/"
new_path = path + 'fixed'
files = glob(path+'*')

if not os.path.isdir(new_path):
    os.mkdir(new_path)
    
# Make the header change
def correction(line):
    
    period_splits = line.split('.')
    key = period_splits[0].split(' ')[-1]
    key = key.replace('_',':')
    new_key = ID_dict[key]

    new_line = ('Protein ' + key + ' ' + new_key 
                + '.' + '.'.join(period_splits[1:]))
    
    return new_line

# Read and correct file
def correct_file(f, new_path):
    with open(f, "r+") as read:
        text=read.readlines() 

    i=0 
    while i < len(text): 
        if text[i][0:7]=='Protein':
            text[i] = correction(text[i])
        i+=1
    
    return text[46:]

def add_file_header(text, filename_dict):
    with open(first, "r+") as read:
        complete_text=read.readlines()[:46]
    complete_text.extend(text)
    return complete_text

# For sorting filenames
filename_dict = {int(f.split('_')[-1]):f 
                 for f in files if f.split('/')[-1][0]=='H'}
first = filename_dict[0]

# Read and correct each file
full_text = []
for idx in sorted(filename_dict):
    print(idx)
    full_text.extend(correct_file(filename_dict[idx], new_path))
    
# Add file header
full_text = add_file_header(full_text, first)

# Write complete file
new_name = (new_path + '/'+ 
            first.split('/')[-1].split('_out_')[0])
with open(new_name,"w") as file: 
    file.writelines(full_text)