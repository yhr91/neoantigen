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

    new_name = new_path + '/'+ f.split('/')[-1]
    with open(new_name,"w") as file: 
        file.writelines(text)
        
for i,f in enumerate(files):
    print(i)
    correct_file(f, new_path)
