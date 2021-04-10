#!/usr/bin/env python
# coding: utf-8

# In this notebook, we will be generating features to use for Deep-PIPE-Sites.

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
import json


# In[2]:


masks_path = "data/yeast_masks_singlesite_area_50_filtered.pkl"
uniprot_file = "data/uniprot-proteome UP000002311.tab"# reference proteome at https://www.uniprot.org/proteomes/UP000002311
pssm_path = "data/yeast_pssms/"

processed_data_path = "data/yeast_processed_norm_area_50/"
PIPE_dir = "data/PIPE_output/landscapes/yeast-yeast/"

# Find the max and min values across entire set and naively normoalize it with the assumption that the domain of the set is the domain of all possible data
NORM = True
# Some landscapes have extreme values at peaks in the PIPE landscape. This normalizes the PIPE landscape to the **average** peak value across all landscapes.
PIPE_AVE_MAX = False 
LOAD_PARAMS = False
load_params_path = processed_data_path + 'min_maxes.json'


# 

# In[ ]:





# # Load data and masks

# ## Load masks

# In[3]:


ppi_masks = pd.read_pickle(masks_path)
ppi_masks.head()


# In[4]:


ppi_masks.shape


# In[5]:


# All the protein sequences we need information for
proteins = ppi_masks['Uniprot ID A'].append(ppi_masks['Uniprot ID B'], ignore_index=True).unique()
proteins, proteins.shape


# ## Load sequences

# In[6]:


uniprot_df = pd.read_csv(uniprot_file,
                        sep = "\t", index_col='Entry')
print("Loaded UniProt proteome")
uniprot_df


# Create naive storage of min/max values for normalization

# In[7]:


if NORM:
    if LOAD_PARAMS:
        with open(load_params_path, 'r') as f:
            min_maxes = json.load(f)
    else:
        min_maxes = {}
    
def update_min_max(min_val, max_val, mask):
    if np.max(mask) > max_val:
        max_val = np.max(mask)
    if np.min(mask) < min_val:
        min_val = np.min(mask)
        
    return min_val, max_val

def norm(min_val, max_val, mask_list):
    return [(mat - min_val) / (max_val - min_val) for mat in mask_list]


# ## ProtDCal

# Generated using https://protdcal.zmb.uni-due.de/pages/form.php

# ## PSSMs

# psiblast -query 'uniprot-proteome UP000002311.fasta' -db ~/sysc4906/s/swissprot -inclusion_ethresh 0.001 -num_iterations 2 -out_ascii_pssm yeast_filtered.pssm -save_pssm_after_last_round
# 
# /home/wma/sysc4906/ncbi-blast-2.11.0+/bin
# 
# for file in *.fasta; do echo "./runPsiBlast.sh '$file'"; done >> blastall.sh
# 
# rename -n 's/sp\|(.*?)\|.*/$1.pssm/' *

# In[8]:


# psiblast -query .\output_sequences.fasta -db nr -out yeast_filtered_psiblast_out -evalue 0.001 -num_iterations 3 -out_pssm yeast_filtered_pssm_checkpoint -out_ascii_pssm yeast_filtered_pssm


# In[9]:


AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


# In[10]:


print("Loading PSSMs")
PSSMs = {}
for swissprot_id in proteins:
    pssm_A = pd.read_csv(pssm_path + swissprot_id + ".pssm", sep = "\s+", skiprows=2).reset_index()[AA].dropna()
    PSSMs[swissprot_id] = pssm_A.to_numpy()


# In[11]:


PSSM_masks = []

if NORM and not LOAD_PARAMS:
    max_val = np.NINF
    min_val = np.inf
print("Generating PSSM masks")
for UA, UB, SM in zip(ppi_masks['Uniprot ID A'], ppi_masks['Uniprot ID B'], ppi_masks['Sites Masks']):
    try:        
        # number according to position    
        seqA = PSSMs[UA]
        seqB = PSSMs[UB]
        
        
        # Do a dot product for each 20 length vector along each position
        mask = np.matmul(seqA, seqB.T)
        
        assert SM.shape == mask.shape
        
        if NORM and not LOAD_PARAMS:
            min_val, max_val = update_min_max(min_val, max_val, mask)
            
        PSSM_masks.append(mask)
    except KeyError as inst:
        print(UA, UB)
        print(f"No uniprot entry found for protein {inst.args}")
        position_landscape.append(np.NaN)


# In[12]:


# np.set_printoptions(threshold=1000)
PSSM_masks[:10]


# In[13]:


if NORM:
    if LOAD_PARAMS:
        PSSM_masks = norm(min_maxes['PSSM'][0], min_maxes['PSSM'][1], PSSM_masks)
    else:
        PSSM_masks = norm(min_val, max_val, PSSM_masks)
        min_maxes['PSSM'] = (min_val, max_val)
        
    print(PSSM_masks[:10])


# In[14]:


ppi_masks['PSSM_masks'] = PSSM_masks
ppi_masks


# In[15]:


# clear from memory
del PSSM_masks


# ## Position

# In[16]:


position_landscape = []
print("Generating position feature landscape")
for UA, UB, SM in zip(ppi_masks['Uniprot ID A'], ppi_masks['Uniprot ID B'], ppi_masks['Sites Masks']):
    try:
        lengths = (uniprot_df['Length'].loc[UA], uniprot_df['Length'].loc[UB])
        
        # number according to position    
        seqA = np.arange(1, lengths[0]+1)[np.newaxis]
        seqB = np.arange(1, lengths[1]+1)[np.newaxis]
        
        # normalize
        seqA = np.divide(seqA, lengths[0])
        seqB = np.divide(seqB, lengths[1])
        
        mask = np.matmul(seqA.T, seqB)
        
        assert SM.shape == mask.shape
        
        position_landscape.append(mask)
    except KeyError as inst:
        print(UA, UB)
        print(f"No uniprot entry found for protein {inst.args}")
        position_landscape.append(np.NaN)
        
position_landscape = np.asarray(position_landscape)
print(lengths, mask.shape)


# In[17]:


position_landscape[:2]


# In[18]:


# # The data are already normalized from 0 to 1
# if NORM:
#     min_maxes['Position'] = (0,1)


# In[19]:


ppi_masks['position_landscape'] = position_landscape
ppi_masks


# In[20]:


# free from memory
del position_landscape


# In[ ]:





# ## PIPE

# In[21]:


# Pandas is a little faster at file access and parsing than using Python

# import time

# PIPE_landscape = []
# PIPE_landscape_SW = []
# start = time.time()
# for protA, protB in zip(ppi_masks['Uniprot ID A'], ppi_masks['Uniprot ID B']):
#     PIPE_landscape_name = protA + '-' + protB
    
#     # load PIPE landscape
#     with open(PIPE_dir + PIPE_landscape_name + '.mat', 'r') as PIPE_file:
#         pipe_lines = PIPE_file.readlines()
    
#     # Last row and column is a sum, we don't need that
#     # get rid of \n too
#     landscape = np.zeros((len(pipe_lines)-1, len(pipe_lines[0].split(' '))-2), dtype=np.float)
#     for j in range(len(pipe_lines)-1):
#         line = np.array(pipe_lines[j].split(' ')[:-2]) # last two are sum, \n respectively
#         line = line.astype(np.float)
        
#         landscape[j] = line
        
#     PIPE_landscape.append(landscape)
    
#     # load PIPE landscape, sw score adjusted
#     with open(PIPE_dir + PIPE_landscape_name + '_SW.mat', 'r') as PIPE_file:
#         pipe_lines = PIPE_file.readlines()
    
#     # Last row and column is a sum, we don't need that
#     # get rid of \n too
#     landscape = np.zeros((len(pipe_lines)-1, len(pipe_lines[0].split(' '))-2), dtype=np.float)
#     for j in range(len(pipe_lines)-1):
#         line = np.array(pipe_lines[j].split(' ')[:-2]) # last two are sum, \n respectively
#         line = line.astype(np.float)
        
#         landscape[j] = line
        
#     PIPE_landscape_SW.append(landscape)
    
# end = time.time()
    
# PIPE_landscape[89].shape, PIPE_landscape_SW[89].shape, ppi_masks.PSSM_masks.iloc[89].shape, end-start       


# In[22]:


# import time
W = 20 # PIPE window size, PIPE matrix is smaller than protein sequence by W-1
pad = ((W-1) // 2, (W-1) - (W-1) //2)

if NORM and not LOAD_PARAMS:
    if PIPE_AVE_MAX:
        max_val = np.empty(len(ppi_masks))
        min_val = np.empty(len(ppi_masks))
    else:
        max_val = np.NINF
        min_val = np.inf
    max_val_sw = np.NINF
    min_val_sw = np.inf

PIPE_landscape = []
PIPE_landscape_SW = []
# start = time.time()
for i, (protA, protB) in enumerate(zip(ppi_masks['Uniprot ID A'], ppi_masks['Uniprot ID B'])):
    PIPE_landscape_name = protA + '-' + protB
    
    # load landscape
    landscape_df = pd.read_csv(PIPE_dir + PIPE_landscape_name + '.mat',
                               delim_whitespace=True,
                               header=None,
                               index_col=None)
    landscape = landscape_df.to_numpy().astype(np.float)
    landscape = np.pad(landscape[:-1, :-1], (pad, pad), 'minimum') # last row/column is a sum
    
    if NORM and not LOAD_PARAMS:
        if PIPE_AVE_MAX:
            max_val[i] = np.max(landscape)
            min_val[i] = np.min(landscape)
        else:
            min_val, max_val = update_min_max(min_val, max_val, landscape)
    
    PIPE_landscape.append(landscape) 
    
    # load sw landscape
    landscape_df = pd.read_csv(PIPE_dir + PIPE_landscape_name + '_SW.mat',
                               delim_whitespace=True,
                               header=None,
                               index_col=None)
    landscape = landscape_df.to_numpy().astype(np.float)
    landscape = np.pad(landscape[:-1, :-1], (pad, pad), 'minimum') # last row/column is a sum
    
    if NORM and not LOAD_PARAMS:
        min_val_sw, max_val_sw = update_min_max(min_val_sw, max_val_sw, landscape)
            
    PIPE_landscape_SW.append(landscape) # last row/column is a sum
    
if PIPE_AVE_MAX and NORM and not LOAD_PARAMS:
    max_val = max_val.mean()
    min_val = min_val.mean()

# end = time.time()
PIPE_landscape[89].shape, PIPE_landscape_SW[89].shape, ppi_masks.PSSM_masks.iloc[89].shape #, end-start


# In[23]:


if NORM:
    if LOAD_PARAMS:
        PIPE_landscape = norm(min_maxes['PIPE'][0],  min_maxes['PIPE'][1], PIPE_landscape)
        PIPE_landscape_SW = norm(min_maxes['PIPE_sw'][0], min_maxes['PIPE_sw'][1], PIPE_landscape_SW)
    else:
        PIPE_landscape = norm(min_val, max_val, PIPE_landscape)
        PIPE_landscape_SW = norm(min_val_sw, max_val_sw, PIPE_landscape_SW)

        min_maxes['PIPE'] = (min_val, max_val)    
        min_maxes['PIPE_sw'] = (min_val_sw, max_val_sw)
    
    print(PIPE_landscape[:2])
    print(PIPE_landscape_SW[:2])


# In[24]:


ppi_masks['PIPE_landscape'] = PIPE_landscape
ppi_masks['PIPE_landscape_SW'] = PIPE_landscape_SW
ppi_masks


# In[25]:


del PIPE_landscape
del PIPE_landscape_SW


# ## Prot2Vec

# In[ ]:





# In[ ]:





# In[26]:


ppi_masks.iloc[839]


# In[27]:


import matplotlib.pyplot as plt
from matplotlib import cm
query=839
fig, axs = plt.subplots(5, figsize=(15,10))
ax1 = axs[0].matshow(ppi_masks.PIPE_landscape.iloc[query])
axs[0].set_title('PIPE Landscape')
axs[0].set_ylabel('Position along Sequence A')
axs[0].set_xlabel('Position along Sequence B')
ax2 = axs[1].matshow(ppi_masks.PIPE_landscape_SW.iloc[query])
axs[1].set_title('PIPE Similarity-Weighted Adjusted Landscape')
axs[1].set_ylabel('Position along Sequence A')
axs[1].set_xlabel('Position along Sequence B')
ax3 = axs[2].matshow(ppi_masks.position_landscape.iloc[query])
axs[2].set_title('Position Landscape')
axs[2].set_ylabel('Position along Sequence A')
axs[2].set_xlabel('Position along Sequence B')
ax4 = axs[3].matshow(ppi_masks.PSSM_masks.iloc[query])
axs[3].set_title('PSSM Landscape')
axs[3].set_ylabel('Position along Sequence A')
axs[3].set_xlabel('Position along Sequence B')
ax5 = axs[4].matshow(ppi_masks['Sites Masks'].iloc[query])
axs[4].set_ylabel('Position along Sequence A')
axs[4].set_title('Interaction Site Mask')
axs[4].set_xlabel('Position along Sequence B')


plt.colorbar(ax1, ax=axs[0])
plt.colorbar(ax2, ax=axs[1])
plt.colorbar(ax3, ax=axs[2])
plt.colorbar(ax4, ax=axs[3])
plt.colorbar(ax5, ax=axs[4])

fig.tight_layout() # Or equivalently,  "plt.tight_layout()"

plt.show()


# # Save

# In[28]:


if not os.path.isdir(processed_data_path):
    os.mkdir(processed_data_path)

if NORM and not LOAD_PARAMS:
    print(min_maxes)
    
    with open(processed_data_path + 'min_maxes.json','w') as f:
        json.dump(min_maxes,f)


# In[31]:


ppi_masks.to_pickle(processed_data_path + 'data_processed_2021-01-25.pkl')


# ## Open and save

# In[32]:


ppi_masks = pd.read_pickle(processed_data_path + 'data_processed_2021-01-25.pkl')


# In[33]:


from sklearn.model_selection import train_test_split
train_df, holdout_df = train_test_split(ppi_masks, test_size=0.2, random_state=86)
test_df, val_df = train_test_split(holdout_df, test_size=0.5, random_state=789)


# In[34]:



# for UA, UB, site_mask, PSSM_mask, position_landscape in zip(ppi_masks['Uniprot ID A'], ppi_masks['Uniprot ID B'], ppi_masks['Sites Masks'], ppi_masks['PSSM_masks'], ppi_masks['position_landscape']):
#     features = np.stack([PSSM_mask, position_landscape])
#     np.save(f'{processed_data_path}masks/{UA}_{UB}.npy', site_mask)
#     np.save(f'{processed_data_path}features/{UA}_{UB}.npy', features)


# In[35]:


train_df


# In[36]:


test_df


# In[37]:


val_df


# In[38]:


if not os.path.isdir(processed_data_path + 'train/'):
    os.mkdir(processed_data_path + 'train/')
if not os.path.isdir(processed_data_path + 'test/'):
    os.mkdir(processed_data_path + 'test/')
if not os.path.isdir(processed_data_path + 'val/'):
    os.mkdir(processed_data_path + 'val/')

if not os.path.isdir(processed_data_path + 'train/masks/'):
    os.mkdir(processed_data_path + 'train/masks/')
if not os.path.isdir(processed_data_path + 'test/masks/'):
    os.mkdir(processed_data_path + 'test/masks/')
if not os.path.isdir(processed_data_path + 'val/masks/'):
    os.mkdir(processed_data_path + 'val/masks/')
if not os.path.isdir(processed_data_path + 'train/features/'):
    os.mkdir(processed_data_path + 'train/features/')
if not os.path.isdir(processed_data_path + 'test/features/'):
    os.mkdir(processed_data_path + 'test/features/')
if not os.path.isdir(processed_data_path + 'val/features/'):
    os.mkdir(processed_data_path + 'val/features/')


for UA, UB, site_mask, PSSM_mask, position_landscape, PIPE_landscape,  PIPE_landscape_SW in zip(test_df['Uniprot ID A'], 
                                                                                                test_df['Uniprot ID B'], 
                                                                                                test_df['Sites Masks'], 
                                                                                                test_df['PSSM_masks'], 
                                                                                                test_df['position_landscape'],
                                                                                                test_df['PIPE_landscape'],
                                                                                                test_df['PIPE_landscape_SW']):
    features = np.stack([PSSM_mask, position_landscape, PIPE_landscape,  PIPE_landscape_SW])
    np.save(f'{processed_data_path}test/masks/{UA}_{UB}.npy', np.expand_dims(site_mask, axis=0))
    np.save(f'{processed_data_path}test/features/{UA}_{UB}.npy', features)


# In[39]:


for UA, UB, site_mask, PSSM_mask, position_landscape, PIPE_landscape,  PIPE_landscape_SW in zip(val_df['Uniprot ID A'], 
                                                                                                val_df['Uniprot ID B'], 
                                                                                                val_df['Sites Masks'], 
                                                                                                val_df['PSSM_masks'], 
                                                                                                val_df['position_landscape'],
                                                                                                val_df['PIPE_landscape'],
                                                                                                val_df['PIPE_landscape_SW']):
    features = np.stack([PSSM_mask, position_landscape, PIPE_landscape,  PIPE_landscape_SW])
    np.save(f'{processed_data_path}val/masks/{UA}_{UB}.npy', np.expand_dims(site_mask, axis=0))
    np.save(f'{processed_data_path}val/features/{UA}_{UB}.npy', features)


# In[40]:




for UA, UB, site_mask, PSSM_mask, position_landscape, PIPE_landscape,  PIPE_landscape_SW in zip(train_df['Uniprot ID A'], 
                                                                                                train_df['Uniprot ID B'], 
                                                                                                train_df['Sites Masks'], 
                                                                                                train_df['PSSM_masks'], 
                                                                                                train_df['position_landscape'],
                                                                                                train_df['PIPE_landscape'],
                                                                                                train_df['PIPE_landscape_SW']):
    features = np.stack([PSSM_mask, position_landscape, PIPE_landscape,  PIPE_landscape_SW])
    np.save(f'{processed_data_path}train/masks/{UA}_{UB}.npy', np.expand_dims(site_mask, axis=0))
    np.save(f'{processed_data_path}train/features/{UA}_{UB}.npy', features)


# In[41]:


np.expand_dims(test_df['Sites Masks'].iloc[0], axis=0)


# In[ ]:





# In[42]:


train_df.to_pickle(processed_data_path + 'train.pkl')
test_df.to_pickle(processed_data_path + 'test.pkl')
val_df.to_pickle(processed_data_path + 'val.pkl')


# ## Create json of just bounding box labels

# In[3]:


train_df = pd.read_pickle(processed_data_path + 'train.pkl')
test_df = pd.read_pickle(processed_data_path + 'test.pkl')
val_df = pd.read_pickle(processed_data_path + 'val.pkl')

train_coords = {}
test_coords = {}
val_coords = {}
for UA, UB, DPA, DPB in zip(test_df['Uniprot ID A'], 
                            test_df['Uniprot ID B'], 
                            test_df['Domain positions A'], 
                            test_df['Domain positions B']):
    test_coords[f'{UA}_{UB}'] = DPA + DPB
    
for UA, UB, DPA, DPB in zip(train_df['Uniprot ID A'], 
                            train_df['Uniprot ID B'], 
                            train_df['Domain positions A'], 
                            train_df['Domain positions B']):
    train_coords[f'{UA}_{UB}'] = DPA + DPB
    
for UA, UB, DPA, DPB in zip(val_df['Uniprot ID A'], 
                            val_df['Uniprot ID B'], 
                            val_df['Domain positions A'], 
                            val_df['Domain positions B']):
    val_coords[f'{UA}_{UB}'] = DPA + DPB
test_coords


# In[4]:


# json doesn't like numpy
def convert(o):
    if isinstance(o, np.int32): return int(o)  
    raise TypeError
    
with open(processed_data_path + 'test/coords.json','w') as f:
    json.dump(test_coords,f, default=convert)
with open(processed_data_path + 'train/coords.json','w') as f:
    json.dump(train_coords,f, default=convert)
with open(processed_data_path + 'val/coords.json','w') as f:
    json.dump(val_coords,f, default=convert)


# In[5]:


type(test_coords['P20459_P14741'][0][0])


# ## TESTING

# In[6]:


train_df = pd.read_pickle(processed_data_path + 'train.pkl')


# In[7]:


train_df


# In[8]:


import torch.nn.functional as F
import torch


def ave_col(series):
    means = np.empty((len(series)))
    for i, landscape in enumerate(series):
        means[i] = landscape.mean()
        
    return means.mean().astype(np.float32)

def std_col(series, mean):
    squared_error_sums = np.empty((len(series)))
    for i, landscape in enumerate(series):
        std_landscape = torch.from_numpy(landscape).float()
        # print(std_landscape.shape)
        std_landscape = F.interpolate(std_landscape.unsqueeze(0).unsqueeze(0), (4910, 4910))
        std_landscape = std_landscape.detach().numpy().flatten()
        std_landscape = (np.abs(std_landscape - mean))**2
                                  
        squared_error_sums[i] = std_landscape.sum()
                                  
    variance = squared_error_sums.sum() / (len(series)*4910*4910 - 1)
        
    return np.sqrt(variance).astype(np.float32)


# In[9]:


avs = [ave_col(train_df.PSSM_masks), ave_col(train_df.position_landscape), ave_col(train_df.PIPE_landscape), ave_col(train_df.PIPE_landscape_SW)]


# In[10]:


stds = [std_col(train_df.PSSM_masks, avs[0]), std_col(train_df.position_landscape, avs[1]), std_col(train_df.PIPE_landscape, avs[2]), std_col(train_df.PIPE_landscape_SW, avs[3])]


# In[11]:


avs, stds


# In[12]:


def convert(o):
    if isinstance(o, np.float32): return float(o)  
    raise TypeError

with open(processed_data_path + 'stds.json','w') as f:
    json.dump([avs, stds], f, default=convert)


# In[ ]:




