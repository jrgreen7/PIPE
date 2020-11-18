#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re


# In[2]:


pfam_file = 'data/pfam_yeast_domains.tsv' # http://pfam.xfam.org/proteome/559292#tabview=tab2
biogrid_file = "data/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.2.191.tab3.txt"
threedid_file = "data/3did_flat_Apr_10_2020.dat" # from https://3did.irbbarcelona.org/download.php
uniprot_file = "data/uniprot-proteome UP000002311.tab"# reference proteome at https://www.uniprot.org/proteomes/UP000002311
# general proteome (reviewed S. cerevisiae) at https://www.uniprot.org/uniprot/?query=taxonomy:%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20/%20S288c)%20(Baker%27s%20yeast)%20[559292]%22&fil=organism%3A%22Saccharomyces+cerevisiae+%28strain+ATCC+204508+%2F+S288c%29+%28Baker%27s+yeast%29+%5B559292%5D%22+AND+reviewed%3Ayes
split_homologous = True
pfam_domains_only = True
single_domain_only = True
max_interaction_area = 0.5

save_file_path = "data/masks_singlesite_domains_only.pkl"


# # Import PFAM data

# In[3]:


# load in header row as single string # names are between angle brackets
pfam_header = re.findall(r'<(.*?)>', pd.read_csv(pfam_file, sep='\n', header=2).columns[0]);
pfam = pd.read_csv(pfam_file, sep='\t', names=pfam_header);
pfam = pfam.drop(range(3)).reset_index(drop=True)
print(f"Pfam data loaded: Shape: {pfam.shape}, Unique proteins: {pfam['seq id'].unique().size}")


# We are choosing to drop all types that are not explicitly domain

# In[4]:


print(f"Pfam subsequence types: {pfam['type'].unique()}")
if pfam_domains_only:
    print('Dropping pfam entries that are not domains')
    pfam = pfam[pfam['type'] == 'Domain'].sort_values('seq id')
    pfam.reset_index(drop=True, inplace=True)
    print(f"Pfam data after dropping: Shape: {pfam.shape}, Unique proteins: {pfam['seq id'].unique().size}")
else:
    print("Keep all types")


# In[5]:


# convert indices to int
pfam['envelope start'] = pfam['envelope start'].round(0).astype(np.int)
pfam['envelope end'] = pfam['envelope end'].round(0).astype(np.int)


# ### Create map of sequence IDs and their indices

# In[6]:


print("Creating map of sequence IDs and their indices in PFAM")
pfam_indices = dict.fromkeys(pfam['seq id'].unique())


# Find all locations in PFAM dataframe where index occurs and add it to map

# In[7]:


for key in pfam_indices:
    pfam_indices[key] = pfam['seq id'][pfam['seq id'] == key].index


# In[8]:


# pfam.loc[list(pfam_indices.values())[0]]


# ## Load BioGRID Interactome for yeast

# In[9]:


biogrid = pd.read_csv(biogrid_file,
                        sep = "\t")


# In[10]:


print(f"BioGRID data loaded: Shape: {biogrid.shape}")


# In[11]:


if split_homologous:
    print("Splitting SWISS-PROT sequences with  \"|\"")
    biogrid['SWISS-PROT Accessions Interactor A'] = biogrid['SWISS-PROT Accessions Interactor A'].apply(lambda s: s.split('|'))
    biogrid['SWISS-PROT Accessions Interactor B'] = biogrid['SWISS-PROT Accessions Interactor B'].apply(lambda s: s.split('|'))
    biogrid = biogrid.explode('SWISS-PROT Accessions Interactor A', ignore_index=True).explode('SWISS-PROT Accessions Interactor B', ignore_index=True)
    print(f"BioGRID data after split: Shape: {biogrid.shape}")


# Drop all for where swiss-prot ID is not available

# In[12]:


biogrid.drop(biogrid[biogrid['SWISS-PROT Accessions Interactor A'] == '-'].index, inplace=True)
biogrid.drop(biogrid[biogrid['SWISS-PROT Accessions Interactor B'] == '-'].index, inplace=True)
print(f"BioGRID data after dropping \"-\": Shape: {biogrid.shape}")


# ## Load domain-domain interactinos

# In[13]:



threeDID = pd.read_csv(threedid_file,
                    sep = "\t", header = None,
                    names=range(7))


# In[14]:


# this is a flat database; find indices with #=ID, that's what we're interested in (domain domain pairs)
id_indices_3did = threeDID[0].loc[threeDID[0] == '#=ID'].index
# strip all that are not those indices in non-pfam columns
threeDID = threeDID[[3,4]].iloc[id_indices_3did]


# In[16]:


# Reformat so it's just pfam IDs
threeDID.columns=["Pfam ID A","Pfam ID B"]
threeDID.reset_index(inplace=True, drop=True)

# Strip extra characters
threeDID['Pfam ID A'] = threeDID['Pfam ID A'].apply(lambda x: re.findall(r'\((.+?)\.', x)[0])
threeDID['Pfam ID B'] = threeDID['Pfam ID B'].apply(lambda x: re.findall(r'.*(?=\.)', x)[0])

print(f"3did Pfam IDs loaded: {len(threeDID)} Unique column A: {threeDID['Pfam ID A'].unique().size} Unique Column B: {threeDID['Pfam ID B'].unique().size}")


# ### Create dictionary of domain pairs

# In[ ]:





# In[17]:


# create a dict using unique PFAM IDs in 3did, initialized with empty lists
print("Creating dict of PFAM pairs")
pfam_pairs = {k : [] for k in pd.concat([threeDID['Pfam ID A'], threeDID['Pfam ID B']]).unique()}

# for all pairs add to dict
for A, B in zip(threeDID['Pfam ID A'], threeDID['Pfam ID B']):
    if B not in pfam_pairs[A]:
        pfam_pairs[A].append(B)
for A, B in zip(threeDID['Pfam ID A'], threeDID['Pfam ID B']):
    if A not in pfam_pairs[B]:
        pfam_pairs[B].append(A)


# ## Find indices of sites for pairs

# In[18]:


# create new columns in biogrid for domain positions
biogrid['domain_a'] = [[] for _ in range(biogrid.shape[0])]
biogrid['domain_b'] = [[] for _ in range(biogrid.shape[0])]
biogrid['domain_seq_a'] = [[] for _ in range(biogrid.shape[0])]
biogrid['domain_seq_b'] = [[] for _ in range(biogrid.shape[0])]


# In[19]:


total_domains_found = 0
total_found = 0

print("Searching for domain pairs")
for PA, PB, DA, DB, DSA, DSB in zip(biogrid['SWISS-PROT Accessions Interactor A'], biogrid['SWISS-PROT Accessions Interactor B'], biogrid['domain_a'], biogrid['domain_b'], biogrid['domain_seq_a'], biogrid['domain_seq_b']):
    found = False

    try:
        # Find locations in PFAM with relevant protein sequences
        pfam_indices_A = pfam_indices[PA] 
        pfam_indices_B = pfam_indices[PB]
        # Get iterable of PFAM IDs for domains for each protein
        # print(pfam_indices_A[0])
        pfam_ids_A = pfam['hmm acc'].loc[pfam_indices_A]
        pfam_ids_B = pfam['hmm acc'].loc[pfam_indices_B]
    except KeyError as inst:
#         print("No pfam entry found for protein")
#         print(inst.args)
        continue
    
    # for each domain PFAM id in A
    for pfam_index_A, pfam_id_A in pfam_ids_A.iteritems():
        # for each domain PFAM id in B
        for pfam_index_B, pfam_id_B in pfam_ids_B.iteritems():
            try:
                # get list of all domain interactions with B
                pfam_id_B_pairs = pfam_pairs[pfam_id_B]
            except KeyError as inst:
#                 print("No pairs found for domain " + inst.args[0])
                continue
            
            # if domain A is in the list of interactions for domain B, we have a match
            if pfam_id_A in pfam_id_B_pairs:

                # save domain starts and ends as tuples
                DA.append(pfam_id_A)
                DB.append(pfam_id_B)
                DSA.append((pfam['envelope start'].loc[pfam_index_A], pfam['envelope end'].loc[pfam_index_A]))
                DSB.append((pfam['envelope start'].loc[pfam_index_B], pfam['envelope end'].loc[pfam_index_B]))
                
                total_domains_found += 1
                found = True
    if found:
        total_found += 1

print(f"Total domains found: {total_domains_found}\nTotal pairs found: {total_found}")


# ## Generate masks

# Drop all where no domains are found

# In[20]:


if single_domain_only:
    print("Selecting PPI pairs with only 1 domain pairs")
    biogrid = biogrid[biogrid.domain_a.str.len() == 1].reset_index()
else:
    print("Selecting PPI pairs with any domain pairs")
    biogrid = biogrid[biogrid.domain_a.str.len() >= 1].reset_index()
    


# In[21]:


# Split out into numpy arrays
Uniprot_id_A = biogrid['SWISS-PROT Accessions Interactor A'].to_numpy()
Uniprot_id_B = biogrid['SWISS-PROT Accessions Interactor B'].to_numpy()
domain_pfam_a = biogrid['domain_a'].to_numpy()
domain_pfam_b = biogrid['domain_b'].to_numpy()
positions_a = biogrid['domain_seq_a'].to_numpy()
positions_b = biogrid['domain_seq_b'].to_numpy()


# Load proteome from UniProt

# In[22]:


uniprot_df = pd.read_csv(uniprot_file,
                        sep = "\t", index_col='Entry')
print("Loaded UniProt proteome")


# Generate site masks

# In[23]:


site_masks = []
print("Generating site masks")
for UA, UB, DA, DB, PA, PB in zip(Uniprot_id_A, Uniprot_id_B, domain_pfam_a, domain_pfam_b, positions_a, positions_b):
    interaction_area = 0
    try:
        lengths = (uniprot_df['Length'].loc[UA], uniprot_df['Length'].loc[UB])
        
        # filter for smaller interaction areas - we don't want to do the whole protein
        if interaction_area <= max_interaction_area * lengths[0]*lengths[1]:
            # initialize mask with dimensions of protein sequence
            mask = np.zeros(lengths, dtype=int)
            
            for pos_A, pos_B in zip(PA, PB):
                
                # calculate area of interaction
                interaction_area += (pos_A[1]-pos_A[0])*(pos_B[1]-pos_B[0])
                # set area of interaction to 1
                mask[(pos_A[0]-1):(pos_A[1]-1), (pos_B[0]-1):(pos_B[1]-1)] = 1 # -1 as protein indexing starts by 1
        
            site_masks.append(mask)
        else:
            site_masks.append(np.NaN)
    except KeyError as inst:
        print(UA, UB)
        print(f"No uniprot entry found for protein {inst.args}")
        site_masks.append(np.NaN)
        
site_masks = np.asarray(site_masks)


# In[24]:


# Turn back into pandas dataframe
masks_domainsOnly = pd.DataFrame({'Uniprot ID A': Uniprot_id_A,
                                            'Uniprot ID B': Uniprot_id_B,
                                            'Domain_id_a': domain_pfam_a,
                                            'Domain_id_b': domain_pfam_b,
                                            'Domain positions A': positions_a,
                                            'Domain positions B': positions_b,
                                            'Sites Masks': site_masks})
# drop all masks for proteins we could not find/area was > 50%
masks_domainsOnly.dropna(inplace=True)
masks_domainsOnly.reset_index(drop=True, inplace=True)

print(f"Created {masks_domainsOnly.shape[0]} domain masks with area <= {max_interaction_area*100} % ")


# In[ ]:


masks_domainsOnly.to_pickle(save_file_path)


# In[ ]:




