import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

silent = False

# In
exp_csv_KEGGDrug_DrugBank = r"..\Exported\exp_csv_KEGGDrug_DrugBank.csv"
exp_csv_interactions = r"..\Exported\exp_csv_interactions.csv"
exp_csv_KEGG_Interaction = r"..\Exported\exp_csv_KEGG_Interaction.csv"
exp_csv_DrugBank_Anvisa = r"..\Exported\exp_csv_DrugBank_Anvisa.csv"
exp_csv_KEGGDrug_Anvisa = r"..\Exported\exp_csv_KEGGDrug_Anvisa.csv"
# Out
exp_csv_KEGGDrug_DrugBank_Analysis = r"..\Exported\exp_csv_KEGGDrug_DrugBank_Analysis.csv"

# Load
df_KEGGDrug_DrugBank = pd.read_csv(exp_csv_KEGGDrug_DrugBank, sep=',')
df_interactions = pd.read_csv(exp_csv_interactions, sep=',')
df_KEGG_interaction = pd.read_csv(exp_csv_KEGG_Interaction, sep=',')
df_DrugBank_Anvisa = pd.read_csv(exp_csv_DrugBank_Anvisa, sep=',')
df_KEGGDrug_Anvisa = pd.read_csv(exp_csv_KEGGDrug_Anvisa, sep=',')


# Coverage Analysis (1st Level)
list_dfs = [df_DrugBank_Anvisa, df_KEGGDrug_Anvisa, df_KEGGDrug_DrugBank]
list_titles = ['Coverage Analysis (1st Level): DrugBank X ANVISA',
			   'Coverage Analysis (1st Level): KEGGDrug X ANVISA', 
			   'Coverage Analysis (1st Level): DrugBank X KEGGDrug']
perf_list = []
for i in range(len(list_dfs)):
	print('### '+list_titles[i]+' ###')
	df_perfect = list_dfs[i][list_dfs[i]['matchingValue']==1.0]
	perf_list.append(df_perfect)
	print('Perfect Match: '+str(len(df_perfect))+' of a total '+str(len(list_dfs[i]))+'.')
	q_index = sum(list_dfs[i]['matchingValue'])/len(list_dfs[i])
	print('Calculated Quality Index: '+str(q_index))
	print()
	# Histogram
	plt.figure(figsize=(12, 12))
	plt.show(block=False)
	plt.hist(list_dfs[i]['matchingValue'], 20)
	plt.xlim(0, 1)
	plt.title(list_titles[i], size=20)
	plt.xlabel("Matching Value", size=20)
	plt.ylabel("Number of Terms", size=20)
	plt.pause(0.01)
df_DrugBank_Anvisa_Perfect = perf_list[0]
df_KEGGDrug_Anvisa_Perfect = perf_list[1]
df_KEGGDrug_DrugBank_Perfect = perf_list[2]
list_dfs = perf_list = [] # free?


# DrugBank Ids With Described Interaction
listDrugBank_interaction = pd.concat([df_interactions['drugbank-id1'], df_interactions['drugbank-id2']], ignore_index=True)
listDrugBank_interaction = listDrugBank_interaction.drop_duplicates()
listDrugBank_interaction.reset_index(inplace=True, drop=True)

# KEGG Ids With Described Interaction
listKEGG_interaction = pd.concat([df_KEGG_interaction['keggdrug-id1'], df_KEGG_interaction['keggdrug-id2']], ignore_index=True)
listKEGG_interaction = listKEGG_interaction.drop_duplicates()
listKEGG_interaction.reset_index(inplace=True, drop=True)


# Coverage Analysis (2nd Level)
df_DrugBank_Anvisa['drugbank-id'].isin(listDrugBank_interaction)
df_reach_interactions = df_DrugBank_Anvisa[df_DrugBank_Anvisa['drugbank-id'].isin(listDrugBank_interaction)]

list_dfs = [df_DrugBank_Anvisa, df_KEGGDrug_Anvisa]
list_dfs_interactions = [listDrugBank_interaction, listKEGG_interaction]
list_col_name = ['drugbank-id', 'keggdrug-id']
list_titles = ['Coverage Analysis (2nd Level): DrugBank X ANVISA Terms Related to Interactions',
			   'Coverage Analysis (2nd Level): KEGGDrug X ANVISA Terms Related to Interactions']
for i in range(len(list_dfs)):
	print('### '+list_titles[i]+' ###')
	df_reachable_interactions = list_dfs[i][list_dfs[i][list_col_name[i]].isin(list_dfs_interactions[i])]
	print('ANVISA can coverage at most: '+str(len(df_reachable_interactions))+' of a total '+str(len(list_dfs_interactions[i]))+'.')
	df_perfect = df_reachable_interactions[df_reachable_interactions['matchingValue']==1.0]
	print('Perfect Match: '+str(len(df_perfect))+' of a total '+str(len(df_reachable_interactions))+'.')
	q_index = sum(df_reachable_interactions['matchingValue'])/len(df_reachable_interactions)
	print('Calculated Quality Index: '+str(q_index))
	print()
	# Histogram
	plt.figure(figsize=(12, 12))
	plt.show(block=False)
	plt.hist(df_reachable_interactions['matchingValue'], 20)
	plt.xlim(0, 1)
	plt.title(list_titles[i], size=20)
	plt.xlabel("Matching Value", size=20)
	plt.ylabel("Number of Terms", size=20)
	plt.pause(0.01)
list_dfs = list_dfs_interactions = [] # free?






# Vector of KEGGDrug Ids
listKEGG = df_KEGGDrug_DrugBank['keggdrug-id']
listKEGG = listKEGG.drop_duplicates()
listKEGG.reset_index(inplace=True, drop=True)
listKEGG = listKEGG.sort_values()

# Maintain only 1 match from KEGG to DrugBank for each KEGG id
for id in listKEGG:

	if not silent: print('KEGG Drug x DrugBank Analysis - Processing KEGG ID: '+str(id)+' \r', end="")

	listIndex = list(df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['keggdrug-id']==id].index)

	if len(listIndex)==1:
		continue
	else:
		max = -1
		maxIndex = -1
		for other in listIndex:
			if df_KEGGDrug_DrugBank.iloc[other,4] > max:
				max = df_KEGGDrug_DrugBank.iloc[other,4]
				maxIndex = other
			if max == 1:
				break
		listIndex.remove(maxIndex)
		df_KEGGDrug_DrugBank.drop(listIndex,inplace=True)
		df_KEGGDrug_DrugBank.reset_index(inplace=True, drop=True)
if not silent: print()

df_KEGGDrug_DrugBank.to_csv(exp_csv_KEGGDrug_DrugBank_Analysis, index = False)


# Check: KEGG Drug Interactions that are presented in DrugBank

df_interactionsAnalysis = df_interactions
df_interactionsAnalysis['onKEGG'] = [0] * len(df_interactionsAnalysis)

df_KEGGDrug_DrugBank_Perfect = df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['matchingValue']==1.0][['keggdrug-id','drugbank-id']]

for i in range(len(df_interactionsAnalysis)):

	if not silent: print('Checking DrugBank Interaction Correspondence with KEGG Drug: '+str(i)+' \r', end="")

	db_id1 = df_interactionsAnalysis.iloc[i,0]
	db_id2 = df_interactionsAnalysis.iloc[i,1]

	kd_id1_candidates = df_KEGGDrug_DrugBank_Perfect[df_KEGGDrug_DrugBank_Perfect['drugbank-id']==db_id1]['keggdrug-id']
	kd_id2_candidates = df_KEGGDrug_DrugBank_Perfect[df_KEGGDrug_DrugBank_Perfect['drugbank-id']==db_id2]['keggdrug-id']

	if len(kd_id1_candidates) < 1 or len(kd_id2_candidates) < 1:
		continue

	if len(kd_id1_candidates) > 1:
		kd_id1 = kd_id1_candidates[kd_id1_candidates['matchingValue']==kd_id1_candidates['matchingValue'].max()].iloc[0]
	else:
		kd_id1 = kd_id1_candidates.iloc[0]

	if len(kd_id2_candidates) > 1:
		kd_id2 = kd_id2_candidates[kd_id2_candidates['matchingValue']==kd_id2_candidates['matchingValue'].max()].iloc[0]
	else:
		kd_id2 = kd_id2_candidates.iloc[0]

	interaction = df_KEGG_interaction[(df_KEGG_interaction['keggdrug-id1'] == kd_id1[0]) & (df_KEGG_interaction['keggdrug-id2'] == kd_id2[0]) | 
								   (df_KEGG_interaction['keggdrug-id1'] == kd_id2[0]) & (df_KEGG_interaction['keggdrug-id2'] == kd_id1[0])]

	if len(interaction) == 0:
		continue
	else:
		df_interactionsAnalysis.iloc[i,2] = (kd_id1[4]+kd_id2[4])/2

if not silent: print()
print('off')



# Check node "trust"
count = 0
for id in listKEGG_interaction:
	matchValue = df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['keggdrug-id']==id].iloc[0,4]
	if matchValue < 1:
		count = count + 1
print('void')


plt.show(block=True) # Deals with block = False (otherwise figures become unresponsive)
print('end') 
