import pandas as pd
import numpy as np

silent = False

exp_csv_KEGGDrug_DrugBank = r"..\Exported\exp_csv_KEGGDrug_DrugBank.csv"
exp_csv_KEGGDrug_DrugBank_Analysis = r"..\Exported\exp_csv_KEGGDrug_DrugBank_Analysis.csv"
exp_csv_interactions = r"..\Exported\exp_csv_interactions.csv"
exp_csv_KEGG_Interaction = r"..\Exported\exp_csv_KEGG_Interaction.csv"

# Load
df_KEGGDrug_DrugBank = pd.read_csv(exp_csv_KEGGDrug_DrugBank, sep=',')
df_interactions = pd.read_csv(exp_csv_interactions, sep=',')
df_KEGG_interaction = pd.read_csv(exp_csv_KEGG_Interaction, sep=',')


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


# KEGG Ids With Described Interaction
listKEGG_interaction = pd.concat([df_KEGG_interaction['keggdrug-id1'], df_KEGG_interaction['keggdrug-id2']], ignore_index=True)
listKEGG_interaction = listKEGG_interaction.drop_duplicates()
listKEGG_interaction.reset_index(inplace=True, drop=True)

# DrugBank Ids With Described Interaction
listDrugBank_interaction = pd.concat([df_interactions['drugbank-id1'], df_interactions['drugbank-id2']], ignore_index=True)
listDrugBank_interaction = listDrugBank_interaction.drop_duplicates()
listDrugBank_interaction.reset_index(inplace=True, drop=True)

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



print('end') 