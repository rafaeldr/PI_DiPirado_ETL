import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utils

silent = False
showFigures = True

# Parameters
runCoverageInteractions = True
runIntersectionInteractions = True

# Time Tracker
exp_csv_ComputingTime = r"..\Exported\exp_csv_ComputingTime.csv"
timeTracker = utils.TimeTracker(exp_csv_ComputingTime)
strSubject = 'Initializing Data Analysis Module'
timeTracker.note(strSubject,'start')

# In
exp_csv_KEGGDrug_DrugBank = r"..\Exported\exp_csv_KEGGDrug_DrugBank.csv"
exp_csv_interactions = r"..\Exported\exp_csv_interactions.csv"
exp_csv_KEGG_Interaction = r"..\Exported\exp_csv_KEGG_Interaction.csv"
exp_csv_DrugBank_Anvisa = r"..\Exported\exp_csv_DrugBank_Anvisa.csv"
exp_csv_KEGGDrug_Anvisa = r"..\Exported\exp_csv_KEGGDrug_Anvisa.csv"
# Out
exp_csv_KEGGDrug_DrugBank_Analysis = r"..\Exported\Analysis\exp_csv_KEGGDrug_DrugBank_Analysis.csv"
exp_csv_DrugBank_Interaction_Analysis = r"..\Exported\Analysis\exp_csv_DrugBank_Interaction_Analysis.csv"
exp_csv_KEGGDrug_Interaction_Analysis = r"..\Exported\Analysis\exp_csv_KEGGDrug_Interaction_Analysis.csv"
exp_csv_DrugBank_KEGGDrug_Intersection_Analysis = r"..\Exported\Analysis\exp_csv_DrugBank_KEGGDrug_Intersection_Analysis.csv"
exp_csv_KEGGDrug_DrugBank_Intersection_Analysis = r"..\Exported\Analysis\exp_csv_KEGGDrug_DrugBank_Intersection_Analysis.csv"

# Load
df_KEGGDrug_DrugBank = pd.read_csv(exp_csv_KEGGDrug_DrugBank, sep=',')
df_interactions = pd.read_csv(exp_csv_interactions, sep=',')
df_KEGG_interaction = pd.read_csv(exp_csv_KEGG_Interaction, sep=',')
df_DrugBank_Anvisa = pd.read_csv(exp_csv_DrugBank_Anvisa, sep=',')
df_KEGGDrug_Anvisa = pd.read_csv(exp_csv_KEGGDrug_Anvisa, sep=',')

# Conditional Load
if not runCoverageInteractions:
	df_interactionsCoverageDrugBank = pd.read_csv(exp_csv_DrugBank_Interaction_Analysis, sep=',')
	df_interactionsCoverageKEGGDrug = pd.read_csv(exp_csv_KEGGDrug_Interaction_Analysis, sep=',')
if not runIntersectionInteractions:
	df_intersectionsDrugBank = pd.read_csv(exp_csv_DrugBank_KEGGDrug_Intersection_Analysis, sep=',')
	df_intersectionsKEGGDrug = pd.read_csv(exp_csv_KEGGDrug_DrugBank_Intersection_Analysis, sep=',')

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
	if showFigures:
		plt.figure(figsize=(12, 12))
		plt.show(block=False)
		plt.hist(list_dfs[i]['matchingValue'], 20, facecolor=np.random.rand(3,))
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


# Coverage Analysis (2nd Level) - Terms/Names
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
	if showFigures:
		plt.figure(figsize=(12, 12))
		plt.show(block=False)
		plt.hist(df_reachable_interactions['matchingValue'], 20, facecolor=np.random.rand(3,))
		plt.xlim(0, 1)
		plt.title(list_titles[i], size=20)
		plt.xlabel("Matching Value", size=20)
		plt.ylabel("Number of Terms", size=20)
		plt.pause(0.01)
list_dfs = list_dfs_interactions = [] # free?

# Coverage Analysis (2nd Level) - Interactions
print('### Coverage Analysis (2nd Level) - Interactions ###')

timeTracker.note(strSubject,'end')
strSubject = 'Data Analysis Module - DrugBank Interactions Score Index Calculation'
timeTracker.note(strSubject,'start')

# DrugBank Interactions Score Index Calculation
if runCoverageInteractions:
	df_interactionsCoverageDrugBank = df_interactions
	df_interactionsCoverageDrugBank['matchingValue'] = 0.0
	for index, row in df_interactionsCoverageDrugBank.iterrows():
		if not silent: print('Checking Interaction: '+str(index+1)+' of ' + str(len(df_interactionsCoverageDrugBank)) + '\r', end="")
		value1 = df_DrugBank_Anvisa[df_DrugBank_Anvisa['drugbank-id']==row['drugbank-id1']]['matchingValue']
		value1 = float(np.mean(value1)) if len(value1)>0 else 0
		value2 = df_DrugBank_Anvisa[df_DrugBank_Anvisa['drugbank-id']==row['drugbank-id2']]['matchingValue']
		value2 = float(np.mean(value2)) if len(value2)>0 else 0
		#row['matchingValue'] = float(value1 * value2) # BUG: read-only
		df_interactionsCoverageDrugBank.at[index,'matchingValue'] = float(value1 * value2)
	if not silent: print()
	df_interactionsCoverageDrugBank.to_csv(exp_csv_DrugBank_Interaction_Analysis, index = False)

timeTracker.note(strSubject,'end')
strSubject = 'Data Analysis Module - KEGGDrug Interactions Score Index Calculation'
timeTracker.note(strSubject,'start')

# KEGGDrug Interactions Score Index Calculation
if runCoverageInteractions:
	df_interactionsCoverageKEGGDrug = df_KEGG_interaction
	df_interactionsCoverageKEGGDrug['matchingValue'] = 0.0
	for index, row in df_interactionsCoverageKEGGDrug.iterrows():
		if not silent: print('Checking Interaction: '+str(index+1)+' of ' + str(len(df_interactionsCoverageKEGGDrug)) + '\r', end="")
		value1 = df_KEGGDrug_Anvisa[df_KEGGDrug_Anvisa['keggdrug-id']==row['keggdrug-id1']]['matchingValue']
		value1 = float(np.mean(value1)) if len(value1)>0 else 0
		value2 = df_KEGGDrug_Anvisa[df_KEGGDrug_Anvisa['keggdrug-id']==row['keggdrug-id2']]['matchingValue']
		value2 = float(np.mean(value2)) if len(value2)>0 else 0
		df_interactionsCoverageKEGGDrug.at[index,'matchingValue'] = float(value1 * value2)
	if not silent: print()
	df_interactionsCoverageKEGGDrug.to_csv(exp_csv_KEGGDrug_Interaction_Analysis, index = False)

# Histograms
list_dfs = [df_interactionsCoverageDrugBank, df_interactionsCoverageKEGGDrug]
list_titles = ['Interactions - DrugBank (Coverage through ANVISA)',
			   'Interactions - KEGGDrug (Coverage through ANVISA)']
for i in range(len(list_dfs)):
	print('### '+list_titles[i]+' ###')
	df_perfect = list_dfs[i][list_dfs[i]['matchingValue']==1.0]
	print('Perfect Match: '+str(len(df_perfect))+' of a total '+str(len(list_dfs[i]))+'.')
	q_index = sum(list_dfs[i]['matchingValue'])/len(list_dfs[i])
	print('Calculated Quality Index: '+str(q_index))
	print()
	# Histogram
	if showFigures:
		plt.figure(figsize=(12, 12))
		plt.show(block=False)
		plt.hist(list_dfs[i]['matchingValue'], 20, facecolor=np.random.rand(3,))
		plt.xlim(0, 1)
		plt.title(list_titles[i], size=20)
		plt.xlabel("Matching Value", size=20)
		plt.ylabel("Interactions", size=20)
		plt.pause(0.01)
list_dfs = list_dfs_interactions = [] # free?

timeTracker.note(strSubject,'end')
strSubject = 'Data Analysis Module - KEGGDrug to DrugBank Pairing Filter'
timeTracker.note(strSubject,'start')

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

timeTracker.note(strSubject,'end')
strSubject = 'Data Analysis Module - Check KEGGDrug Interactions that are presented in DrugBank'
timeTracker.note(strSubject,'start')

# Check: KEGG Drug Interactions that are presented in DrugBank
if runIntersectionInteractions:
	df_intersectionsDrugBank = df_interactions
	df_intersectionsDrugBank['onKEGG'] = 0.0
	
	for i in range(len(df_intersectionsDrugBank)):
		if not silent: print('Checking DrugBank Interaction Correspondence with KEGG Drug: '+str(i+1)+' of ' + str(len(df_intersectionsDrugBank)) +' \r', end="")

		db_id1 = df_intersectionsDrugBank.iloc[i,0]
		db_id2 = df_intersectionsDrugBank.iloc[i,1]

		kd_id1_candidates = df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['drugbank-id']==db_id1]
		kd_id2_candidates = df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['drugbank-id']==db_id2]

		if len(kd_id1_candidates) < 1 or len(kd_id2_candidates) < 1:
			continue

		kd_id1 = kd_id1_candidates[kd_id1_candidates['matchingValue']==kd_id1_candidates['matchingValue'].max()].iloc[0]
		kd_id2 = kd_id2_candidates[kd_id2_candidates['matchingValue']==kd_id2_candidates['matchingValue'].max()].iloc[0]

		interaction = df_KEGG_interaction[(df_KEGG_interaction['keggdrug-id1'] == kd_id1[0]) & (df_KEGG_interaction['keggdrug-id2'] == kd_id2[0]) | 
									   (df_KEGG_interaction['keggdrug-id1'] == kd_id2[0]) & (df_KEGG_interaction['keggdrug-id2'] == kd_id1[0])]

		if len(interaction) == 0:
			continue
		else:
			df_intersectionsDrugBank.iloc[i,2] = float(kd_id1[4]*kd_id2[4])
	if not silent: print()
	df_intersectionsDrugBank.to_csv(exp_csv_DrugBank_KEGGDrug_Intersection_Analysis, index = False)

timeTracker.note(strSubject,'end')
strSubject = 'Data Analysis Module - Check DrugBank Interactions that are presented in KEGGDrug'
timeTracker.note(strSubject,'start')

# Check: DrugBank Interactions that are presented in KEGG Drug
if runIntersectionInteractions:
	df_intersectionsKEGGDrug = df_KEGG_interaction
	df_intersectionsKEGGDrug['onDrugBank'] = 0.0

	for i in range(len(df_intersectionsKEGGDrug)):
		if not silent: print('Checking KEGGDrug Interaction Correspondence with DrugBank: '+str(i+1)+' of ' + str(len(df_intersectionsKEGGDrug)) +' \r', end="")

		kd_id1 = df_intersectionsKEGGDrug.iloc[i,0]
		kd_id2 = df_intersectionsKEGGDrug.iloc[i,1]

		db_id1_candidates = df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['keggdrug-id']==kd_id1]
		db_id2_candidates = df_KEGGDrug_DrugBank[df_KEGGDrug_DrugBank['keggdrug-id']==kd_id2]

		if len(db_id1_candidates) < 1 or len(db_id2_candidates) < 1:
			continue
		if len(db_id1_candidates) > 1 or len(db_id2_candidates) > 1:
			print('Unexpected error: KEGGDrug pairing with DrugBank returned more than 1 id.')
			exit(1)

		db_id1 = int(db_id1_candidates['drugbank-id'])
		db_id2 = int(db_id2_candidates['drugbank-id'])

		interaction = df_interactions[(df_interactions['drugbank-id1'] == db_id1) & (df_interactions['drugbank-id2'] == db_id2) | 
									   (df_interactions['drugbank-id1'] == db_id2) & (df_interactions['drugbank-id2'] == db_id1)]

		if len(interaction) == 0:
			continue
		else:
			db_id1_matchValue = float(db_id1_candidates['matchingValue'])
			db_id2_matchValue = float(db_id2_candidates['matchingValue'])
			df_intersectionsKEGGDrug.iloc[i,2] = float(db_id1_matchValue*db_id2_matchValue)
	if not silent: print()
	df_intersectionsKEGGDrug.to_csv(exp_csv_KEGGDrug_DrugBank_Intersection_Analysis, index = False)

timeTracker.note(strSubject,'end')

timeTracker.export() # Remove when calling from ETL main module!
if showFigures:
	print('Done! Waiting for figures to be closed...')
	plt.show(block=True) # Deals with block = False (otherwise figures become unresponsive)
print('End of Data Analysis Module') 