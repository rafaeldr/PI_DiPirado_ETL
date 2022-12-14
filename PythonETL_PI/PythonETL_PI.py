import xmltodict
import pandas as pd
import numpy as np
import collections
import unidecode
import TranslationModule
import TermMatchingModule
import SQLModule as sql
import KEGGDrugModule
import utils
import os

# Parameters
callTranslator = False  # Keep false unless required (implies in costs from GoogleCloud)
callTermMatchingDrugBank = True # Keep false unless required (implies in high computation time) [Requires exp_csv_DrugBank_Anvisa]
callTermMatchingKEGGDrug = True # Keep false unless required (implies in high computation time) [Requires exp_csv_KEGGDrug_Anvisa]
callKEGGDrugModule = True
callTermMatchingKEGGDrugDrugBank = True # Keep false unless required (implies in high computation time) [Requires exp_csv_KEGGDrug_DrugBank]
mode_DrugBank_OR_KEGGDrug = 'kegg' # for KEGG Drug: 'kegg', for DrugBank: 'drugbank'  (imply in BigTable only)
prodEnvironment = False # False for "development/test"; true for "production" execution
silent = False          # Display track of progress info (when False)
TermMatchingModule.silent = silent
TranslationModule.silent = silent

# File locations
if prodEnvironment:
    drugbank_file = r"..\DataSources\full database.xml"    # Full DrugBank - Takes ~6 min to parse
    KEGGDrugModule.file = r"..\DataSources\drug"
    anvisa_file = r"..\DataSources\data_anvisa.csv"
    exp_csv_pAtivos_Traducoes = r"..\Exported\exp_csv_pAtivos_Traducoes.csv"
else:
    drugbank_file = r"..\DataSources\drugbank_sample6.xml"  # Sample DrugBank with only 6 drugs
    KEGGDrugModule.file = r"..\DataSources\drug_sample"
    anvisa_file = r"..\DataSources\data_anvisa_sample.csv"
    exp_csv_pAtivos_Traducoes = r"..\DataSources\exp_csv_pAtivos_Traducoes_sample.csv"
# Export
isExist = os.path.exists(r"..\Scripts")
if not isExist: os.makedirs(r"..\Scripts")
isExist = os.path.exists(r"..\Exported")
if not isExist: os.makedirs(r"..\Exported")
exp_csv_drugs = r"..\Exported\exp_csv_drugs.csv"
exp_csv_drugsSynonyms = r"..\Exported\exp_csv_drugsSynonyms.csv"
exp_csv_interactions = r"..\Exported\exp_csv_interactions.csv"
exp_csv_pAtivos = r"..\Exported\exp_csv_pAtivos.csv"
exp_csv_pAtivosAccented = r"..\Exported\exp_csv_pAtivosAccented.csv"
exp_csv_Nomes = r"..\Exported\exp_csv_Nomes.csv"
exp_csv_Nomes_pAtivos = r"..\Exported\exp_csv_Nomes_pAtivos.csv"
exp_csv_Analysis_Nomes_pAtivos = r"..\Exported\exp_csv_Analysis_Nomes_pAtivos.csv"
exp_csv_DrugBank_Anvisa = r"..\Exported\exp_csv_DrugBank_Anvisa.csv" # REQUIRED When callTermMatchingDrugBank is False
exp_csv_KEGG_Drugs = r"..\Exported\exp_csv_KEGG_Drugs.csv"
exp_csv_KEGG_DrugsSynonyms = r"..\Exported\exp_csv_KEGG_DrugsSynonyms.csv"
exp_csv_KEGG_Interaction = r"..\Exported\exp_csv_KEGG_Interaction.csv"
exp_csv_KEGGDrug_Anvisa = r"..\Exported\exp_csv_KEGGDrug_Anvisa.csv" # REQUIRED When callTermMatchingKEGGDrug is False
exp_csv_KEGGDrug_DrugBank = r"..\Exported\exp_csv_KEGGDrug_DrugBank.csv" # REQUIRED When callTermMatchingKEGGDrugDrugBank is False
exp_csv_ComputingTime = r"..\Exported\exp_csv_ComputingTime.csv"

timeTracker = utils.TimeTracker(exp_csv_ComputingTime)
KEGGDrugModule.timeTracker = timeTracker

# Importing Data Sources (AS-IS) - ANVISA
strSubject = 'Importing Data Sources (AS-IS) - ANVISA'
if not silent: print(strSubject)
timeTracker.note(strSubject,'start')
df_anvisa = pd.read_csv(anvisa_file, sep=';')
timeTracker.note(strSubject,'end')
# Importing Data Sources (AS-IS) - DrugBank
strSubject = 'Importing Data Sources (AS-IS) - DrugBank'
if not silent: print(strSubject) 
timeTracker.note(strSubject,'start')
with open(drugbank_file, encoding='utf-8') as fd:
    drugbank_dict = xmltodict.parse(fd.read())
timeTracker.note(strSubject,'end')

# region DrugBank
strSubject = 'DrugBank - Data Processing (Transform)'
timeTracker.note(strSubject,'start')

# Drugs - Extraction
if not silent: print('DrugBank - Extracting Names') 
data = []
dataSynonyms = []
for drug in drugbank_dict['drugbank']['drug']:
    if type(drug['drugbank-id']) == list:
        id = str(drug['drugbank-id'][0]['#text'])
    else:
        id = str(drug['drugbank-id']['#text'])
    id = id[2:] # Adjust ID to integer (not varchar/string)
    data.append({'drugbank-id':int(id),
                 'name':str(drug['name']).strip().upper()
                 })

    # Drugs - Synonyms
    if drug['synonyms']!=None:
        if type(drug['synonyms']['synonym']) == collections.OrderedDict: # only 1 registry
            synonymName = str(drug['synonyms']['synonym']['#text']).strip().upper()
            dataSynonyms.append({'drugbank-id':int(id),
                         'name':synonymName
                         })
        elif type(drug['synonyms']['synonym']) == list:
            for synonym in drug['synonyms']['synonym']:
                synonymName = str(synonym['#text']).strip().upper()
                dataSynonyms.append({'drugbank-id':int(id),
                             'name':synonymName
                             })
        else:
            print('Unexpected error: Unexpected condition during DrugBank Synonyms extraction.')
            exit(1)

df_drugs = pd.DataFrame(data)
df_drugsSynonyms = pd.DataFrame(dataSynonyms)

# Interactions - Extraction
if not silent: print('DrugBank - Extracting Interactions') 
data = []
for drugOrigin in drugbank_dict['drugbank']['drug']:
    if type(drugOrigin['drugbank-id']) == list:
        drugOrigin_id = str(drugOrigin['drugbank-id'][0]['#text'])
    else:
        drugOrigin_id = str(drugOrigin['drugbank-id']['#text'])

    if drugOrigin['drug-interactions']!=None:
        if type(drugOrigin['drug-interactions']['drug-interaction']) == collections.OrderedDict: # only 1 registry
            drugDestiny_id = str(drugOrigin['drug-interactions']['drug-interaction']['drugbank-id'])
            data.append([int(drugOrigin_id[2:]),
                         int(drugDestiny_id[2:])
                         ])
# Changed for execution time improvement
#            data.append({'drugbank-id1':drugOrigin_id,
#                         'drugbank-id2':drugDestiny_id})
        elif type(drugOrigin['drug-interactions']['drug-interaction']) == list:
            for drugDestiny in drugOrigin['drug-interactions']['drug-interaction']:
                drugDestiny_id = str(drugDestiny['drugbank-id'])
                data.append([int(drugOrigin_id[2:]),
                             int(drugDestiny_id[2:])
                             ])
        else:
            print('Unexpected error: Unexpected condition during DrugBank drug-interaction extraction.')
            exit(1)

# Removing reversed duplicates
if not silent: print('DrugBank - Interactions - Removing Reversed Duplicates') 
data = {tuple(sorted(item)) for item in data}
df_interactions = pd.DataFrame(data, columns=['drugbank-id1','drugbank-id2'])

# Test Constraint id1 < id2 always
if sum(df_interactions['drugbank-id1']<df_interactions['drugbank-id2']) != len(df_interactions):
    print("Unexpected error: DrugBank interactions should respect 'id1 < id2'.")
    exit(1)

# Check Referential Integrity
diff_set = set(set(df_interactions['drugbank-id1']).union(df_interactions['drugbank-id2'])).difference(df_drugs['drugbank-id'])
if len(diff_set) > 0:
    if not silent: print('Warning: Found '+str(len(diff_set))+' DrugBank ids that do not respect referential integrity.')
    if not silent: print('         These are '+str(diff_set)[1:-1]+'.')
    if not silent: print('         Trying to force the inclusion of these ids!')
    if not silent: print('         An extra computing power that could be avoided.')
    
    # Force Adjustments
    found = len(diff_set)
    counter = 0
    for drugOrigin in drugbank_dict['drugbank']['drug']:
        if drugOrigin['drug-interactions']!=None:
            if drugOrigin['drug-interactions']['drug-interaction'] == collections.OrderedDict: # only 1 registry
                drugDestiny_id = str(drugOrigin['drug-interactions']['drug-interaction']['drugbank-id'])
                drugDestiny_id = int(drugDestiny_id[2:])
                if drugDestiny_id in diff_set:
                    drugDestiny_name = str(drugOrigin['drug-interactions']['drug-interaction']['name'])
                    df_drugs = df_drugs.append({'drugbank-id': drugDestiny_id, 'name': drugDestiny_name.strip().upper()}, ignore_index=True)
                    diff_set.remove(drugDestiny_id)
                    counter += 1
            elif type(drugOrigin['drug-interactions']['drug-interaction']) == list:
                for drugDestiny in drugOrigin['drug-interactions']['drug-interaction']:
                    drugDestiny_id = str(drugDestiny['drugbank-id'])
                    drugDestiny_id = int(drugDestiny_id[2:])
                    if drugDestiny_id in diff_set:
                        drugDestiny_name = str(drugDestiny['name'])
                        df_drugs = df_drugs.append({'drugbank-id': drugDestiny_id, 'name': drugDestiny_name.strip().upper()}, ignore_index=True)
                        diff_set.remove(drugDestiny_id)
                        counter += 1
    
    if counter<found or len(diff_set)>0:
        print("Unexpected error: Some ids referenced by interactions don't exist as DrugBank drug names.")
        exit(1)

# DrugBank Drug Synonyms (Concat)
df_drugsSynonyms = pd.concat([df_drugs, df_drugsSynonyms], ignore_index=True)
df_drugsSynonyms = df_drugsSynonyms.drop_duplicates()
df_drugsSynonyms.reset_index(inplace=True, drop=True)

# Exporting
if not silent: print('DrugBank - Exporting CSVs') 
df_drugs.to_csv(exp_csv_drugs, index = False)
df_drugsSynonyms.to_csv(exp_csv_drugsSynonyms, index = False)
df_interactions.to_csv(exp_csv_interactions, index = False)

timeTracker.note(strSubject,'end')
# endregion DrugBank

# region ANVISA

strSubject = 'ANVISA - Data Processing (Transform)'
timeTracker.note(strSubject,'start')

# Names (Anvisa_Name)

s_Names = df_anvisa['NOME_PRODUTO']  # Series
s_Registry = df_anvisa['NUMERO_REGISTRO_PRODUTO']  
s_pAtivos = df_anvisa['PRINCIPIO_ATIVO']  

# Test
if len(s_Names) != len(s_Registry) or len(s_Names) != len(s_pAtivos):
    print('Unexpected error: Data Series s_names and s_Registry and s_pAtivos differ in length!')
    exit(1)

dictNames = dict()                  # ANVISA Commercial Names
dictNamesAccented = dict()                  # Same WITH accentuation
listRegistry = list()               # ANVISA Registry Number (optional/debugging purposes)
dictPrinciples = dict()             # ANVISA Active Principles
dictPrinciplesAccented = dict()             # Same WITH accentuation
list_Names_Principles = list()      # ANVISA Relation (Commercial Names <--> Active Principles)
new_id_name = 0
new_id_principle = 0
for i in range(len(s_Names)):
    # Names & Registry
    if not silent: print('ANVISA - Processing Names and Active Principles: '+str(i+1)+' of '+str(len(s_Names))+'\r', end="")
    name_accented = str(s_Names[i]).strip().upper()
    name_accented = " ".join(name_accented.split())  # Normalize White Spaces
    name = unidecode.unidecode(name_accented)
    if name not in dictNames:
        new_id_name += 1
        dictNames[name] = new_id_name
        dictNamesAccented[name_accented] = new_id_name
        listRegistry.append((new_id_name, int(s_Registry[i]))) # append tuple in a list (previous bug: "extend" flat the elements)
    else:
        listRegistry.append((int(dictNames[name]), int(s_Registry[i])))

    # Extracting Active Principles
    rowStr = s_pAtivos[i]
    if type(rowStr)==str:
        principleList = list(map(str.upper,list(map(str.strip, rowStr.split('+')))))
        
        for principle in principleList:
            principle = " ".join(principle.split())  # Normalize White Spaces
            principle_u = unidecode.unidecode(principle)
            # Active Principles Entity (Keep Unicity)
            if principle_u not in dictPrinciples:
                new_id_principle += 1
                dictPrinciples[principle_u] = new_id_principle
                dictPrinciplesAccented[principle] = new_id_principle
                # Active Principles - Relation with Name
                list_Names_Principles.append((int(dictNames[name]), new_id_principle)) # Here name is always on its dict  // Tuple inside List
            else:
                if (int(dictNames[name]), int(dictPrinciples[principle_u])) not in list_Names_Principles:
                    list_Names_Principles.append((int(dictNames[name]), int(dictPrinciples[principle_u]))) # Same
    else:
        continue # just ignore
if not silent: print()


# Search for Products With Exact Same Name of Action Principles (Analysis Task)
list_Equal_Names_Principles = list()
for i in reversed(range(len(list(dictNames.keys())))): # Reversed cause size changes over iterations
    if not silent: print('ANVISA - Analyzing Names and Active Principles: (reversed) '+str(i+1)+' of '+str(len(list(dictNames.keys())))+'     \r', end="")
    nameStr = list(dictNames.keys())[i]
    nameStrList = list(map(str.upper,list(map(str.strip, nameStr.split('+')))))
    if len(nameStrList) == 1:
        # Case A : Product Name = 1 Active Principle    
        if nameStrList[0] in dictPrinciples:

            idx_number = dictNames[nameStrList[0]]
            del dictNames[nameStrList[0]] # Step 1: Remove from dictNames
            key_value = list(dictNamesAccented.keys())[list(dictNamesAccented.values()).index(idx_number)] # Required since key_value can have different accentuation
            del dictNamesAccented[key_value] # Step 2: Remove from dictNamesAccented

            for item in reversed(list_Names_Principles):
                if item[0] == idx_number:
                    list_Names_Principles.remove(item)
    # Case B : Product Name = 2-More Active Principles
    elif set(nameStrList).issubset(dictPrinciples): # Only the case where ALL exists as principles
        idx_number = dictNames[nameStr] # original full name
        del dictNames[nameStr]  # remove Name (composed)
        key_value = list(dictNamesAccented.keys())[list(dictNamesAccented.values()).index(idx_number)] 
        del dictNamesAccented[key_value] 

        for item in reversed(list_Names_Principles):
            if item[0] == idx_number:       # Remember: Can be already removed during Case A?
                list_Names_Principles.remove(item)
    else:
        list_Equal_Names_Principles.append((int(dictNames[nameStr]), nameStr))  # No match (but multiple names)
if not silent: print()

# Manual Cleanup Section (after visual inspection)

if prodEnvironment:
    if not silent: print('ANVISA - Executing Cleanup Procedures') 
    # 1. Renaming
    dictPrinciples['HIDROBROMETO DE CITALOPRAM'] = dictPrinciples.pop('HIDROBROMETO DE CITALOPRAM (PORT. 344/98 LISTA C 1)')
    dictPrinciplesAccented['HIDROBROMETO DE CITALOPRAM'] = dictPrinciplesAccented.pop('HIDROBROMETO DE CITALOPRAM (PORT. 344/98 LISTA C 1)')
    dictPrinciples['OXANDROLONA'] = dictPrinciples.pop('OXANDROLONA (PORT. 344/98 LISTA C 5)')
    dictPrinciplesAccented['OXANDROLONA'] = dictPrinciplesAccented.pop('OXANDROLONA (PORT. 344/98 LISTA C 5)')

    # 2. Change Codes and Delete
    change_dict = {3739 : 23, 2974 : 23, 1331 : 23, 1843 : 1873, 1299 : 639, 571 : 570,
                   3872 : 2962, 1180 : 211, 2610 : 2450, 3404 : 1761, 1387 : 414, 2603 : 768, 
                   2723 : 1994, 2561 : 17, 2727 : 311, 1974 : 19, 3561 : 105, 3821 : 105, 3210 : 2924, 
                   3215 : 1239, 3782 : 169, 1837 : 37, 2928 : 20, 1424 : 20, 2716 : 1256, 3417 : 1256, 
                   9 : 1265, 3061 : 1938, 1714 : 1487, 2607 : 416, 2040 : 501, 418 : 501, 3687 : 501, 
                   2046 : 3397, 2350 : 2944}

    change_list = [3739, 2974, 1331, 1843, 1299, 571, 3872, 1180, 2610, 3404,
                   1387, 2603, 2723, 2561, 2727, 1974, 3561, 3821, 3210, 3215,
                   3782, 1837, 2928, 1424, 2716, 3417, 9, 3061, 1714, 2607,
                   2040, 418, 3687, 2046, 2350]

    # Adjust relationship
    for item in reversed(list_Names_Principles):
        if item[1] in change_list:  #[0] 'idProduto' [1] 'idPrincipio'
            item_temp = list(item)
            item_temp[1] = change_dict[item_temp[1]]
            item_temp = tuple(item_temp)
            list_Names_Principles.remove(item)
            # Needs to check if the "new relation" doesn't already exist (with the other id)
            if item_temp not in list_Names_Principles:
                list_Names_Principles.append(item_temp)

    # Remove principle
    for item in change_list:
        idx_number = item
        key_value = list(dictPrinciples.keys())[list(dictPrinciples.values()).index(idx_number)]
        del dictPrinciples[key_value]
        key_value = list(dictPrinciplesAccented.keys())[list(dictPrinciplesAccented.values()).index(idx_number)]
        del dictPrinciplesAccented[key_value] 

    # 3. Just Delete
    delete_list = [3954, 3878, 3867, 3823, 3799, 3727, 3472, 3396, 3348, 3284, 2751, 
                   2742, 2665, 2562, 2418, 2156, 2115, 1977, 1739, 1610, 936, 632, 1281]

    for item in delete_list:
        idx_number = item
        key_value = list(dictPrinciples.keys())[list(dictPrinciples.values()).index(idx_number)]
        del dictPrinciples[key_value]
        key_value = list(dictPrinciplesAccented.keys())[list(dictPrinciplesAccented.values()).index(idx_number)]
        del dictPrinciplesAccented[key_value] 

    # Remove relationship
    for item in reversed(list_Names_Principles):
        if item[1] in delete_list:
            list_Names_Principles.remove(item)

# Prepare DataFrames
if not silent: print('ANVISA - Creating DataFrames') 
df_Anvisa_PrinciplesAccented = pd.DataFrame(dictPrinciplesAccented.items(), columns=['nome_pAtivo','id_pAtivo'])
df_Anvisa_Principles = pd.DataFrame(dictPrinciples.items(), columns=['nome_pAtivo','id_pAtivo'])
df_Anvisa_Names = pd.DataFrame(dictNames.items(), columns=['nomeProduto','id'])
df_Anvisa_Names_Principles = pd.DataFrame(list_Names_Principles, columns=['idProduto','idPrincipio'])
df_Equal_Names_Principles = pd.DataFrame(list_Equal_Names_Principles, columns=['idProduto','nomeProduto'])


# Check Referential Integrity
if len(df_Anvisa_PrinciplesAccented) != len(df_Anvisa_Principles):
    print('Unexpected error: Data Frames for Active Principles (Accented and Clean) do not match in size.')
    exit(1)
# Loose Act. Principle ids?
if len(set(df_Anvisa_Names_Principles['idPrincipio']).difference(df_Anvisa_Principles['id_pAtivo'])) > 0:
    print('Unexpected error: ANVISA Data Frame for Names <-> Active Principles contains wrong Active Principles ids.')
    exit(1)
# Loose Anvisa Product Name ids?
if len(set(df_Anvisa_Names_Principles['idProduto']).difference(df_Anvisa_Names['id'])) > 0:
    print('Unexpected error: ANVISA Data Frame for Names <-> Active Principles contains wrong Product ids.')
    exit(1)

timeTracker.note(strSubject,'end') # ANVISA

# Translation Section Call (Run Once) - "Limited Resource" [Google Translator API]
if callTranslator:
    strSubject = 'ANVISA - Calling Translation Module'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    s_pAtivos = df_Anvisa_PrinciplesAccented['nome_pAtivo']  # Series
    s_pAtivosTranslated = TranslationModule.BatchTranslate(s_pAtivos)
    s_pAtivosTranslated.to_csv(exp_csv_pAtivos_Traducoes, index = False)
    df_PrinciplesTraducoes = s_pAtivosTranslated
    timeTracker.note(strSubject,'end')
else:
    strSubject = 'ANVISA - Loading Translations'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    df_PrinciplesTraducoes = pd.read_csv(exp_csv_pAtivos_Traducoes, sep=',')
    timeTracker.note(strSubject,'end')


# Adjust Translations for Upper Case
# For some strange (and external) reason, google translator send some terms with mixed case
df_PrinciplesTraducoes['0'] = df_PrinciplesTraducoes['0'].str.upper()

# Test
if len(df_Anvisa_PrinciplesAccented) != len(df_PrinciplesTraducoes):
    print("Unexpected error: Dataframes doesn't match in length! Translations failed or data changed shape since last batch translation.")
    exit(1)


# Adjust DataFrames
df_Anvisa_PrinciplesAccented['translated_pAtivo'] = df_PrinciplesTraducoes
df_Anvisa_Principles['translated_pAtivo'] = df_PrinciplesTraducoes

# Export DataFrames
if not silent: print('ANVISA - Exporting CSVs')
df_Anvisa_PrinciplesAccented.to_csv(exp_csv_pAtivosAccented, encoding="utf-8", index = False)
df_Anvisa_Principles.to_csv(exp_csv_pAtivos, encoding="utf-8", index = False)
df_Anvisa_Names.to_csv(exp_csv_Nomes, encoding="utf-8", index = False)
df_Anvisa_Names_Principles.to_csv(exp_csv_Nomes_pAtivos, encoding="utf-8", index = False)
df_Equal_Names_Principles.to_csv(exp_csv_Analysis_Nomes_pAtivos, encoding="utf-8", index = False)


# endregion

# region KEGG Drug

# Importing KEGG Drug
if callKEGGDrugModule:
    df_KEGG_drugs, df_KEGG_drugsSynonyms, df_KEGG_interaction = KEGGDrugModule.importKEGGDrug()
    df_KEGG_drugs.to_csv(exp_csv_KEGG_Drugs, encoding="utf-8", index = False)
    df_KEGG_drugsSynonyms.to_csv(exp_csv_KEGG_DrugsSynonyms, encoding="utf-8", index = False)
    df_KEGG_interaction.to_csv(exp_csv_KEGG_Interaction, encoding="utf-8", index = False)
else:
    strSubject = 'KEGG Drug - Loading Preprocessed Data'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    df_KEGG_drugs = pd.read_csv(exp_csv_KEGG_Drugs, sep=',')
    df_KEGG_drugsSynonyms = pd.read_csv(exp_csv_KEGG_DrugsSynonyms, sep=',')
    df_KEGG_interaction = pd.read_csv(exp_csv_KEGG_Interaction, sep=',')
    timeTracker.note(strSubject,'end')

# endregion


# Validation: Unicity of Names 
# (does not guarantee semantic unicity)
if len(df_drugs['name']) != len(set(df_drugs['name'])):
    print("Unexpected error: DrugBank Names are not unique!")
    exit(1)
#if len(df_Anvisa_PrinciplesAccented['translated_pAtivo']) != len(set(df_Anvisa_PrinciplesAccented['translated_pAtivo'])):
#    print("Unexpected error: Anvisa Active Principle Names are not unique.")
#    exit(1) # IT FAILS
# -> Slightly different terms in portuguese generate the same (and correct) terminology in english
# --> This will be considered in Pass 1 for term matching


# Nomenclature Pair Matching Module Call - DrugBank
if callTermMatchingDrugBank:
    strSubject = 'DrugBank + ANVISA - Calling Term Matching Module (With Synonyms)'
    if not silent: print(strSubject)
    timeTracker.note(strSubject,'start')
    #df_DrugBank_Anvisa = TermMatchingModule.match(df_drugs['name'], df_drugs['drugbank-id'], df_Anvisa_PrinciplesAccented['translated_pAtivo'], df_Anvisa_PrinciplesAccented['id_pAtivo']) # Order matters
    #df_DrugBank_Anvisa = TermMatchingModule.match(df_Anvisa_PrinciplesAccented['translated_pAtivo'], df_Anvisa_PrinciplesAccented['id_pAtivo'], df_drugs['name'], df_drugs['drugbank-id'], 'drugbank')
    df_DrugBank_Anvisa = TermMatchingModule.match(df_Anvisa_PrinciplesAccented['translated_pAtivo'], df_Anvisa_PrinciplesAccented['id_pAtivo'], df_drugsSynonyms['name'], df_drugsSynonyms['drugbank-id'], 'drugbank')
    df_DrugBank_Anvisa.to_csv(exp_csv_DrugBank_Anvisa, index = False)
    timeTracker.note(strSubject,'end')
else:
    strSubject = 'DrugBank + ANVISA - Loading Preprocessed Term Matching'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    df_DrugBank_Anvisa = pd.read_csv(exp_csv_DrugBank_Anvisa, sep=',')
    timeTracker.note(strSubject,'end')


# Nomenclature Pair Matching Module Call - KEGGDrug
if callTermMatchingKEGGDrug:
    strSubject = 'KEGG Drug + ANVISA - Calling Term Matching Module (With Synonyms)'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    #df_KEGGDrug_Anvisa = TermMatchingModule.match(df_Anvisa_PrinciplesAccented['translated_pAtivo'], df_Anvisa_PrinciplesAccented['id_pAtivo'], df_KEGG_drugs['name'], df_KEGG_drugs['keggdrug-id'], 'keggdrug')
    df_KEGGDrug_Anvisa = TermMatchingModule.match(df_Anvisa_PrinciplesAccented['translated_pAtivo'], df_Anvisa_PrinciplesAccented['id_pAtivo'], df_KEGG_drugsSynonyms['name'], df_KEGG_drugsSynonyms['keggdrug-id'], 'keggdrug')
    df_KEGGDrug_Anvisa.to_csv(exp_csv_KEGGDrug_Anvisa, index = False)
    timeTracker.note(strSubject,'end')
else:
    strSubject = 'KEGG Drug + ANVISA - Loading Preprocessed Term Matching'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    df_KEGGDrug_Anvisa = pd.read_csv(exp_csv_KEGGDrug_Anvisa, sep=',')
    timeTracker.note(strSubject,'end')


# Nomenclature Pair Matching Module Call - KEGGDrug + DrugBank
if callTermMatchingKEGGDrugDrugBank:
    strSubject = 'KEGG Drug + DrugBank - Calling Term Matching Module (With Synonyms)'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    df_KEGGDrug_DrugBank = TermMatchingModule.match(df_KEGG_drugsSynonyms['name'], df_KEGG_drugsSynonyms['keggdrug-id'], df_drugsSynonyms['name'], df_drugsSynonyms['drugbank-id'], 'keggdrug-drugbank')
    df_KEGGDrug_DrugBank.to_csv(exp_csv_KEGGDrug_DrugBank, index = False)
    timeTracker.note(strSubject,'end')
else:
    strSubject = 'KEGG Drug + DrugBank - Loading Preprocessed Term Matching'
    if not silent: print(strSubject) 
    timeTracker.note(strSubject,'start')
    df_KEGGDrug_DrugBank = pd.read_csv(exp_csv_KEGGDrug_DrugBank, sep=',')
    timeTracker.note(strSubject,'end')


# KEGGDrug + DrugBank Interaction Analysis

# TODO


# BigTable Section
strSubject = 'Populating BigTable Data Structure'
if not silent: print(strSubject) 
timeTracker.note(strSubject,'start')

# TODO: Find a way to combine both (Next Release)
if mode_DrugBank_OR_KEGGDrug == 'kegg':
    df_BigTable = df_KEGG_drugs.rename(columns={'keggdrug-id': 'id_principal', 'name' : 'nome'})
    df_BigTable['tipo_origem'] = [1] * len(df_KEGG_drugs)
elif mode_DrugBank_OR_KEGGDrug == 'drugbank':
    df_BigTable = df_drugs.rename(columns={'drugbank-id': 'id_principal', 'name' : 'nome'})
    df_BigTable['tipo_origem'] = [1] * len(df_drugs)
else:
    print('Unexpected error: Unexpected mode for Big Table generation.')
    exit(1)

df_aux = df_Anvisa_Names.rename(columns={'id': 'id_principal', 'nomeProduto' : 'nome'})
df_aux['tipo_origem'] = [2] * len(df_Anvisa_Names)
df_BigTable = df_BigTable.append(df_aux, ignore_index=True)

df_aux = df_Anvisa_PrinciplesAccented.rename(columns={'id_pAtivo': 'id_principal', 'nome_pAtivo' : 'nome'})
df_aux = df_aux[['id_principal','nome']]
df_aux['tipo_origem'] = [3] * len(df_Anvisa_PrinciplesAccented)
df_BigTable = df_BigTable.append(df_aux, ignore_index=True)

timeTracker.note(strSubject,'end')

# SQL Scripts
strSubject = 'Generating All SQL Scripts'
timeTracker.note(strSubject,'start')

# ['drugbank-id', 'name']
if not silent: print('SQL Scripts - Creating CREATE and INSERT scripts for database tables')
df_drugs.rename(columns={'drugbank-id': 'drugbank_id'}, inplace = True)
sqlDrugBank_Name = sql.SQLScripting(df_drugs,'DrugBank_Nome', [], [], ['drugbank_id'])
# ['drugbank-id1', 'drugbank-id2']
df_interactions.rename(columns={'drugbank-id1': 'drugbank_id1', 'drugbank-id2': 'drugbank_id2'}, inplace = True)
sqlDrugBank_Interactions = sql.SQLScripting(df_interactions, 'DrugBank_Interacao', [], [], ['drugbank_id1','drugbank_id2'], 
                                            ['drugbank_id1','drugbank_id2'], ['DrugBank_Nome','DrugBank_Nome'], ['drugbank_id','drugbank_id']) # FK | Referenced Tables & Columns
# ['nomeProduto', 'id']
sqlAnvisa_Name = sql.SQLScripting(df_Anvisa_Names, 'Anvisa_Nome', ['id','nomeProduto'], [], ['id'])
# ['nome_pAtivo', 'id_pAtivo', 'translated_pAtivo']
sqlAnvisa_Principles = sql.SQLScripting(df_Anvisa_PrinciplesAccented, 'Anvisa_PrincipioAtivo', ['id_pAtivo', 'nome_pAtivo', 'translated_pAtivo'], [], ['id_pAtivo'],)
# ['idProduto', 'idPrincipio']
sqlAvisa_NameActPrinciple = sql.SQLScripting(df_Anvisa_Names_Principles, 'Anvisa_Nome_PrincipioAtivo', [], [], ['idProduto', 'idPrincipio'],
                                             ['idProduto', 'idPrincipio'], ['Anvisa_Nome', 'Anvisa_PrincipioAtivo'], ['id', 'id_pAtivo'])
# ['drugbank-id', 'name_drugbank', 'id_pAtivo', 'name_anvisa','matchingValue']
df_DrugBank_Anvisa.rename(columns={'drugbank-id': 'drugbank_id'}, inplace = True)
sqlDrugBank_Anvisa = sql.SQLScripting(df_DrugBank_Anvisa, 'DrugBank_Anvisa', ['drugbank_id','id_pAtivo','matchingValue'], [], ['drugbank_id','id_pAtivo'],
                                      ['drugbank_id', 'id_pAtivo'], ['DrugBank_Nome', 'Anvisa_PrincipioAtivo'], ['drugbank_id', 'id_pAtivo'])
# ['id_principal', 'nome', 'tipo_origem']
sqlBigTable = sql.SQLScripting(df_BigTable, 'BigTable_Nomes', [], [], [])


# KEGG Drug Specific SQL Scripts
# ['keggdrug-id', 'name']
if not silent: print('SQL Scripts - Creating CREATE and INSERT scripts for KEGG Drug database tables')
df_KEGG_drugs.rename(columns={'keggdrug-id': 'keggdrug_id'}, inplace = True)
sqlKEGGDrug_Name = sql.SQLScripting(df_KEGG_drugs,'KEGGDrug_Nome', [], [], ['keggdrug_id'])
# ['keggdrug-id1', 'keggdrug-id2']
df_KEGG_interaction .rename(columns={'keggdrug-id1': 'keggdrug_id1', 'keggdrug-id2': 'keggdrug_id2'}, inplace = True)
sqlKEGGDrug_Interactions = sql.SQLScripting(df_KEGG_interaction, 'KEGGDrug_Interacao', [], [], ['keggdrug_id1','keggdrug_id2'], 
                                            ['keggdrug_id1','keggdrug_id2'], ['KEGGDrug_Nome','KEGGDrug_Nome'], ['keggdrug_id','keggdrug_id']) # FK | Referenced Tables & Columns
# ['keggdrug-id', 'id_pAtivo', 'matchingValue']
df_KEGGDrug_Anvisa.rename(columns={'keggdrug-id': 'keggdrug_id'}, inplace = True)
sqlKEGGDrug_Anvisa = sql.SQLScripting(df_KEGGDrug_Anvisa, 'KEGGDrug_Anvisa', ['keggdrug_id','id_pAtivo','matchingValue'], [], ['keggdrug_id','id_pAtivo'],
                                      ['keggdrug_id', 'id_pAtivo'], ['KEGGDrug_Nome', 'Anvisa_PrincipioAtivo'], ['keggdrug_id', 'id_pAtivo'])


# Exporting Scripts (Scripts Directory)
if not silent: print('SQL Scripts - Exporting script files')
sqlDrugBank_Name.exportSQLScripts()
sqlDrugBank_Interactions.exportSQLScripts()
sqlAnvisa_Name.exportSQLScripts()
sqlAnvisa_Principles.exportSQLScripts()
sqlAvisa_NameActPrinciple.exportSQLScripts()
sqlDrugBank_Anvisa.exportSQLScripts()
sqlBigTable.exportSQLScripts()
# KEGG
sqlKEGGDrug_Name.exportSQLScripts()
sqlKEGGDrug_Interactions.exportSQLScripts()
sqlKEGGDrug_Anvisa.exportSQLScripts()

timeTracker.note(strSubject,'end')

timeTracker.export()
print('ETL Successfully Executed!')