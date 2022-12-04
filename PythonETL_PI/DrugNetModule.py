import pandas as pd
import igraph as ig
import networkx as nx

# In
exp_csv_interactions = r"..\Exported\exp_csv_interactions.csv"
exp_csv_KEGG_Interaction = r"..\Exported\exp_csv_KEGG_Interaction.csv"

# Out
exp_graph_DrugBank = r"..\Exported\Analysis\exp_graph_DrugBank.net"
exp_graph_KEGGDrug = r"..\Exported\Analysis\exp_graph_KEGGDrug.net"
exp_graph_DrugBank_dot = r"..\Exported\Analysis\exp_graph_DrugBank.dot"
exp_graph_KEGGDrug_dot = r"..\Exported\Analysis\exp_graph_KEGGDrug.dot"

# Load DataFrames
df_DrugBank_edges = pd.read_csv(exp_csv_interactions, sep=',')
df_KEGGDrug_edges = pd.read_csv(exp_csv_KEGG_Interaction, sep=',')

graph_DrugBank = ig.Graph.DataFrame(df_DrugBank_edges)
graph_DrugBank.write_pajek(exp_graph_DrugBank)
graph_DrugBank.write_dot(exp_graph_DrugBank_dot)

graph_KEGGDrug = ig.Graph.DataFrame(df_KEGGDrug_edges)
graph_KEGGDrug.write_pajek(exp_graph_KEGGDrug)
graph_KEGGDrug.write_dot(exp_graph_KEGGDrug_dot)

#ig.plot(graph_DrugBank)
