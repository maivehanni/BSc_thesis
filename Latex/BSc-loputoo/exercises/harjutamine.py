# Simulating deletions

import pandas
from time import time

from cobra.io import load_model
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

cobra_model = load_model("textbook")
ecoli_model = load_model("iJO1366")

# Knocking out single genes and reactions
print('complete model: ', cobra_model.optimize())
# with cobra_model:
#     cobra_model.reactions.PFK.knock_out()
#     print('pfk knocked out: ', cobra_model.optimize())
# print('complete model: ', cobra_model.optimize())

## Knocked out genes can also have no effect in case of redundancy, or some gene knock outs can affect more reactions
# with cobra_model:
#     cobra_model.genes.b1723.knock_out()
#     print('pfkA knocked out: ', cobra_model.optimize())
#     cobra_model.genes.b3916.knock_out()
#     print('pfkB knocked out: ', cobra_model.optimize()) 

# Performing all single gene deletions on a model or for a subset of genes
#all_single_gene_deletions_results = single_gene_deletion(cobra_model)
subset_deletions = single_gene_deletion(cobra_model, cobra_model.genes[:20])
print(subset_deletions)

#This can also be done for reactions
single_reaction_deletion(cobra_model, cobra_model.reactions[:20])
