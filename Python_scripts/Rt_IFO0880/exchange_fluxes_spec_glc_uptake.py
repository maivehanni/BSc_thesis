import cobra
import pandas as pd
import os
from os.path import join
from cobra import Model, Reaction, Metabolite
from cobra.sampling import sampling
import numpy as np
os.environ["R_HOME"] = f"{os.environ['CONDA_PREFIX']}\\Lib\\R"
import rpy2.robjects
from plotnine import *
import matplotlib.pyplot as plt

# Importing the model
model=cobra.io.read_sbml_model("C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Model_files\\Rt_IFO0880.xml")

# Change glucose bounds
model.reactions.get_by_id("EX_glc__D_e").upper_bound = 9999
model.reactions.get_by_id("EX_glc__D_e").lower_bound = -9999

# Changing the objective function to glucose min
model.objective = "BIOMASS_RT" 
solution = model.optimize() 

# Glucose uptakes, lab data
glucose_uptakes= [-0.476, -1.114, -1.648, -2.305, -2.6619851, -3.1]

# Get all fluxes on different glucose uptake
all_fluxes = solution.fluxes.to_frame(name='Flux')

flux_values_specific_glucose_uptake = pd.DataFrame(columns=['Glucose uptake', *all_fluxes.index], index=range(len(glucose_uptakes))) #flux_values.index gives the row names column, * extracts the list of strings
biomass_GR = []

for i in range(len(glucose_uptakes)):
    model.reactions.EX_glc__D_e.bounds = glucose_uptakes[i], glucose_uptakes[i]
    solution = model.optimize('maximize')
    flux_values_specific_glucose_uptake.loc[i] = solution.fluxes[['EX_glc__D_e', *all_fluxes.index]].values
    biomass_GR += [solution.objective_value]
flux_values_specific_glucose_uptake.insert(0, 'Biomass growth rate', biomass_GR, True)

# Example for finding a certain reaction flux: flux_values_specific_glucose_uptake['EX_o2_e']

# Adding new column on specified location: df.insert(column_nr, "column_name", data=[21, 23, 24, 21], allow_duplicates=True)

print(flux_values_specific_glucose_uptake)

# # Get all fluxes to excel
# with pd.ExcelWriter('C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Results\\Simulated_fluxes\\flux_values_specific_glucose_uptake_all.xlsx') as excel_writer:
#     flux_values_specific_glucose_uptake.to_excel(excel_writer, sheet_name='Sheet1', index=False)

# #  Get all flux values separately for dif growth rates, make them to a csv file
# for i in range(len(glucose_uptakes)):
#     flux_values_specific_glucose_uptake.loc[i].t

# Getting exchange fluxes: glucose, CO2, ammonium and other minerals, others, like glycerol ['glc__D_e', 'o2_e', 'glyc_e', 'nh4_e',	'so4_e',	'pi_e', 'co2_e']  

# Make a pd dataframe with all exchange fluxes that are not zero, then make a pivot table with wanted metabolites fluxes on different growth rates
# biomass_GR = []

for i in range(1, len(glucose_uptakes)):
    model.reactions.EX_glc__D_e.bounds = glucose_uptakes[i], glucose_uptakes[i]
    solution = model.optimize()
    model_summary = model.summary().to_frame()
    # biomass_GR += [solution.objective_value]

exchange_fluxes_all = model_summary
exchange_fluxes_all.insert(0, 'Biomass growth rate', biomass_GR, True)
    
exchange_fluxes_all = exchange_fluxes_all[(exchange_fluxes_all['flux']) != 0.0] # for getting non-zero fluxes only
exchange_fluxes_all['flux'] = abs(exchange_fluxes_all['flux'])

# Get all non-zero exchange fluxes to a pivot table
exchange_fluxes_table_all = pd.pivot_table(exchange_fluxes_all, values='flux', index=['GR'], columns=['metabolite'])

# Get specific metabolites with their fluxes
exchange_fluxes_table = pd.pivot_table(exchange_fluxes_all[exchange_fluxes_all.metabolite.isin(['glc__D_e', 'o2_e', 'glyc_e', 'nh4_e',	'so4_e',	'pi_e', 'co2_e'])], 
                                       values='flux', index=['GR'], columns=['metabolite'])

# Make plots for exchange reactions

plt.plot(biomass_GR, exchange_fluxes_table[exchange_fluxes_table.columns], 'o', label = exchange_fluxes_table.columns)

plt.xlabel('Biomass growth rate')
plt.ylabel('Flux')
plt.title("Exchange fluxes")
plt.legend(fontsize=6, loc='upper left', bbox_to_anchor=(1.05, 1))
plt.yticks(range(0, 20, 1))
plt.show()

