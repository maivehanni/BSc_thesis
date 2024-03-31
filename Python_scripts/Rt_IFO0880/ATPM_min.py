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
import copy

model_data1=cobra.io.read_sbml_model("C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Model_files\\Rt_IFO0880.xml")
# # FLuxes when optimized for ATPM
model_data1.objective = "ATPM" 
solution1 = model_data1.optimize('minimize')

model_data=cobra.io.read_sbml_model("C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Model_files\\Rt_IFO0880.xml")

# Lab data
growth_rates = [0.03, 0.08, 0.12, 0.17, 0.20, 0.23] #[0.049, 0.100, 0.151, 0.203, 0.25, 0.301]
glucose_uptakes= [-0.476, -1.114, -1.648, -2.305, -2.7, -3.1] # the fifth glc uptake value was calculated

# All fluxes
all_fluxes = solution1.fluxes.to_frame(name='Flux')

all_fluxes_dif_GR_ATPM_min = pd.DataFrame(columns=['Growth rate', *all_fluxes.index], index=range(len(growth_rates))) #flux_values.index gives the row names column, * extracts the list of strings

for i in range(len(growth_rates)):

    model = copy.deepcopy(model_data)
    model.objective = "ATPM"
    model.reactions.BIOMASS_RT.bounds = growth_rates[i], growth_rates[i]

    medium = model.medium
    medium["EX_glc__D_e"] = -(glucose_uptakes[i])
    model.medium = medium
    solution = model.optimize('minimize')
    
    all_fluxes_dif_GR_ATPM_min.loc[i] = solution.fluxes[['BIOMASS_RT', *all_fluxes.index]].values

print(all_fluxes_dif_GR_ATPM_min[['BIOMASS_RT', 'EX_glc__D_e', 'G6PDH2r', 'XPK', 'ATPM', 'ACITL']])
