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
model.objective = "BIOMASS_RT" 

# Lab data
glucose_uptakes= [0.476, 1.114, 1.648, 2.305, 2.6619851, 3.1] # the fifth glc uptake value was calculated


# Changing glucose uptake
medium = model.medium
medium["EX_glc__D_e"] = glucose_uptakes[0]
model.medium = medium

# Optimize the model
solution = model.optimize()   

all_fluxes = solution.fluxes.to_frame(name='Flux')
ACITL_flux = all_fluxes.loc['BIOMASS_RT'] # example for getting specific flux
print(ACITL_flux)
# model.summary()

# # Get all fluxes to excel
# with pd.ExcelWriter('C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Results\\Simulated_fluxes\\Biomass_maximization\\all_fluxes_glc_uptake_005.xlsx') as excel_writer:
#     all_fluxes.to_excel(excel_writer, sheet_name='Glucose uptake = 0.476', index=True)

# #Make a csv file of all fluxes
# # Get all flux values separately for dif growth rates, make them to a csv file
# all_fluxes.to_csv(f'C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Results\\Simulated_fluxes\\Biomass_maximization\\all_fluxes_glc_uptake_05.csv', index=True)
    
# # Get json file
# from pathlib import Path
# from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
# import logging

# save_json_model(model, f"C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Model_files\\edited_Rt_IFO0880_model_glc_uptake_005.json")