import cobra
from cobra.sampling import sampling
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
import pandas as pd
import os
from os.path import join
from cobra import Model, Reaction, Metabolite
import numpy as np
os.environ["R_HOME"] = f"{os.environ['CONDA_PREFIX']}\\Lib\\R"
import rpy2.robjects
from plotnine import *

# Importing the model
model=cobra.io.read_sbml_model("C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Model_files\\Rt_IFO0880.xml")

# Description of the model
print(f'{len(model.reactions)} reactions initially')
print('%i metabolites initially' % len(model.metabolites))
print('%i genes initially' % len(model.genes))

#Print the reaction of biomasss
solution = model.optimize() #his will maximize or minimize (maximizing is the default) flux through the objective reactions 
print(solution) # gives the flux for the objective function
# print(model.summary()) # Displays information on the input and output behavior of the model, along with the optimized objective

# print(model.metabolites.nadh_c.summary()) # For inspecting the input-output behavior of individual metabolites
# print(model.metabolites.atp_c.summary()) # For seeing the main energy production and consumption reactions

# Change glucose uptake
medium = model.medium 
print(medium) # shows current growth medium: returns a dictionary that contains the upper flux bounds for all active exchange fluxes (the ones having non-zero flux bounds
medium['EX_glc__D_e'] = 2.0 # modifying growth medium to their respective upper import bounds
model.medium = medium

# print(model.medium)
print(model.optimize())

# Change upper and lower bounds
# model.reactions.EX_glc__D_e.bounds = -1000, 1000

print(model.objective.expression) # shows the objectives mathematical exspression

# Changing the objective function, ex to ATPM
model.objective = "ATPM"  # change objective function to ATP minimization, set bounds for biomass (lower and upper bound the same) max upper 0.3
#print all fluxes, find biomass function, and NADPH functions (in PPP)
print(model.objective.expression) 
print(model.optimize('minimize'))

# The upper bound should be 1000, so that we get
# the actual optimal value
# model.reactions.get_by_id("ATPM").upper_bound = 1000
# linear_reaction_coefficients(model)

# from pathlib import Path
# from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
# import logging
# save_json_model(model, "edited_model.json")

