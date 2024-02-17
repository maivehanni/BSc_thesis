# Set objective function to ATP minimization
# Set bounds for biomass: 0.05, 0.10, 0.15, 0.20, 0.25, 0.30 (max Rhodotorula growth rate) - if there're problems, set lower and upper bound as the same
# Print all fluxes, find biomass function and NADPH functions (in PPP)
# 'TALA'=Transaldolase, 'TKT1'=Transketolase, 'TKT2'=Transketolase, 'RPI'=Ribose-5-phosphate isomerase, 
# 'G6PDH2rp'=Glucose 6-phosphate dehydrogenase, 'GND'=Phosphogluconate dehydrogenase, 'PGLp'=6-phosphogluconolactonase


# Make a graph where specific growth rate is on x-axis and the flux is on y-axis

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
import matplotlib.pyplot as plt

# Importing the model
model=cobra.io.read_sbml_model("C:\\Users\\Maive\\Desktop\\BSc_loputoo\\Model_files\\Rt_IFO0880.xml")

# Changing the objective function
#model.objective = "ATPM" 

solution = model.optimize() # use model.optimize('minimize') to minimize (maximizing is the default) flux through the objective reactions 

# Analyzing FBA solutions: model.summary() - displays information on the input and output behaviour along with the objective function (but not about the intracellular reactions)
print(model.summary())

# The input-output behavior of individual metabolites can also be inspected using summary methods.
# print(solution.fluxes.r_0466)
# print(solution.fluxes.r_0659)
# print(solution.fluxes.r_0889)
# print(solution.fluxes.r_2140)
# print(solution.fluxes.r_2141)



# Show all fluxes for objevtive funstion
all_fluxes = solution.fluxes

# Change bounds for biomass: 0.05, 0.10, 0.15, 0.20, 0.25, 0.30
# Different ways of changing the bounds
# biomass_reaction = model.reactions.get_by_id('BIOMASS_RT')
# biomass_reaction.upper_bound = 0.05
# biomass_reaction.lower_bound = 0.05
# biomass_reaction_upper_bound = model.reactions.get_by_id("BIOMASS_RT").upper_bound = 0.05
# biomass_reaction_bounds = model.reactions.BIOMASS_RT.bounds = 0.05, 0.05
biomass_reaction = model.reactions.get_by_id('BIOMASS_RT')
biomass_reaction.upper_bound = 0.05
biomass_reaction.lower_bound = 0.05


# Extract specific flux values 
ppp_enzymes = ['TALA', 'TKT1', 'TKT2', 'RPI', 'G6PDH2rp', 'GND', 'PGLp']

flux_values_of_interest = []
for i in ppp_enzymes:
    flux_values_of_interest += [solution.fluxes[i]]

# print(flux_values_of_interest)

# Plot the flux values

# x = range(1, 8)
# y = flux_values_of_interest[:]

# plt.scatter(x, y)
# plt.xlabel('Reactions')
# plt.ylabel('Flux Values')
# plt.title('Flux Values for Metabolites of Interest')
# plt.show()


