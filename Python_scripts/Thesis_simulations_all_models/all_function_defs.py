# Functions used for model simulations and visualization of results 
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
from matplotlib.pyplot import figure
from cobra.flux_analysis.loopless import loopless_solution
from matplotlib import colormaps
import matplotlib
from cobra.sampling import sample


def all_fluxes_biomass_max_df(model_path: str, glucose_uptakes: list, biomass_rxn_ID: str, glc_ID: str):
    # Creating an empty dataframe, which has enzymes as columns and as many rows as glucose uptakes
    model1=cobra.io.read_sbml_model(model_path)
    solution1 = model1.optimize()
    metabolites = solution1.fluxes.to_frame(name='Flux')
    all_fluxes_biomass_max = pd.DataFrame(columns=[*metabolites.index], index=range(len(glucose_uptakes))) #flux_values.index gives the row names column, * extracts the list of strings
 
    # Calculating flux data on different glucose uptakes
    for i in range(len(glucose_uptakes)):
       
        model=cobra.io.read_sbml_model(model_path) # It's important to load the original model again each time before optimizing again (otherwise the solution differs as the initial conditions would differ)
        model.objective = biomass_rxn_ID 
        
        model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  # NB! It's also possible to change glucose uptake from model.medium dictionary, but that gives differences in solution fluxes
         
        model.optimize()
        solution = loopless_solution(model)
        
        all_fluxes_biomass_max.loc[i] = solution.fluxes[[*metabolites.index]].values

    return all_fluxes_biomass_max 



def all_fluxes_NGAM_min_df(model_path: str, glucose_uptakes: list, growth_rates: list, NGAM_rxn_ID: str, glc_ID: str, biomass_rxn_ID: str):
    # Creating an empty dataframe, which has enzymes as columns and as many rows as glucose uptakes
    model1=cobra.io.read_sbml_model(model_path)
    solution1 = model1.optimize()
    enzymes = solution1.fluxes.to_frame(name='Flux')
    all_fluxes_NGAM_min = pd.DataFrame(columns=[*enzymes.index], index=range(len(glucose_uptakes))) #flux_values.index gives the row names column, * extracts the list of strings
 
    # Calculating flux data with all glucose uptakes
    for i in range(len(glucose_uptakes)):
       
        model=cobra.io.read_sbml_model(model_path)
        model.objective = NGAM_rxn_ID 
        model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  
        model.reactions.get_by_id(biomass_rxn_ID).bounds = growth_rates[i], growth_rates[i]

         
        model.optimize('minimize') #
        solution = loopless_solution(model)
        
        all_fluxes_NGAM_min.loc[i] = solution.fluxes[[*enzymes.index]].values

    return all_fluxes_NGAM_min 



def metabolites_fluxes(model_path: str, all_fluxes_df, metabolites: list):
    
    model = cobra.io.read_sbml_model(model_path)
    metabolites_fluxes = all_fluxes_df[metabolites]
    
    # Sum two phosphoketolases in Rt_IFO0880 models together for easier comparison with models that have only one phosphoketolase
    if 'XPK' in metabolites:
        metabolites_fluxes.loc[:, 'XPK'] = metabolites_fluxes.loc[:, 'XPK'] + all_fluxes_df.loc[:, 'FPK']

    # Change metabolites IDs to their name for easier understanding
    
    # Because in Rt_IFO0880 models TKT1 and TKT2 have the same name (Transketolase), this models namehave to be changed differently
    if 'TKT1' in metabolites:
        for i in range(len(metabolites)):
            if metabolites_fluxes.columns[i] != 'TKT1' and metabolites_fluxes.columns[i] != 'TKT2':
                metabolites_fluxes = metabolites_fluxes.rename(columns = {metabolites_fluxes.columns[i]: getattr(model.reactions, metabolites_fluxes.columns[i]).name})
    else:
        for i in range(len(metabolites)):
            metabolites_fluxes = metabolites_fluxes.rename(columns = {metabolites_fluxes.columns[i]: getattr(model.reactions, metabolites_fluxes.columns[i]).name})

    # Change some names 
    metabolites_fluxes = metabolites_fluxes.rename(columns = {'Glucose 6-phosphate dehydrogenase': 'oxPPP', 'glucose 6-phosphate dehydrogenase': 'oxPPP', 'non-growth associated maintenance reaction': 'NGAM',
                                                              'ATP maintenance requirement': 'NGAM', 'Xylulose-5-phosphate phosphoketolase': 'Phosphoketolase', 
                                                              'TKT1': 'Transketolase 1', 'TKT2': 'Transketolase 2', 'phosphoketolase (fructose 6-phosphate)': 'phosphoketolase'})

    return metabolites_fluxes



def plot_ex_intr_fluxes(all_fluxes_df, exchange_fluxes, intracellular_fluxes, ACL_phosphoketolase, title: str, biomass_rxn_ID: str): 
    # Plot exchange and intracellular fluxes, add lab data to exchange fluxes plot
    
    # Lab data
    GR = [0.049, 0.100, 0.151, 0.203, 0.301]
    glc = [0.476, 1.114, 1.648, 2.305, 3.1]
    co2 = [1.171, 2.521, 3.854, 5.834, 7.415]
    o2 = [1.083, 2.521, 3.851, 4.352, 6.327]

    data = {'GR': GR, 'glc': glc, 'o2': o2, 'co2': co2}
    lab_data = pd.DataFrame(data)
    
    
    fig, ax = plt.subplots(1, 2, figsize=(12,6)) #, 
    fig.suptitle(title)

    # Exp data
    x1 = lab_data['GR']
    y1 = lab_data[['glc', 'o2', 'co2']]
    # Simulations data
    x = all_fluxes_df[biomass_rxn_ID]
    y = np.abs(exchange_fluxes)
    
    ax[0].plot(x1, y1, '--', label= ['Exp glucose exchange', 'Exp oxygen exchange', 'Exp carbon dioxide exchange'], linewidth=0.85) #
    ax[0].plot(x, y, '-', label= y.columns) #
    
    ax[0].set_xlim([0.024, 0.24])  
    ax[0].legend(fontsize=10, loc='upper left')
    ax[0].set_title("(a) Exchange fluxes", fontsize = 13) # fluxes biomass maximization
    ax[0].set_xlabel('Biomass growth rate $(1/h)$')
    ax[0].set_ylabel('Flux $(mmol/gDW/h)$')

    y3 = np.abs(intracellular_fluxes)
    y4 = np.abs(ACL_phosphoketolase)

    ax[1].plot(x, y3, '--', label= y3.columns, linewidth=0.85) #
    ax[1].plot(x, y4, '-', label= y4.columns, linewidth=2) # ACL and phosphoketolase

    ax[1].legend(fontsize=10, loc='upper left')
    ax[1].set_title("(b) Intracellular fluxes", fontsize = 13)
    ax[1].set_xlabel('Biomass growth rate $(1/h)$')
    ax[1].set_ylabel('Flux $(mmol/gDW/h)$')
    
    return fig    
    
def cofactor_balances_biomass_max(model_path: str, cofactor_list: list, glucose_uptakes: list, i: int, biomass_rxn_ID: str, glc_ID: str):
    producing_fluxes = pd.DataFrame() 
    consuming_fluxes = pd.DataFrame()
    
    model=cobra.io.read_sbml_model(model_path)
    model.objective = biomass_rxn_ID 
    model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  
        
    model.optimize()
    solution = loopless_solution(model)
        
    for metabolite in cofactor_list:
        producing_fluxes = pd.concat([producing_fluxes, getattr(model.metabolites, metabolite).summary().producing_flux])
        consuming_fluxes = pd.concat([consuming_fluxes, getattr(model.metabolites, metabolite).summary().consuming_flux])
    
    for reaction in producing_fluxes.index:
        if reaction in consuming_fluxes.index and abs(round(producing_fluxes.loc[reaction, 'flux'], 3)) == abs(round(consuming_fluxes.loc[reaction, 'flux'], 3)):
            producing_fluxes = producing_fluxes.drop([reaction])
            consuming_fluxes = consuming_fluxes.drop([reaction])
            
    cofactor_fluxes = pd.concat([producing_fluxes, consuming_fluxes])

    cofactor_fluxes = cofactor_fluxes.sort_values(by='flux', ascending=False).drop(columns = ['percent']) # drop percent column, bc these percents are not for nadph sum (the percent is for specific compartment)
    cofactor_fluxes = cofactor_fluxes[(cofactor_fluxes['flux']) != 0.0] # for getting non-zero fluxes only

    cofactor_sum_producing_flux = sum(cofactor_fluxes[cofactor_fluxes['flux'] > 0]['flux']) # for getting the sum of producing fluxes
    cofactor_sum_consuming_flux = sum(cofactor_fluxes[cofactor_fluxes['flux'] < 0]['flux']) # for getting the sum of consumed fluxes
    print(f'SUM produced: {cofactor_sum_producing_flux}, SUM consumed: {cofactor_sum_consuming_flux}')

    if round(cofactor_sum_producing_flux, 3) == round(abs(cofactor_sum_consuming_flux), 3):
        cofactor_fluxes['percent'] = abs(cofactor_fluxes['flux']/cofactor_sum_producing_flux) # add percent column

    return cofactor_fluxes 



def cofactor_balances_NGAM_min(model_path: str, cofactor_list: list, glucose_uptakes: list, growth_rates: list, NGAM_rxn_ID: str, glc_ID: str, biomass_rxn_ID: str, i: int):
    producing_fluxes = pd.DataFrame() 
    consuming_fluxes = pd.DataFrame()
    
    model=cobra.io.read_sbml_model(model_path)
    
    model.objective = NGAM_rxn_ID 
    model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  
    model.reactions.get_by_id(biomass_rxn_ID).bounds = growth_rates[i], growth_rates[i]

        
    model.optimize('minimize') #
    solution = loopless_solution(model)
        
    for metabolite in cofactor_list:
        producing_fluxes = pd.concat([producing_fluxes, getattr(model.metabolites, metabolite).summary().producing_flux])
        consuming_fluxes = pd.concat([consuming_fluxes, getattr(model.metabolites, metabolite).summary().consuming_flux])
    
    for reaction in producing_fluxes.index:
        if reaction in consuming_fluxes.index and abs(round(producing_fluxes.loc[reaction, 'flux'], 3)) == abs(round(consuming_fluxes.loc[reaction, 'flux'], 3)):
            producing_fluxes = producing_fluxes.drop([reaction])
            consuming_fluxes = consuming_fluxes.drop([reaction])
            
    cofactor_fluxes = pd.concat([producing_fluxes, consuming_fluxes])

    cofactor_fluxes = cofactor_fluxes.sort_values(by='flux', ascending=False).drop(columns = ['percent']) # drop percent column, bc these percents are not for nadph sum (the percent is for specific compartment)
    cofactor_fluxes = cofactor_fluxes[(cofactor_fluxes['flux']) != 0.0] # for getting non-zero fluxes only

    cofactor_sum_producing_flux = sum(cofactor_fluxes[cofactor_fluxes['flux'] > 0]['flux']) # for getting the sum of producing fluxes
    cofactor_sum_consuming_flux = sum(cofactor_fluxes[cofactor_fluxes['flux'] < 0]['flux']) # for getting the sum of consumed fluxes
    print(f'SUM produced: {cofactor_sum_producing_flux}, SUM consumed: {cofactor_sum_consuming_flux}')

    if round(cofactor_sum_producing_flux, 3) == round(abs(cofactor_sum_consuming_flux), 3):
        cofactor_fluxes['percent'] = abs(cofactor_fluxes['flux']/cofactor_sum_producing_flux) # add percent column

    return cofactor_fluxes 




def cofactor_fluxes_pie_chart(model_path: str, cofactor_fluxes, **fig_kw):
    threshold = 0.025 # threshold shows the percent of the flux for including in others sector on pie chart 
    # The three lines below are for grouping together reactions with low fluxes in producing
    producing_cofactor_fluxes_draw = cofactor_fluxes[(cofactor_fluxes['flux'] > 0).copy()]    
    producing_cofactor_fluxes_draw.loc[producing_cofactor_fluxes_draw['percent'] < threshold, 'reaction'] = 'Other producing'
    producing_cofactor_fluxes_draw = producing_cofactor_fluxes_draw.groupby('reaction')[['percent', 'flux']].sum()        
    
    # The three lines below are for grouping together reactions with low fluxes in consuming 
    consuming_cofactor_fluxes_draw = cofactor_fluxes[(cofactor_fluxes['flux'] < 0).copy()]
    consuming_cofactor_fluxes_draw.loc[consuming_cofactor_fluxes_draw['percent'] < threshold, 'reaction'] = 'Other consuming'
    consuming_cofactor_fluxes_draw = consuming_cofactor_fluxes_draw.groupby('reaction')[['percent', 'flux']].sum()
    
    # y_producing = producing_cofactor_fluxes_draw['percent']
    # labels_producing = producing_cofactor_fluxes_draw[['reaction', 'flux']] 
    
    # y_consuming = abs(consuming_cofactor_fluxes_draw['percent'])
    # labels_consuming = consuming_cofactor_fluxes_draw[['reaction', 'flux']]
    
    producing_and_consuming_fluxes = pd.concat([producing_cofactor_fluxes_draw, consuming_cofactor_fluxes_draw])
    
    # reaction_IDs = pd.concat([labels_producing, labels_consuming])
    
    model = cobra.io.read_sbml_model(model_path)
    solution = model.optimize()   
    
    reaction_names_w_flux = []
    for reaction in producing_and_consuming_fluxes.index:
        if reaction != 'Other producing' and reaction != 'Other consuming':
            reaction_names_w_flux += [''.join([''.join([str(round((producing_and_consuming_fluxes.loc[reaction, 'percent'])*100, 1)),'% ']), getattr(model.reactions, reaction).name, ' (', str(round(producing_and_consuming_fluxes.loc[reaction, 'flux'], 2)), ')'])]
        elif reaction == 'Other producing':
            reaction_names_w_flux += [''.join([''.join([str(round((producing_and_consuming_fluxes.loc[reaction, 'percent'])*100, 1)),'% ']), 'Other producing', ' (', str(round(producing_and_consuming_fluxes.loc[reaction, 'flux'], 2)), ')'])] 
        elif reaction == 'Other consuming':
            reaction_names_w_flux += [''.join([''.join([str(round((producing_and_consuming_fluxes.loc[reaction, 'percent'])*100, 1)),'% ']), 'Other consuming', ' (', str(round(producing_and_consuming_fluxes.loc[reaction, 'flux'], 2)), ')'])] 
    fig = plt.figure()
         
    pie_chart = plt.pie(producing_and_consuming_fluxes.loc[:, 'percent'], labels = reaction_names_w_flux)  #autopct='%1.1f%%' pd.concat([producing_cofactor_fluxes_draw, consuming_cofactor_fluxes_draw])[['flux', 'percent']]
    # plt.tight_layout()
    
    fig.set_size_inches(15, 6)


    # plt.legend(producing_and_consuming_fluxes, reaction_names_w_flux, title = 'Reaction names', loc="center left",  bbox_to_anchor=(1, 0, 0.5, 1))
    
    # plt.title(title)

    return pie_chart, fig


# Get all fluxes to excel
def all_fluxes_to_excel(path: str, all_fluxes_df):
    with pd.ExcelWriter(path) as excel_writer:
        all_fluxes_df.to_excel(excel_writer, sheet_name='Sheet 1', index=True)
        
# Get csv file of fluxes (from 0-4)

def fluxes_to_csv(path: str, all_fluxes_df, i: int):
    all_fluxes_df.loc[i].to_csv(path, index=True)
    
    
def flux_sampling():
    s = sample(model, 100)
    s.DESCR()
