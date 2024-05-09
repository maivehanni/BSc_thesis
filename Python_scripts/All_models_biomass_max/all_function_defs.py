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
from cobra.flux_analysis.loopless import loopless_solution

def all_fluxes_biomass_max_df(model_path, glucose_uptakes, biomass_rxn_ID, glc_ID):
    # Creating an empty dataframe, which has enzymes as columns and as many rows as glucose uptakes
    model1=cobra.io.read_sbml_model(model_path)
    solution1 = model1.optimize()
    enzymes = solution1.fluxes.to_frame(name='Flux')
    all_fluxes_biomass_max = pd.DataFrame(columns=[*enzymes.index], index=range(len(glucose_uptakes))) #flux_values.index gives the row names column, * extracts the list of strings
 
    # Calculating flux data with all glucose uptakes
    for i in range(len(glucose_uptakes)):
       
        model=cobra.io.read_sbml_model(model_path)
        model.objective = biomass_rxn_ID 
        
        model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  
         
        model.optimize()
        solution = loopless_solution(model)
        
        all_fluxes_biomass_max.loc[i] = solution.fluxes[[*enzymes.index]].values

    return all_fluxes_biomass_max 



def all_fluxes_NGAM_min_df(model_path, glucose_uptakes, growth_rates, obj_func_ID, glc_ID, biomass_rxn_ID):
    # Creating an empty dataframe, which has enzymes as columns and as many rows as glucose uptakes
    model1=cobra.io.read_sbml_model(model_path)
    solution1 = model1.optimize()
    enzymes = solution1.fluxes.to_frame(name='Flux')
    all_fluxes_NGAM_min = pd.DataFrame(columns=[*enzymes.index], index=range(len(glucose_uptakes))) #flux_values.index gives the row names column, * extracts the list of strings
 
    # Calculating flux data with all glucose uptakes
    for i in range(len(glucose_uptakes)):
       
        model=cobra.io.read_sbml_model(model_path)
        model.objective = obj_func_ID 
        model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  
        model.reactions.get_by_id(biomass_rxn_ID).bounds = growth_rates[i], growth_rates[i]

         
        model.optimize('minimize')
        solution = loopless_solution(model)
        
        all_fluxes_NGAM_min.loc[i] = solution.fluxes[[*enzymes.index]].values

    return all_fluxes_NGAM_min 



def metabolites_fluxes(model_path, all_fluxes_df, metabolites):
    
    model = cobra.io.read_sbml_model(model_path)
    metabolites_fluxes = all_fluxes_df[metabolites]
    
    # Sum two phosphoketolases together for easier comparison with models that have only one phosphoketolase
    if 'XPK' in metabolites:
        metabolites_fluxes.loc[:, 'XPK'] = metabolites_fluxes.loc[:, 'XPK'] + all_fluxes_df.loc[:, 'FPK']

    if 'TKT1' in metabolites:
        # Change metabolites IDs with their name for easier understanding
        for i in range(len(metabolites)):
            if metabolites_fluxes.columns[i] != 'TKT1' and metabolites_fluxes.columns[i] != 'TKT2':
                metabolites_fluxes = metabolites_fluxes.rename(columns = {metabolites_fluxes.columns[i]: getattr(model.reactions, metabolites_fluxes.columns[i]).name})
    else:
        for i in range(len(metabolites)):
            metabolites_fluxes = metabolites_fluxes.rename(columns = {metabolites_fluxes.columns[i]: getattr(model.reactions, metabolites_fluxes.columns[i]).name})

    # Change some names 
    metabolites_fluxes = metabolites_fluxes.rename(columns = {'Glucose 6-phosphate dehydrogenase': 'oxPPP', 'ATP maintenance requirement': 'NGAM', 'Xylulose-5-phosphate phosphoketolase': 'Phosphoketolase', 
                                                              'TKT1': 'Transketolase 1', 'TKT2': 'Transketolase 2'})

    return metabolites_fluxes



def plot_ex_intr_fluxes(all_fluxes_df, exchange_fluxes, intracellular_fluxes, ACL_phosphoketolase, title, biomass_rxn_ID):
    # Plot exchange and intracellular fluxes
    fig, ax = plt.subplots(1, 2, figsize=(12,6)) #, 
    fig.suptitle(title)

    # Sample data

    x1 = all_fluxes_df[biomass_rxn_ID]
    y1 = np.abs(exchange_fluxes)

    ax[0].plot(x1, y1, '-', label= y1.columns) #
    ax[0].legend(fontsize=10, loc='upper left')
    ax[0].set_title("Exchange fluxes") #fluxes biomass maximization
    ax[0].set_xlabel('Biomass growth rate $(1/h)$')
    ax[0].set_ylabel('Flux $(mmol/gDW/h)$')

    x2 = all_fluxes_df[biomass_rxn_ID]
    y2 = np.abs(intracellular_fluxes)
    y3 = np.abs(ACL_phosphoketolase)

    ax[1].plot(x2, y2, '-', label= y2.columns) #
    ax[1].plot(x2, y3, '--', label= y3.columns) # ACL and phosphoketolase

    ax[1].legend(fontsize=10, loc='upper left')
    ax[1].set_title("Intracellular fluxes")
    ax[1].set_xlabel('Biomass growth rate $(1/h)$')
    ax[1].set_ylabel('Flux $(mmol/gDW/h)$')
    
    
    
def cofactor_balances_biomass_max(model_path, cofactor_list, glucose_uptakes, i, biomass_rxn_ID, glc_ID):
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



def cofactor_balances_NGAM_min(model_path, cofactor_list, glucose_uptakes, growth_rates, NGAM_rxn_ID, glc_ID, biomass_rxn_ID, i):
    producing_fluxes = pd.DataFrame() 
    consuming_fluxes = pd.DataFrame()
    
    model=cobra.io.read_sbml_model(model_path)
    
    model.objective = NGAM_rxn_ID 
    model.reactions.get_by_id(glc_ID).bounds = -(glucose_uptakes[i]), -(glucose_uptakes[i])  
    model.reactions.get_by_id(biomass_rxn_ID).bounds = growth_rates[i], growth_rates[i]

        
    model.optimize('minimize')
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




def cofactor_fluxes_pie_chart(model_path, cofactor_fluxes, **fig_kw):
    threshold = 0.05# threshold shows the percent of the flux for including in others sector on pie chart 
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
    plt.tight_layout()

    # plt.legend(producing_and_consuming_fluxes, reaction_names_w_flux, title = 'Reaction names', loc="center left",  bbox_to_anchor=(1, 0, 0.5, 1))
    
    # plt.title(title)

    return pie_chart, fig


# Get all fluxes to excel
def all_fluxes_to_excel(path, all_fluxes_df):
    with pd.ExcelWriter(path) as excel_writer:
        all_fluxes_df.to_excel(excel_writer, sheet_name='Sheet 1', index=True)
        
# Get csv file of fluxes (from 0-4)

def fluxes_to_csv(path, all_fluxes_df, i):
    all_fluxes_df.loc[i].to_csv(path, index=True)