import os
import pandas as pd
import numpy as np
import scipy
import math
import plotnine
import matplotlib.pyplot as plt
import warnings
from scipy.integrate import solve_ivp
from datetime import datetime
from datetime import timedelta

## file paths
p_half_sats = "~/data/arianna/2026-krinos-diatom-model/2026_08_29_complete_workflow/output_estimation_logN_linP_20261009/half_sats/phosphorus_half_sats.csv"
n_half_sats = "~/data/arianna/2026-krinos-diatom-model/2026_08_29_complete_workflow/output_estimation_logN_linP_20261009/half_sats/nitrogen_half_sats.csv"
temp_relation_file = "~/data/arianna/2026-krinos-diatom-model/2026_08_29_complete_workflow/output_estimation_logN_linP_20261009/temp_relation/temperature_mean_std.csv"

## read in allometric scaling relationships
scaling_relationships = pd.read_csv("../../data/scaling_relationships_calc.csv")

## read in estimated half-saturation constants from prior workflow
estimated_half_sats_p = pd.read_csv(p_half_sats)
estimated_half_sats_p = estimated_half_sats_p.loc[estimated_half_sats_p.exponential_p_val < 0.11,:]
estimated_half_sats_p=estimated_half_sats_p.loc[estimated_half_sats_p.exponential_slope<=0,:]
estimated_half_sats_n = pd.read_csv(n_half_sats)
estimated_half_sats_n = estimated_half_sats_n.loc[estimated_half_sats_n.exponential_p_val < 0.11,:]
estimated_half_sats_n=estimated_half_sats_n.loc[estimated_half_sats_n.exponential_slope<=0,:]

def get_half_sat(Taxon,nutrient="P",metaT_or_allom="metaT"):
    if metaT_or_allom=="allom":
        if nutrient=="P":
            return scaling_relationships.loc[scaling_relationships.taxon==Taxon,
                                   "half_sat_p"].iloc[0]
        return scaling_relationships.loc[scaling_relationships.taxon==Taxon,
                               "half_sat_n"].iloc[0]
    
    if nutrient=="P":
        if Taxon in list(estimated_half_sats_p.Taxon):
            return float(estimated_half_sats_p.loc[estimated_half_sats_p.Taxon==Taxon,
                                   "estimated_parameter"].iloc[0])
        return scaling_relationships.loc[scaling_relationships.taxon==Taxon,
                                   "half_sat_p"].iloc[0]
    
    if Taxon in list(estimated_half_sats_n.Taxon):
        return float(estimated_half_sats_n.loc[estimated_half_sats_n.Taxon==Taxon,
                                   "estimated_parameter"].iloc[0])
    return scaling_relationships.loc[scaling_relationships.taxon==Taxon,
                               "half_sat_n"].iloc[0]


temp_relation=pd.read_csv(temp_relation_file)
def get_temp_params(Taxon,param="mean"):
    if len(temp_relation.loc[temp_relation.Order==Taxon,"mean_temp"])==0 & (param=="mean"):
        return 8
    if len(temp_relation.loc[temp_relation.Order==Taxon,"mean_temp"])==0 & (param=="std"):
        return 5
    if param=="mean":
        if len(temp_relation.loc[temp_relation.Order==Taxon,"mean_temp"])==0 | \
            (temp_relation.loc[temp_relation.Order==Taxon,"mean_temp"].iloc[0]\
           !=temp_relation.loc[temp_relation.Order==Taxon,"mean_temp"].iloc[0]):
            return 8
        return float(temp_relation.loc[temp_relation.Order==Taxon,"mean_temp"])
    if param=="std":
        if len(temp_relation.loc[temp_relation.Order==Taxon,"std_temp"])==0 | \
           (temp_relation.loc[temp_relation.Order==Taxon,"std_temp"].iloc[0]\
           !=temp_relation.loc[temp_relation.Order==Taxon,"std_temp"].iloc[0]):
            return 5
        return float(temp_relation.loc[temp_relation.Order==Taxon,"std_temp"])
    
def get_half_sats_model(estimated_half_sats_n,estimated_half_sats_p):
    return list(set(estimated_half_sats_n.Taxon).union(set(estimated_half_sats_p.Taxon))) + ["Coscinodiscales"]

## use the metaT-informed parameters to decide the list of organisms in the model
organisms_for_model = get_half_sats_model(estimated_half_sats_n,estimated_half_sats_p) 


### standard parameter values
consistent_k_light = 500
consistent_T_opt = 8
consistent_std_T = 5
consistent_max_gr = 1
consistent_half_sat_silicon = 1
consistent_half_sat_metab = 0.5
consistent_half_sat_nitr = 1
consistent_half_sat_phos = 0.1

## MAIN MODEL FUNCTIONS

def identify_limiting_nutrient(mm_n,mm_p,mm_s,mm_m):
    return min(mm_n,mm_p,mm_s,mm_m)

## notes - now doing this per day - but continue to include light variability over the year
## by averaging the total light over the course of the day (lower average during darker months)
def calculate_growth_rate(k_light, mu_max, T_opt, std_T, 
                          limiting_resource, Temp, Light, mortality_rate):
    light_sensitivity = 1 #Light / (k_light + Light)
    temp_sensitivity = np.exp(-((Temp - T_opt)/(std_T*2))**2)
    GR = mu_max * limiting_resource * light_sensitivity * temp_sensitivity
    return GR

def calculate_limiting_original(half_sat, Resource, maxResource):
    if Resource < 0:
        Resource = 0
    U = Resource / (half_sat + Resource)
    return U

def calculate_limiting(half_sat, Resource):
    
    if Resource < 0:
        Resource = 0
    beta = 0.5
    inhibition_constant = half_sat * 4 # inhibition starts at 2x the saturation val.
    U = Resource / (Resource + half_sat)
    #U = (Resource)/ (half_sat + Resource + (Resource**2) / inhibition_constant)

    if (U < 0)|(Resource<0.001):
        U = 0
    return U

def calculate_all_limiting(organism_list, parameters, resource_code,
                         Resource, temp_curr, light_curr):
    limiting_list=[]
    for organism in organism_list:
        if resource_code == "Nitr":
            half_sat_nitr=parameters[organism]["half_sat_nitr"]
            limiting_list.append(calculate_limiting(half_sat_nitr, 
                                                Resource))
        elif resource_code == "Phos":
            half_sat_phos=parameters[organism]["half_sat_phos"]
            limiting_list.append(calculate_limiting(half_sat_phos, 
                                                Resource))
   
    return limiting_list

def holling_type_three(half_sat_graz, phyto_density, refuge_phyto, graz_max):
    return graz_max * (phyto_density**2) / \
        (half_sat_graz**2+(phyto_density**2))
  
def ivlev_grazing(half_sat_graz, saturation_rate, 
                  grazing_preference,
                  phyto_density, graz_max):
    return graz_max  * grazing_preference * (1-np.exp(-saturation_rate*phyto_density))

# P = current phyto conc
# graz_max = max grazing 
# saturation_rate = Ivlev grazing saturation constant
# alpha = switching exponent
def kill_the_winner_ivlev(P, graz_max, saturation_rate, alpha=2.0):
    
    P = np.array(P, dtype=float)
    
    # Avoid division by zero if all phytoplankton are completely depleted
    if np.sum(P) == 0:
        return np.zeros_like(P)
    
    # 1. Calculate density dependence of predation
    numerator = (P ** alpha)
    p_dynamic = numerator / np.sum(numerator)
    
    # 2. Calculate total effective food that zoops have access to
    F = np.sum(p_dynamic * P)
    
    # 3. Calculate multi-prey Ivlev ingestion rate for each species
    g_i = graz_max * p_dynamic * (1.0 - np.exp(-saturation_rate * F))
    
    return g_i

def transporter_model_all_resource(time, state, organism_list, parameters, time_vect, sub_frame_temp):
    if not np.isfinite(time):
        return np.zeros_like(state)  # Return zero derivatives to safely stop or reject the step
    curr_index = int(time)
        
    ## read in environmental variables
    light_curr = 1 # light no longer incorporated
    temp_curr = sub_frame_temp.\
        loc[sub_frame_temp["index"]==curr_index,\
            "temperature_daily"].iloc[0] # temperature
    
    nitr_curr=state[-2]
    phos_curr=state[-1]
    
    chl_curr = 10 # ug/L
    mortality_rate = 0.2
    zoop_mortality = 0.01
    graz_max = 0.5
    half_sat_graz = 0.5
    saturation_rate = 0.693/half_sat_graz
    refuge_phyto = 0.25 # population at which grazing is alleviated (they're hiding)
    
    alpha = 1.25 # scaling to control "kill the winner" dynamics
    
    half_sat_Phy = 0.5
    assimilation = 0.4

    limiting_vect_Nitr = calculate_all_limiting(organism_list,
                                                parameters, resource_code = "Nitr",
                             Resource=nitr_curr,
                             temp_curr=temp_curr, light_curr=light_curr)
    
    limiting_vect_Phos = calculate_all_limiting(organism_list,
                                                parameters, resource_code = "Phos",
                             Resource=phos_curr,
                             temp_curr=temp_curr, light_curr=light_curr)
    
    growth_rates=[]
    grazing_amount=[]
    growth_multipliers=[]
    increase_phytos=[]
    d_phytos=[]
    overall_phyto_density = sum(state[0:-3])
    grazing_rates = kill_the_winner_ivlev(state[0:-3],graz_max,saturation_rate,alpha=alpha)
    for curr_org,ind_org in zip(organism_list,range(len(organism_list))):
        pref_denom = [state_curr**alpha for state_curr in state[0:-3]]
        grazing_preferences = state[ind_org]**alpha / sum(pref_denom)
        mu_max = parameters[curr_org]["max_gr"]
        
        limiting_resource = min(limiting_vect_Nitr[ind_org],
                                limiting_vect_Phos[ind_org])
            
        if limiting_vect_Nitr[ind_org] == limiting_resource:
            growth_multiplier = limiting_vect_Nitr[ind_org] 
        else:
            growth_multiplier = limiting_vect_Phos[ind_org]
        
        growth_multipliers.append(growth_multiplier)
        #limiting_resources.append([time,limiting_resource])
        k_light = consistent_k_light
        T_opt = parameters[curr_org]["T_opt"]
        std_T = parameters[curr_org]["std_T"]
        gr_curr=calculate_growth_rate(k_light, mu_max, T_opt, 
                                      std_T, limiting_resource, temp_curr,
                                      light_curr, mortality_rate)
        
        growth_rates.append(gr_curr)
        increase_phytos.append(gr_curr*(state[ind_org]+0.01))
        d_phytos.append(gr_curr*(state[ind_org]+0.01)-\
                        grazing_rates[ind_org]*state[-3]-\
                        mortality_rate*state[ind_org]**2)
        grazing_amount.append(grazing_rates[ind_org]*state[-3])

        
    
    d_zoop = assimilation * sum(grazing_amount)-\
             zoop_mortality * (state[-3]+0.01)
    d_nitr = ((mortality_rate * sum(state[0:-3]**2) +
              (1-assimilation) * sum(grazing_amount) + \
              zoop_mortality * state[-3] - \
              sum(increase_phytos)))*16/17
    d_phos = ((mortality_rate * sum(state[0:-3]**2) +
              (1-assimilation) * sum(grazing_amount) + \
              zoop_mortality * state[-3] - \
              sum(increase_phytos)))*1/17
    
             #sum(grazing_rates) * (state[-3]+0.01) - \
               #sum(grazing_rates) * state[-3] + 
    return d_phytos + [d_zoop,d_nitr,d_phos]

### FUNCTION FOR INFILLING DATA ON DAYS TO GET NEAREST MONTH ###

def diff_month(d1, d2):
    return (d1.year - d2.year) * 12 + d1.month - d2.month

def infill_data(ptown_data_input,variable_in):
    for curr in range(0,len(ptown_data_input.index)):
        if ptown_data_input[variable_in][curr]!=ptown_data_input[variable_in][curr]: ## if data missing
            variable_differences = pd.DataFrame({"var_diff": ptown_data_input["DateMonthYr"]-ptown_data_input.DateTime[curr],
                                                 "date_month_yr": ptown_data_input.DateMonthYr,
                                                 "month": ptown_data_input.Month,
                                                 "variable": ptown_data_input[variable_in],
                                                 "index": list(range(0,len(ptown_data_input.index)))})
            variable_differences=variable_differences.drop_duplicates(subset=["date_month_yr"])
            variable_differences = variable_differences.loc[variable_differences["variable"]==variable_differences["variable"],:]
            variable_differences_1=variable_differences.copy(deep=True)
            variable_differences_1["abs_diff"] = np.abs(variable_differences_1.var_diff)
            variable_differences_1 = variable_differences_1.loc[[curr>timedelta(0) for curr in variable_differences_1.var_diff],:]
            variable_differences_1 = variable_differences_1.sort_values("abs_diff",ascending=True).head(n=1) ## get the closest values
           
            variable_differences_2=variable_differences.copy(deep=True)
            variable_differences_2["abs_diff"] = np.abs(variable_differences_2.var_diff)
            variable_differences_2 = variable_differences_2.loc[[curr<timedelta(0) for curr in variable_differences_2.var_diff],:]
            variable_differences_2 = variable_differences_2.sort_values("abs_diff",ascending=True).head(n=1)
            variable_differences_2 = variable_differences_2.sort_values("date_month_yr")
            
            variable_differences_1 = pd.concat([variable_differences_1,variable_differences_2]).reset_index()

            value_out = variable_differences_1["variable"][0] +\
                (variable_differences_1["variable"][1] - variable_differences_1["variable"][0]) / \
                diff_month(variable_differences_1.date_month_yr[0],
                           variable_differences_1.date_month_yr[1])
            
            if (variable_in=="temperature")&(ptown_data_input.Month[curr]):
                value_out = variable_differences_1["variable"][0] -\
                (variable_differences_1["variable"][1] - variable_differences_1["variable"][0]) / \
                diff_month(variable_differences_1.date_month_yr[0],
                           variable_differences_1.date_month_yr[1])

            ptown_data_input[variable_in][curr] = value_out
    return ptown_data_input