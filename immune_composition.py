import os
import numpy as np 
import pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests


#-----------------------------------------
#Recurrence

cellprevalence_path = "intermediate_data/cellprevalence_df.csv"

cellprevalence_df = pd.read_csv(cellprevalence_path, index_col="ID")

#Get a list of all the cell types to test
each_celltype = cellprevalence_df.columns[:15]

univariate_results = []

for celltype in each_celltype:
    
    this_cell = cellprevalence_df[[celltype, "Recurrence_time", "Recurrence"]]
    #Perform univariate Cox regression for this cell type
    cox = CoxPHFitter()
    cox.fit(this_cell, duration_col="Recurrence_time", event_col="Recurrence")
    summary_df = cox.summary
    coefficients = summary_df["coef"].values[0]
    hazards = cox.hazard_ratios_[0]
    p_values = summary_df["p"].values[0]
    
    univariate_results.append([celltype, coefficients, hazards, p_values])

univariate_results_df = pd.DataFrame(univariate_results, columns=["Celltype", "Coef", "HR", "P"])
univariate_results_df.sort_values(by="P", ascending=True, inplace=True)

p_values = univariate_results_df["P"]

corrected = multipletests(p_values, method="fdr_bh")[1]

univariate_results_df["BH-Corrected FDR"] = corrected

univariate_results_df.to_csv("results/immune_composition_recurrence.csv", index=False)

#-----------------------------------------
#Recurrence

cellprevalence_path = "intermediate_data/cellprevalence_df.csv"

cellprevalence_df = pd.read_csv(cellprevalence_path, index_col="ID")

#Get a list of all the cell types to test
each_celltype = cellprevalence_df.columns[:15]

univariate_results = []

for celltype in each_celltype:
    
    this_cell = cellprevalence_df[[celltype, "Survival_time", "Survival"]]
    #Perform univariate Cox regression for this cell type
    cox = CoxPHFitter()
    cox.fit(this_cell, duration_col="Survival_time", event_col="Survival")
    summary_df = cox.summary
    coefficients = summary_df["coef"].values[0]
    hazards = cox.hazard_ratios_[0]
    p_values = summary_df["p"].values[0]
    
    univariate_results.append([celltype, coefficients, hazards, p_values])

univariate_results_df = pd.DataFrame(univariate_results, columns=["Celltype", "Coef", "HR", "P"])
univariate_results_df.sort_values(by="P", ascending=True, inplace=True)

p_values = univariate_results_df["P"]

corrected = multipletests(p_values, method="fdr_bh")[1]

univariate_results_df["BH-Corrected FDR"] = corrected

univariate_results_df.to_csv("results/immune_composition_survival.csv", index=False)