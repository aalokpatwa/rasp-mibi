import os
import numpy as np 
import pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

#-----------------------------------------
#Recurrence

expression_path = "intermediate_data/proteinexpression_df.csv"

proteinexpression_df = pd.read_csv(expression_path, index_col="ID")

#indices of the proteins considered to be "functional" proteins, and the clinical information
functional_protein_indices = [2,6,14,16,22,23,24,25,26,27,28,29,30,31,35,37,38,39, 44, 45, 46, 47]

#slice the functional proteins only
functional_df = proteinexpression_df.iloc[:, functional_protein_indices]

#Leave out the clinical data to iterate over the proteins
proteins = functional_df.columns[:-4]

univariate_results = []

for protein in proteins:
    this_protein = functional_df[[protein, "Recurrence_time", "Recurrence"]]
    #Fit the Cox PH model and extract information
    cox = CoxPHFitter()
    cox.fit(this_protein, duration_col="Recurrence_time", event_col="Recurrence")
    summary_df = cox.summary
    coefficients = summary_df["coef"].values[0]
    hazards = cox.hazard_ratios_[0]
    p_values = summary_df["p"].values[0]
    
    univariate_results.append([protein, coefficients, hazards, p_values])

univariate_results_df = pd.DataFrame(univariate_results, columns=["Protein", "Coef", "HR", "P"])
univariate_results_df.sort_values(by="P", ascending=True, inplace=True)

p_values = univariate_results_df["P"]

corrected = multipletests(p_values, method="fdr_bh")[1]

univariate_results_df["BH-Corrected FDR"] = corrected

univariate_results_df.to_csv("results/protein_expression_recurrence.csv", index=False)

#-----------------------------------------
#Survival


univariate_results = []

for protein in proteins:
    this_protein = functional_df[[protein, "Survival_time", "Survival"]]
    #Fit the Cox PH model and extract information
    cox = CoxPHFitter()
    cox.fit(this_protein, duration_col="Survival_time", event_col="Survival")
    summary_df = cox.summary
    coefficients = summary_df["coef"].values[0]
    hazards = cox.hazard_ratios_[0]
    p_values = summary_df["p"].values[0]
    
    univariate_results.append([protein, coefficients, hazards, p_values])

univariate_results_df = pd.DataFrame(univariate_results, columns=["Protein", "Coef", "HR", "P"])
univariate_results_df.sort_values(by="P", ascending=True, inplace=True)

p_values = univariate_results_df["P"]

corrected = multipletests(p_values, method="fdr_bh")[1]

univariate_results_df["BH-Corrected FDR"] = corrected

univariate_results_df.to_csv("results/protein_expression_survival.csv", index=False)