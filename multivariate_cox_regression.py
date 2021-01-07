#Imports
import numpy as np
import pandas as pd
import lifelines
import seaborn as sns

CSV_PATH = "intermediate_data/covariate_rsf_data.csv" #put path to csv summarizing cluster features, clinical variables, architecutre distinction, and clinical outcome

#Read csv summarizing cluster features, clinical variables, morphology distinction, and clinical outcome
df = pd.read_csv(CSV_PATH, index_col="ID")

#define which cluster feature to examine
#options: [coexpression_cluster, functional_proteins_cluster, immunoregulatory_protein_cluster]
#can either examine each one-at-a-time or iterate through them
cluster_choice = "coexpression_cluster" 

#manipulate the architecture distinction the variable which is originally dtype: str into a quantitative categorical variable
df["Architecture"] = df["Architecture"].astype("category").cat.codes

#Create two separate feature matrices for recurrence and survival
recurrence_df = df.drop(columns=["Survival", "Survival_time"])[[cluster_choice, "grade", "age", "Architecture", "Recurrence", "Recurrence_time"]]
survival_df = df.drop(columns=["Recurrence", "Recurrence_time"])[[cluster_choice, "grade", "age", "Architecture", "Survival", "Survival_time"]]

#Define and fit Cox PH Fitter for recurrence
recurrence_cph = lifelines.CoxPHFitter()
recurrence_cph.fit(recurrence_df, duration_col='Recurrence_time', event_col='Recurrence')
recurrence_summary = recurrence_cph.print_summary()

#Define and fit Cox PH Fitter for survival
survival_cph = lifelines.CoxPHFitter()
survival_cph.fit(survival_df, duration_col='Survival_time', event_col='Survival')
survival_summary = survival_cph.print_summary()