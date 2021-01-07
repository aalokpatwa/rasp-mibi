#Imports
from pysurvival.models import survival_forest
import numpy as np
import pandas as pd
import seaborn as sns
import os
from pysurvival.models import survival_forest
from pysurvival.utils.metrics import concordance_index
from sklearn.model_selection import train_test_split
import shap
import matplotlib.pyplot as plt

#--------------------------------------------------------------------
#Calculating Shapley values

#Provide the path to the csv with the multivariate data
CSV_PATH = "intermediate_data/covariate_rsf_data.csv"

#Import CSV into a dataframe
data = pd.read_csv(CSV_PATH, index_col="ID")

#Manipulate the architecture feature from a string into a quantitative categorical variable
data["Architecture"] = data["Architecture"].astype("category").cat.codes

#Define the random survival forest model
rf = survival_forest.RandomSurvivalForestModel(num_trees = 200)

#Define predictors and response for inputting into the model
X = data[["functional_proteins_cluster", "immunoregulatory_protein_cluster", "coexpression_cluster", "age", "grade", "Architecture"]]
T = data["Survival_time"]
E = data["Survival"]

#Fit the RSF
fitted = rf.fit(X=X, T=T, E=E, max_depth=5)

#Train the Shap KernelExplainer and calculate Shapley values
explainer = shap.KernelExplainer(fitted.predict_risk, data=X)
shap_values = explainer.shap_values(X)

#Plot the Shapley values as a bar plot
fig = shap.summary_plot(shap_values, features=X, plot_type='bar')

#--------------------------------------------------------------------
#Evaluating model performance

#Define a number of iterations to fit the model on different seeds
iterations = 1000

importance_array = np.zeros((6))
concordances = []


for iter in range(iterations):
    data_train, data_test = train_test_split(data, test_size=0.2, random_state=iter)
    X_train = data_train[["functional_proteins_cluster", "immunoregulatory_protein_cluster", "coexpression_cluster", 
              "age", "grade", "Architecture"]]
    T_train = data_train["Recurrence_time"]
    E_train = data_train["Recurrence"]

    X_test = data_test[["functional_proteins_cluster", "immunoregulatory_protein_cluster", "coexpression_cluster", 
              "age", "grade", "Architecture"]]
    T_test = data_test["Recurrence_time"]
    E_test = data_test["Recurrence"]

    #Fit the RSF according to the training set
    fitted = rf.fit(X=X_train, T=T_train, E=E_train, max_depth=5, seed=iter)
    concordance = concordance_index(rf, X, T, E, include_ties = True, additional_results=False)
    concordances.append(concordance)

mean_concordance = np.mean(np.array(concordances))

print ("Mean concordance index: " + str(round(mean_concordance, 4)))
