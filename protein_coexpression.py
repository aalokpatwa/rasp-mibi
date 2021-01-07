import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import seaborn


#Path to the clinical data about recurrence. 
clinical_path = "rawdata/clinical_data.csv"
clinical_df = pd.read_csv(clinical_path, index_col=["ID"])

#Path to interaction matrices
matrices_path = "/Users/aalokpatwa/Desktop/AalokMacPro/MIBI/nonvoronoi/within_cell_com/cellsize_norm_background/"

# Out of all proteins, only include the functional proteins
markers_to_include = [2,6,14,16,22,23,24,25,26,27,28,29,30,31,35,37,38,39]

columns = []
feature_list = []

#Iterate over all of the co-occurrence matrices in this certain radius. 
for patient_index in range(len(os.listdir(matrices_path))):
    patient_glcm = os.listdir(matrices_path)[patient_index]

    #Skip over the pesky .DS_Store file that shows up in Mac file systems.
    if (patient_glcm[0] == "."):
        continue

    #Find the internal_ID of the current patient. 
    identifier = patient_glcm.split(".")[0]

    #Read the co-occurrence matrix of the current patient and cast to a numpy array excluding the first col.
    current_patient_glcm = pd.read_csv(os.path.join(matrices_path, patient_glcm), index_col="Unnamed: 0")
    current_patient_glcm.set_index(current_patient_glcm.columns, inplace=True)

    #feature_name = current_patient_glcm.columns[chosen_feature]
    np_glcm = current_patient_glcm.to_numpy()    

    patient_features = []

    #Flatten the co-occurrence matrix into a feature vector.
    #Cannot simply run np.flatten because the matrix is symmetrical - I only want the top-right triangle
    #But, have to include the diagonal as well.
    for row in range(np_glcm.shape[0]):
        for column in range(row+1, np_glcm.shape[1]):
            if row in markers_to_include and column in markers_to_include:
                if patient_index == 0:
                    feature_name = current_patient_glcm.index[row] + " " + current_patient_glcm.columns[column]
                    columns.append(feature_name)
                patient_features.append(np_glcm[row][column])

    #Find the recurrence outcome of this current patient and how long it took for them to recur
    #Add these values to the feature list for eventual use in the Kaplan-Meier plot.
    try:
        patient_features.append(clinical_df.at[int(identifier), "Recurrence"])
    except KeyError:
        continue
    patient_features.append(clinical_df.at[int(identifier), "Recurrence_time"])
    patient_features.append(int(identifier))
    feature_list.append(patient_features)

#Determine the names of the columns in the DataFrame for easier future access.
columns.append("Recurrence")
columns.append("Recurrence_time")
columns.append("ID")

#Create a dataframe using the features, the recurrence events, and the time taken to recur.

features_df = pd.DataFrame(feature_list, columns=columns)

#Obtain a versino of this dataframe with only the features.
data_only = features_df.drop(columns=["Recurrence_time", "Recurrence", "ID"])
data_only["Duplicate"] = features_df[feature_name]
data_only.set_index(features_df["ID"], inplace=True)

#Create the dendrogram.
clustergram = seaborn.clustermap(data_only, method="weighted",
                            metric="canberra", standard_scale = 1, cmap="viridis", figsize=(10,8), cbar_pos=None)

#Number of clusters to take from the dendrogram
k = 2

#Use scipy fcluster to find clusters from the clustered dendrograms.
clusters = list(fcluster(clustergram.dendrogram_row.linkage, k, criterion='maxclust'))

unique_clusters = len(np.unique(np.array(clusters)))

if (unique_clusters < 2):
    print ("Not able to find more than one cluster. ") #If there was no clustering, this method failed.

#Create a new column in the DataFrame that includes what cluster each patient falls into
features_df["clust"] = clusters

first_cluster_count = clusters.count(1)
second_cluster_count = clusters.count(2)

#Define the KaplanMeierFitter
kmf = KaplanMeierFitter()

#T = Time, E=Event. These are the two parameters that go into Kaplan-Meier curves. 
T = features_df["Recurrence_time"]
E = features_df["Recurrence"]

group1 = (features_df["clust"] == 1)
group2 = (features_df["clust"] == 2)

T1 = T[group1]
E1 = E[group1]
T2 = T[group2]
E2 = E[group2]

color_clust1 = "#F39B7FFF"
color_clust2 = "#4DBBD5FF"

# Just for visualization purposes, make the worse-outcome cluster orange and the other blue
if E1.mean() < E2.mean():
    T1 = T[group2]
    E1 = E[group2]
    T2 = T[group1]
    E2 = E[group1]

features_df["clust"][group1] = 2
features_df["clust"][group2] = 1

first_cluster_count = len(T1.index)
second_cluster_count = len(T2.index)

results_first = logrank_test(T1,T2, event_observed_A=E1, event_observed_B=E2)

p1 = round(results_first.p_value, 4)

print ("Gave two clusters with counts (" + 
    str(first_cluster_count) + ", " + 
    str(second_cluster_count) + 
    ") and yielded a log-rank of " + 
    str(p1))

plt.figure(figsize=(10,8))
kmf.fit(T1, E1, label='Cluster 1: n=' + str(first_cluster_count))
ax = kmf.plot(ci_show=False, show_censors=True, color=color_clust1, lw=5)

kmf.fit(T2, E2, label='Cluster 2: n=' + str(second_cluster_count))
ax = kmf.plot(ax=ax, ci_show=False, show_censors=True, color=color_clust2, lw=5)
plt.annotate("Log-rank p: " + str(p1), xy=(0.6, 0.15), xycoords="figure fraction", fontsize=25,
        bbox=dict(facecolor='none', edgecolor='black', alpha=0.3, boxstyle="Round, pad=0.5, rounding_size=0.2"))
plt.legend(fontsize=25)
plt.xlabel("Time (days)", fontsize=25)
plt.xticks(fontsize=15)
plt.ylabel("Proportion Alive", fontsize=25)
plt.yticks(fontsize=15)
plt.ylim(0,1)
plt.tight_layout()

plt.savefig("results/coexpression_km_recurrence.png", dpi=300)
plt.close()