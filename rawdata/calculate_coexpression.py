import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import seaborn
from itertools import combinations
from itertools import product

binary_infopath = "intermediate_data/protein_positivity/"
biomarker_frames = pd.read_csv("rawdata/proteins_by_frame.csv")

biom_columns = biomarker_frames["Biomarker"].values

for patient in os.listdir(binary_infopath):
    print (patient)
    com = np.zeros((44,44))
    
    infodf = pd.read_csv(binary_infopath + patient)
        
    n_cells = len(infodf.index)
    
    for cell in range(n_cells):
        this_cell = infodf.iloc[[cell]]
        pos_columns = []
        for column in this_cell.columns[3:]:
            if this_cell[column].values[0] == 1:
                pos_columns.append(int(column))
        combs = combinations(pos_columns, 2)
        for combination in combs:
            com[combination[0], combination[1]] += 1
            com[combination[1], combination[0]] += 1
    com_df = pd.DataFrame(com, columns=biom_columns)
    com_df.set_index(pd.Index(biom_columns), inplace=True)
    com_df.to_csv("intermediate_data/created_coexpression_matrices/" + patient)