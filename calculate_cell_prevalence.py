# NOTE: This code cannot be run without the raw MIBI data. The data is not included with this repo.
# As such, this code is simply to understand the flow of the analysis, not to be run.

import os
import cv2
import numpy as np 
import PIL 
import pandas as pd
from skimage import measure
import scipy

imagepath = "rawdata/TNBCcelltypes"

#Keep track of the counts using a matrix of size n_patients X n_cells
all_info = []

#Iterate through all patients
for patient in os.listdir(imagepath):

    #Find out the recurrence label of this patient
    identifier = int(patient.split("_")[0][1:])
    
    if identifier == 30 or identifier == 22 or identifier == 38:
        continue
    
    #Read image
    image = cv2.imread(os.path.join(imagepath,patient), 0)
    
    #Connect components
    labeled, total_cells = measure.label(image, return_num=True)
    
    #Find the number of tumor cells in the image (for future reference)
    nr_tumor = measure.label((image == 4), return_num=True)[1]
    
    information = []
    information.append(identifier)
    
    information.append(total_cells)
    
    #Find the counts of all celltypes
    for celltype in range(2,17):
        selected_cell = (image == celltype)
        labeled, nr = measure.label(selected_cell, return_num=True)
        nr /= total_cells
        information.append(nr)
    
    
    all_info.append(information)
   
#Create a DataFrame of the counts with the names of the cells rather than numbers 
countdf = pd.DataFrame(all_info, columns=["ID", "Total Cells", "Endothelial", "Mesenchyme",
                                         "Tumor", "Treg", "CD4_T", "CD8_T", "CD3_T", "NK",
                                         "B", "Neutrophil", "Macrophage", "DC", "DC_Mono",
                                         "Mono_neutrophil", "Other", "Label", "Days"])

countdf.set_index("ID", inplace=True)
countdf.sort_index(ascending=True, inplace=True)
countdf.to_csv("intermediate_data/created_cellprevalence_df.csv")