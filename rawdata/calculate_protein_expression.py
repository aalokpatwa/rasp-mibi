# NOTE: This code cannot be run without the raw MIBI data. The data is not included with this repo.
# As such, this code is simply to understand the flow of the analysis, not to be run.


from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from PIL import Image, ImageDraw
from scipy.cluster.hierarchy import fcluster
import seaborn
from itertools import combinations
from itertools import product

'''
celltype_dir = "rawdata/TNBCcelltypes/"
marker_dir = "rawdata/MIBITIFF/"
marker_frames = "rawdata/proteins_by_frame.csv"
savedir = "intermediate_data/protein_expression/"

frames = pd.read_csv(marker_frames)

columns = frames.index.values.tolist()
columns.insert(0, "Celltype")
columns.insert(0, "Column")
columns.insert(0, "Row")


for patient in os.listdir(celltype_dir):
    identifier = patient.split("_")[0][1:]
    
    marker_path = marker_dir + "p" + identifier + ".tif"
    combined_im = create_multimarker_image(marker_path)
    cv2_im = cv2.imread(os.path.join(celltype_dir, patient), cv2.IMREAD_GRAYSCALE)

    labeled, nr_objects = skimage.measure.label(cv2_im, return_num=True)
    props = skimage.measure.regionprops(labeled)
    
    all_cells = []
    for cell in range(nr_objects):
        centroid_x = int(props[cell].centroid[0])
        centroid_y = int(props[cell].centroid[1])
        
        celltype = cv2_im[centroid_x, centroid_y]

        binary_mask = (labeled == cell)
        biom_expression = combined_im[binary_mask]
        size = biom_expression.shape[0]
        expression_vector = np.sum(biom_expression, axis=0) / size
        expression_vector = expression_vector.tolist()
        expression_vector.insert(0, celltype)
        expression_vector.insert(0, centroid_y)
        expression_vector.insert(0, centroid_x)
        all_cells.append(expression_vector)
    
    dataframe = pd.DataFrame(all_cells, columns=columns)
    dataframe.to_csv(savedir + identifier + ".csv", index=False)



#----------------------------------------
# Assign cells as "positive" or negative for each protein based on their expression level

#Calculate the expression levels in the background (interstitial space)
celltype_dir = "rawdata/TNBCcelltypes/"
marker_dir = "rawdata/MIBITIFF/"

expression_vector = np.zeros(44)
total_pixels = 0

for patient in os.listdir(celltype_dir):
    identifier = patient.split("_")[0][1:]
    if identifier == "30" or identifier == "22" or identifier == "38":
        continue
    
    marker_path = marker_dir + "p" + identifier + ".tif"
    combined_im = create_multimarker_image(marker_path)
    cv2_im = cv2.imread(os.path.join(celltype_dir, patient), cv2.IMREAD_GRAYSCALE)
    dim = cv2_im.shape[0]

    #The TIFFs have a 29-pixel black border surrounding the image which is not a part of the sample. Exclude it.
    cells_only = cv2_im[29:dim-29, 29:dim-29]
    background_mask = (cells_only == 0)
    combined_im = combined_im[29:dim-29, 29:dim-29]
    
    #Extract just one cell
    selection = combined_im[background_mask]
    size = selection.shape[0]
    total_pixels += size
    expression_vector += np.sum(biom_expression, axis=0)

final_vector = expression_vector / total_pixels

'''


#Binarize the matrices according to the threshold calculated across all the images.

expressiondir = "intermediate_data/protein_expression/"
savedir = "intermediate_data/created_protein_positivity/"

threshold = final_vector
for patient in os.listdir(expressiondir):
    patient_glcm = pd.read_csv(expressiondir + patient, index_col="Unnamed: 0")
    for column in range(len(patient_glcm.columns[3:])):
        col_threshold = threshold[column]
        patient_glcm[patient_glcm.columns[3+column]] = (patient_glcm[patient_glcm.columns[3+column]] > col_threshold).astype("int32")
    patient_glcm.to_csv(savedir + patient, index=False)
