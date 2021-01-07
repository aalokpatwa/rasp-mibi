from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
import lifelines
import matplotlib.pyplot as plt
import numpy as np
import sklearn
import glob, os
import pandas as pd
from PIL import Image, ImageDraw
from scipy.stats import ttest_ind
import scipy
from scipy.cluster.hierarchy import fcluster
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import cv2
import skimage
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import time
import seaborn
from itertools import combinations
from itertools import product
import networkx as nx
import json
import karateclub
import pickle
import keras
import random

def create_voronoi(centdf): #returns the Voronoi diagram object
    vor = Voronoi(np.c_[centdf.column.values, centdf.row.values])
    return vor
def plot_voronoi(vor, identifier): #Void - just plots the graph
    fig, ax = plt.subplots(figsize=(16, 16))

    fig = voronoi_plot_2d(vor, ax, show_vertices=False, line_colors='blue', 
                      line_width=2, point_size=2)
    ax.set_xlim([0, 2048])
    ax.set_ylim([2048, 0])
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("/Users/aalokpatwa/Desktop/MIBI/voronoi/voronoi_plotsv2/" + identifier + ".png", dpi=300)
    plt.close()
    #plt.show()

def make_voronoi_dataframe(voronoi, centroiddf):
    count_skipped = 0
    number_centroids = voronoi.points.shape[0]
    information_list = []
    for centroid in range(number_centroids):
        centroid_list = []
        centroid_y = voronoi.points[centroid][0]
        centroid_x = voronoi.points[centroid][1]
        centroid_list.append(centroid)
        centroid_list.append((centroid_x, centroid_y))
        region_index = voronoi.point_region[centroid]
        centroid_list.append(region_index)
        vertices = voronoi.regions[region_index]
        if -1 in vertices:
            count_skipped += 1
            continue
        else:
            vertex_list = []
            for vertex in vertices:
                y_coord = int(round(voronoi.vertices[vertex][0]))
                x_coord = int(round(voronoi.vertices[vertex][1]))
                vertex_list.append((y_coord, x_coord))
            centroid_list.append(vertex_list)
        centroid_list.append(centroiddf.at[centroid, "celltype"])
        information_list.append(centroid_list)       
    infodf = pd.DataFrame(information_list, columns=["CentroidIndex","CentroidCoord", "RegionIndex",
                                                     "Vertices", "Celltype"])
    return infodf, count_skipped

def extract_current_mask(vertexlist):
    img = Image.new("L", (2048,2048), 0)
    ImageDraw.Draw(img).polygon(vertexlist, fill=1)
    mask = np.array(img)
    binary_mask = np.where(mask==1)
    coordinates = list(zip(binary_mask[0], binary_mask[1]))
    return coordinates, len(coordinates)

def expression_within_cell(biomarker_image, pixellist, size):
    total_vector = np.zeros(44)
    for pixel in pixellist:
        expression_vector = biomarker_image[pixel[0], pixel[1]]
        total_vector += expression_vector
    total_vector = total_vector / size
    return total_vector

def create_multimarker_image(imagepath):
    biomarker_pil = Image.open(imagepath)
    n_frames = biomarker_pil.n_frames
    
    first_level = np.array(biomarker_pil, dtype="uint8")
    first_level = np.reshape(first_level, (2048,2048,1))
        
    combined_image = first_level
    
    for frame in range(1, n_frames):
        biomarker_pil.seek(frame)
        current_level = np.array(biomarker_pil, dtype="uint8").reshape((2048,2048,1))
        combined_image = np.concatenate((combined_image, current_level), axis=2)
    return combined_image

def create_neighbor_matrix(voronoi, voronoi_df):
    number_regions = len(voronoi_df.index)
    adjacency_list = []
    for region in range(number_regions):
        adjacency_list.append([])
    ridge_points = voronoi.ridge_points
    for edge in ridge_points:
        first_centroid = edge[0]
        second_centroid = edge[1]
        first_cell = voronoi_df[voronoi_df["CentroidIndex"] == first_centroid]
        second_cell = voronoi_df[voronoi_df["CentroidIndex"] == second_centroid]
        if (first_cell.empty or second_cell.empty):
            continue
        first_index = int(first_cell.index[0])
        second_index = int(second_cell.index[0])
        adjacency_list[first_index].append(second_centroid)
        adjacency_list[second_index].append(first_centroid)
    new_df = voronoi_df.copy()
    new_df["Adjacency"] = adjacency_list
    return new_df


# Iterate through the ridge points because they define a certain combination
# Find the indices of the cells that the border indicates
# Find their Binary Expression matrix using the within_cell_binary_median
# Find the possible combinations from two distinct sets, and then add one to each of the corresponding selections.

#IMPORTANT - THE PURPOSE OF THIS CELL IS TO USE NON VORONOI MARKER EXPRESSION TO RUN THIS

binary_infopath = "/Users/aalokpatwa/Desktop/MIBI/nonvoronoi/marker_expression/cellsize_binary_background/"
centroid_path = "/Users/aalokpatwa/Desktop/MIBI/centroids/"
biomarker_frames = pd.read_csv("/Users/aalokpatwa/Desktop/MIBI/biomarker_frames.csv")

biom_columns = biomarker_frames["Biomarker"].values

for patient in os.listdir(binary_infopath):
    print (patient)
    if patient[0] == ".":
        continue
    com = np.zeros((44,44))
    centroid_df = pd.read_csv(centroid_path + patient)
    biom_df = pd.read_csv(binary_infopath + patient)
    vor = create_voronoi(centroid_df)
    edges = vor.ridge_points
    for edge in edges:
        first_centroid = int(edge[0])
        second_centroid = int(edge[1])
        
        first_cell = biom_df.loc[[first_centroid]]
        second_cell = biom_df.loc[[second_centroid]]
        
            
        
        first_pos = []
        second_pos = []
        for column in first_cell.columns[3:]:
            if first_cell[column].values[0] == 1:
                first_pos.append(int(column))
            if second_cell[column].values[0] == 1:
                second_pos.append(int(column))
        combinations = product(first_pos, second_pos)
        for comb in combinations:
            first_marker = comb[0]
            second_marker = comb[1]
            com[first_marker, second_marker] += 1
            if (first_marker != second_marker):
                com[second_marker, first_marker] += 1
    com_df = pd.DataFrame(com, columns=biom_columns)
    com_df.set_index(pd.Index(biom_columns), inplace=True)
    com_df.to_csv("/Users/aalokpatwa/Desktop/MIBI/voronoi/marker_cooccurrence/biom_com/cellsize_nonvoronoi_background/" + patient)