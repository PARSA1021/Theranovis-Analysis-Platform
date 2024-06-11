import anndata
import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import os
from typing import Tuple

# 가장 가까운 이웃을 계산하는 함수 정의
def calculate_NearestNeighbor(subAD: anndata.AnnData, phenotype_1: str, phenotype_2: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    static_df = pd.DataFrame(index=["Min", "Mean", "Median", "Max", "Std"])

    centroidsP1 = subAD.obs[subAD.obs['phenotype'] == phenotype_1][['spatial_X', 'spatial_Y']].values
    centroidsP2 = subAD.obs[subAD.obs['phenotype'] == phenotype_2][['spatial_X', 'spatial_Y']].values

    if centroidsP1.size == 0 or centroidsP2.size == 0:
        print(f"No data found for phenotype_1: {phenotype_1} or phenotype_2: {phenotype_2}")
        return static_df, pd.DataFrame(columns=['Start_index', "Start_point", "End_index", 'End_point', 'distance', "Start -> End"])

    p1p2DF = pd.DataFrame(columns=['Start_index', "Start_point", "End_index", 'End_point', 'distance', "Start -> End"])

    print(f"Calculating distance {phenotype_1} -> {phenotype_2}")
    for ix, row in tqdm.tqdm(subAD.obs.loc[subAD.obs['phenotype'] == phenotype_1].iterrows(), total=subAD.obs.loc[subAD.obs['phenotype'] == phenotype_1].shape[0]):
        c = row[['spatial_X', 'spatial_Y']].values.astype(float)
        distance_array = sp.spatial.distance.cdist([c], centroidsP2)
        distance = distance_array.min()
        nearest_end = np.where(distance_array == distance.min())[1]
        end_point = centroidsP2[nearest_end][0]
        end_index = subAD.obs[(subAD.obs.spatial_X == end_point[0]) & (subAD.obs.spatial_Y == end_point[1])].index[0]
        p1p2DF = p1p2DF.append({"Start_index": ix, "Start_point": c, "End_index": end_index, "End_point": end_point, 'distance': distance, "Start -> End": f"{phenotype_1} -> {phenotype_2}"}, ignore_index=True)

    # 여기에서 phenotype_2에서 phenotype_1로의 거리도 계산합니다.
    print(f"Calculating distance {phenotype_2} -> {phenotype_1}")
    for ix, row in tqdm.tqdm(subAD.obs.loc[subAD.obs['phenotype'] == phenotype_2].iterrows(), total=subAD.obs.loc[subAD.obs['phenotype'] == phenotype_2].shape[0]):
        c = row[['spatial_X', 'spatial_Y']].values.astype(float)
        distance_array = sp.spatial.distance.cdist([c], centroidsP1)
        distance = distance_array.min()
        nearest_end = np.where(distance_array == distance.min())[1]
        end_point = centroidsP1[nearest_end][0]
        end_index = subAD.obs[(subAD.obs.spatial_X == end_point[0]) & (subAD.obs.spatial_Y == end_point[1])].index[0]
        p1p2DF = p1p2DF.append({"Start_index": ix, "Start_point": c, "End_index": end_index, "End_point": end_point, 'distance': distance, "Start -> End": f"{phenotype_2} -> {phenotype_1}"}, ignore_index=True)

    static_df[f"[ {phenotype_1} -> {phenotype_2} ]"] = [p1p2DF.distance.mean(),
                                                        p1p2DF.distance.min(),
                                                        p1p2DF.distance.median(),
                                                        p1p2DF.distance.max(),
                                                        p1p2DF.distance.std()]

    return static_df, p1p2DF

# 화살표를 그리는 함수 정의
def draw_arrow(start, end):
    plt.arrow(start[0], start[1], end[0] - start[0], end[1] - start[1],
              shape='full', lw=0.3, length_includes_head=True, head_width=10, color='black')

# 표적 가장 가까운 이웃을 그리는 함수 정의
def draw_phenotype_nearest_neighbor(newAD: anndata.AnnData, pair_df: pd.DataFrame, pair_start_end: str, cMap: dict, output_dir_path: str, type_string: str, dpi=600, transparent=True):
    plt.figure(figsize=(16, 12))
    plt.axis('off')
    data = newAD.obs[newAD.obs.phenotype != 'Unknown']
    sns.scatterplot(data=data, x='spatial_X', y='spatial_Y', hue='phenotype', palette=cMap, edgecolor="none", alpha=0.9)
    vis = pair_df[pair_df["Start -> End"] == pair_start_end]
    for i in tqdm.tqdm(vis.index):
        start = vis.loc[i, "Start_point"]
        end = vis.loc[i, "End_point"]
        draw_arrow(start, end)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, frameon=False)
    plt.axis('scaled')
    plt.savefig(os.path.join(output_dir_path, f'{type_string}.png'), dpi=dpi, transparent=transparent)
    plt.close()

# 주석을 찾는 함수 정의
def find_leaf(df, sections):
    annotations = []
    for sec in sections:
        for i, row in df.iterrows():
            if row["Parent"] == sec:
                annotations.append(row["Name"])
    return annotations

# QuPath 데이터 읽기
def read_data(qupath_annotation_path: str, qupath_detection_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, list]:
    # 주석 측정 읽기
    annot_df = pd.read_csv(qupath_annotation_path, sep='\t')
    annot_df.reset_index(drop=True, inplace=True)

    # 감지 측정 읽기
    detect_df = pd.read_csv(qupath_detection_path, sep='\t')
    detect_df.reset_index(drop=True, inplace=True)

    # 계층 구조 확인
    sections = annot_df.loc[annot_df["Parent"] == "Root object (Image)"]["Name"]
    annotations = find_leaf(annot_df, sections)
    try:
        if set(detect_df.Parent.unique()) != set(annotations):
            raise ValueError("Annotation file and Detection file did not match. Please check your files.")
    except ValueError as ve:
        print(ve)

    return annot_df, detect_df, annotations
