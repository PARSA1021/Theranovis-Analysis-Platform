import pandas as pd
import numpy as np
import scanpy as sc
import scipy as sp
import anndata as ad
import tqdm
from typing import Tuple
import streamlit as st
import h5py
from anndata import read_h5ad
import base64

def computeTop20Btm10(ad):
    '''
    Compute the ratio of top 20th percentile to bottom 10th percentile for each protein in each sample
    Input: anndata object
    Output: a dataframe with sampleID, Protein, ratio of top20/btm10
    '''
    top20btm10DF = pd.DataFrame(columns=['ImageID','Protein','top20btm10'])
    # for each sample
    for sID in ad.obs.ImageID.sort_values().unique():
        subAD = ad[ad.obs.ImageID == sID]
        for x in subAD.var_names:
            aX = subAD[:,x].X.flatten()
            # compute 20 largest values in aX
            top20 = np.sort(aX)[-20:]
            # compute the mean of bottom 10th percentile of aX
            btm10 = np.sort(aX)[:int(len(aX)*0.1)]
            #print(aX.shape, top20.shape, btm10.shape)
            top20btm10 = np.mean(top20)/np.mean(btm10)
            temp = pd.DataFrame({'ImageID':sID,'Protein':x,'top20btm10':top20btm10}, index = [0])
            top20btm10DF = pd.concat([top20btm10DF,temp],axis=0)
    # pandas remove append function so 
    #            top20btm10DF = top20btm10DF.append({'ImageID':sID,'Protein':x,'top20btm10':top20btm10}, ignore_index=True)
    return top20btm10DF

# for each marker clip value to mean of top 20 values
def mediantop20(subAD):
    outAD = subAD.copy()
    for ix,x in enumerate(outAD.var_names):
        aX = subAD[:,x].X.flatten()
        top20 = np.sort(aX)[-20:]
        outAD.X[:,ix] = np.clip(subAD[:,x].X.flatten(),0,np.median(top20))
    return outAD

# remove expression outliers from the data
def removeOutliers(ad):
    # separate each sample
    s = {}
    for sID in ad.obs.ImageID.sort_values().unique():
        print(sID)
        s[sID] = ad[ad.obs['ImageID']==sID]
        print(s[sID].X.max(axis=0))
        s[sID] = mediantop20(s[sID])
        print(s[sID].X.max(axis=0))
    outAD = sc.concat(s.values())
    return outAD 

def computeSpatialScore(subAD, p1,p2,p3):
    """
    # given three labels compute spatial score
    # as distance between closest centroids of p1 and p2 and p1 and p3
    # for each cell in p2, compute the distance between the closet cell in p1 and p3
    """
    centroidsP1 = subAD.obs[subAD.obs['phenotype']==p1][['spatial_X','spatial_Y']].values
    centroidsP3 = subAD.obs[subAD.obs['phenotype']==p3][['spatial_X','spatial_Y']].values
    # create a new DF for p2
    
    p2DF = pd.DataFrame(columns=['objID','distance','distP2P3'])

    for ix,row in tqdm.tqdm(subAD.obs.loc[subAD.obs['phenotype']==p2].iterrows(), total=subAD.obs.loc[subAD.obs['phenotype']==p2].shape[0]):
        c = row[['spatial_X','spatial_Y']].values.astype(float)
        distance = sp.spatial.distance.cdist([c],centroidsP1).min()
        distP2P3 = sp.spatial.distance.cdist([c],centroidsP3).min()
        p2DF = p2DF.append({'objID': ix,'distance':distance,'distP2P3':distP2P3},ignore_index=True)
    p2DF['spatialScore'] = p2DF['distance']/p2DF['distP2P3']
    return p2DF

def calculate_Density(phenotype_area_matrix : pd.DataFrame, phenotype_count_matrix :  pd.DataFrame, scale=1e-6) ->  pd.DataFrame:
    """
    # Calculate cell density [cell count]/[cell area] and Change scale µm^2 to mm^2
    # You can modify scale using scale parameter for multiply µm^2 by what you want
    """
    tissue_area = phenotype_area_matrix["Total Cells"]*scale
    density = phenotype_count_matrix.div(tissue_area,axis=0)
    density["Tissue Area"] = tissue_area
    columnd= list([density.columns[-1]]) + list(density.columns[:-1])
    return density[columnd]

def calculate_NearestNeighbor(subAD : ad.AnnData, phenotype_1 : str, phenotype_2 : str) -> Tuple[pd.DataFrame,pd.DataFrame]:
    """
    # as distance between closest centroids of p1 and p2 
    # for each cell in p2, compute the distance between the closet cell in p1 
    it has bring one and one match start to end. it is means that start object has only one end object even if start object has many end objects have same distance.
    Maybe in that case, this function raise error about array shape, so if you got a array shape error, we suggest you check that case
    please report error, for making better service 
    """
    # 여기에 단순히 출발 하는 오브젝트 말고 도착 오브젝트도 구할 수 있으면 좋겠어 그래야 나중에 그림 그릴 수 있을 거 같아

    static_df = pd.DataFrame(index=["Min","Mean","Median","Max","Std"])

    centroidsP1 = subAD.obs[subAD.obs['phenotype']==phenotype_1][['spatial_X','spatial_Y']].values
    centroidsP2 = subAD.obs[subAD.obs['phenotype']==phenotype_2][['spatial_X','spatial_Y']].values
    # create a new DF for p2
    
    p1p2DF = pd.DataFrame(columns=['Start_index',"Start_point","End_index",'End_point','distance',"Start -> End"])
    # P1에서 가장 가까운 P2의 거리
    print("Calculate distatance {} -> {}".format(phenotype_1,phenotype_2))
    for ix,row in tqdm.tqdm(subAD.obs.loc[subAD.obs['phenotype']==phenotype_1].iterrows(), total=subAD.obs.loc[subAD.obs['phenotype']==phenotype_1].shape[0]):
        c = row[['spatial_X','spatial_Y']].values.astype(float)
        distance_array= sp.spatial.distance.cdist([c],centroidsP2)
        distance = distance_array.min()
        nearest_end = np.where(distance_array == distance.min())[1]
        end_point= centroidsP2[nearest_end][0]
        end_index = subAD.obs[subAD.obs.spatial_X == end_point[0]][subAD.obs.spatial_Y == end_point[1]].index[0]
        p1p2DF = p1p2DF.append({"Start_index": ix, "Start_point" : c,"End_index":end_index, "End_point":end_point,'distance':distance,"Start -> End" :"{} -> {}".format(phenotype_1,phenotype_2)},ignore_index=True)
#    print("[ {} -> {} ] mean__: {}".format(phenotype_1, phenotype_2, p1p2DF.distance.mean()))
#    print("[ {} -> {} ] min___: {}".format(phenotype_1, phenotype_2, p1p2DF.distance.min()))
#    print("[ {} -> {} ] median: {}".format(phenotype_1, phenotype_2, p1p2DF.distance.median()))
#    print("[ {} -> {} ] max___: {}".format(phenotype_1, phenotype_2, p1p2DF.distance.max()))
#    print("[ {} -> {} ] std___: {}".format(phenotype_1, phenotype_2, p1p2DF.distance.std()))
    static_df["[ {} -> {} ]".format(phenotype_1, phenotype_2)]=[p1p2DF.distance.mean(),
                                            p1p2DF.distance.min(),
                                            p1p2DF.distance.median(),
                                            p1p2DF.distance.max(),
                                            p1p2DF.distance.std()]



    p2p1DF = pd.DataFrame(columns=['Start_index',"Start_point","End_index",'End_point','distance',"Start -> End"])
    # P2에서 가장 가까운 P1의 거리
    print("Calculate distatance {} -> {}".format(phenotype_2,phenotype_1))
    for ix,row in tqdm.tqdm(subAD.obs.loc[subAD.obs['phenotype']==phenotype_2].iterrows(), total=subAD.obs.loc[subAD.obs['phenotype']==phenotype_2].shape[0]):
        c = row[['spatial_X','spatial_Y']].values.astype(float)
        distance_array= sp.spatial.distance.cdist([c],centroidsP1)
        distance = distance_array.min()
        nearest_end = np.where(distance_array == distance.min())[1]
        end_point= centroidsP1[nearest_end][0]
        end_index = subAD.obs[subAD.obs.spatial_X == end_point[0]][subAD.obs.spatial_Y == end_point[1]].index[0]
        p2p1DF = p2p1DF.append({"Start_index": ix, "Start_point" : c,"End_index":end_index, "End_point":end_point,'distance':distance, "Start -> End": "{} -> {}".format(phenotype_2,phenotype_1)},ignore_index=True)
#    print("[ {} -> {} ] mean__: {}".format(phenotype_2, phenotype_1, p2p1DF.distance.mean()))
#    print("[ {} -> {} ] min___: {}".format(phenotype_2, phenotype_1, p2p1DF.distance.min()))
#    print("[ {} -> {} ] median: {}".format(phenotype_2, phenotype_1, p2p1DF.distance.median()))
#    print("[ {} -> {} ] max___: {}".format(phenotype_2, phenotype_1, p2p1DF.distance.max()))
#    print("[ {} -> {} ] std___: {}".format(phenotype_2, phenotype_1, p2p1DF.distance.std()))
    static_df["[ {} -> {} ]".format(phenotype_2, phenotype_1)]=[p2p1DF.distance.mean(),
                                            p2p1DF.distance.min(),
                                            p2p1DF.distance.median(),
                                            p2p1DF.distance.max(),
                                            p2p1DF.distance.std()]
    nearest_pair_df = pd.concat([p1p2DF,p2p1DF])
    return static_df, nearest_pair_df


def generate_cell_Table(data_groupby: object, type_string: str) -> pd.DataFrame:
    phenotype_matrix = data_groupby.index.to_frame(False)
    phenotype_matrix[type_string] = list(data_groupby)
    phenotype_matrix = phenotype_matrix.sort_values(by=['Target', 'phenotype'], ascending=True)

    unknown_matrix = phenotype_matrix[phenotype_matrix["phenotype"] == "Unknown"]
    phenotype_matrix = phenotype_matrix[phenotype_matrix["phenotype"] != "Unknown"]

    unknown_matrix = unknown_matrix[["Target", type_string]].rename(columns={"Target": "phenotype", type_string: "Unknown"}).set_index("phenotype")

    arranged_table = pd.DataFrame(index=phenotype_matrix.phenotype.unique())
    for id in phenotype_matrix["Target"].unique():
        temp = phenotype_matrix[phenotype_matrix.Target == id][["phenotype", type_string]].rename(columns={type_string: id}).set_index("phenotype")
        arranged_table = pd.merge(arranged_table, temp, left_index=True, right_index=True, how='left')

    arranged_table = arranged_table.T
    arranged_table = pd.merge(arranged_table, unknown_matrix, left_index=True, right_index=True, how='left').fillna(0)
    arranged_table["Total Cells"] = arranged_table.sum(axis=1)
    arranged_table.loc['Total'] = arranged_table.sum(axis=0)
    
    return arranged_table



def main():
    st.title("Phenotype Analysis")

    # 파일 업로드
    uploaded_file = st.file_uploader("Upload your H5AD file", type=["h5ad"])

    if uploaded_file is not None:
        # H5AD 파일 로드
        adata = read_h5ad(uploaded_file)
        st.success("File uploaded successfully.")

        # Outlier 제거
        cleaned_data = removeOutliers(adata)
        st.success("Data cleaned from outliers.")

        # Phenotype Counts 및 Phenotype Areas 계산 및 결과 출력
        df = cleaned_data.obs[['Target', 'phenotype', 'Area']]  # Area 열은 예시로 사용함, 필요에 따라 다른 열로 변경
        phenotype_counts = df.groupby('Target')['phenotype'].value_counts()
        phenotype_count_matrix = generate_cell_Table(phenotype_counts, "count")

        phenotype_area = df.groupby(['Target', 'phenotype'])['Area'].sum()
        phenotype_area_matrix = generate_cell_Table(phenotype_area, "area")

        st.subheader("Phenotype Counts:")
        st.dataframe(phenotype_count_matrix)

        st.subheader("Phenotype Areas:")
        st.dataframe(phenotype_area_matrix)

        # Cell Density 계산 및 결과 출력
        density_df = calculate_Density(phenotype_area_matrix, phenotype_count_matrix)
        st.subheader("Cell Density:")
        st.dataframe(density_df)

        # 사용자로부터 Phenotype 1 및 Phenotype 2를 선택받아 Nearest Neighbor Distances 계산하고 결과 출력
        phenotype_1 = st.selectbox("Select phenotype 1:", df['phenotype'].unique(), key="phenotype_1")
        phenotype_2 = st.selectbox("Select phenotype 2:", df['phenotype'].unique(), key="phenotype_2")
        p1p2DF, static_df = calculate_NearestNeighbor(cleaned_data, phenotype_1, phenotype_2)

        st.subheader("Nearest Neighbor Distances:")
        st.dataframe(p1p2DF)

        st.subheader("Static Information:")
        st.dataframe(static_df)

        # Static Information 데이터 테이블을 CSV 파일로 다운로드
        csv = static_df.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="static_info.csv">Download Static Information CSV File</a>'
        st.markdown(href, unsafe_allow_html=True)

        # Top 20th to Bottom 10th Ratio 계산 및 결과 출력
        top20btm10DF = computeTop20Btm10(cleaned_data)
        st.subheader("Top 20th to Bottom 10th Ratio:")
        st.dataframe(top20btm10DF)

        # Phenotype 1, Phenotype 2, Phenotype 3을 선택받아 Spatial Score를 계산하고 결과 출력
        p1 = st.selectbox("Select phenotype 1:", cleaned_data.obs['phenotype'].unique(), key="phenotype_3")
        p2 = st.selectbox("Select phenotype 2:", cleaned_data.obs['phenotype'].unique(), key="phenotype_4")
        p3 = st.selectbox("Select phenotype 3:", cleaned_data.obs['phenotype'].unique(), key="phenotype_5")
        spatial_score_df = computeSpatialScore(cleaned_data, p1, p2, p3)
        st.subheader("Spatial Score:")
        st.dataframe(spatial_score_df)

if __name__ == "__main__":
    main()