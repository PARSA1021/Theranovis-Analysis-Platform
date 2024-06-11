import streamlit as st
import anndata
import pandas as pd
import json
import tempfile
import os

from Utils import calculate_NearestNeighbor, draw_phenotype_nearest_neighbor, read_data

def main():
    st.title("Phenotype Nearest Neighbor Analysis")

    st.sidebar.header("Upload your files")

    uploaded_h5ad_file = st.sidebar.file_uploader("Upload your .h5ad file", type=["h5ad"])
    uploaded_json_file = st.sidebar.file_uploader("Upload your config JSON file", type=["json"])
    uploaded_annotation_file = st.sidebar.file_uploader("Upload your annotation file", type=["txt"])
    uploaded_detection_file = st.sidebar.file_uploader("Upload your detection file", type=["txt"])

    if uploaded_h5ad_file and uploaded_json_file and uploaded_annotation_file and uploaded_detection_file:
        with tempfile.NamedTemporaryFile(delete=False) as tmp_h5ad, \
            tempfile.NamedTemporaryFile(delete=False) as tmp_json, \
            tempfile.NamedTemporaryFile(delete=False) as tmp_annotation, \
            tempfile.NamedTemporaryFile(delete=False) as tmp_detection:

            tmp_h5ad.write(uploaded_h5ad_file.read())
            tmp_json.write(uploaded_json_file.read())
            tmp_annotation.write(uploaded_annotation_file.read())
            tmp_detection.write(uploaded_detection_file.read())

            tmp_h5ad_path = tmp_h5ad.name
            tmp_json_path = tmp_json.name
            tmp_annotation_path = tmp_annotation.name
            tmp_detection_path = tmp_detection.name

        # JSON 파일 로드
        with open(tmp_json_path, "r") as json_file:
            config_data = json.load(json_file)

        # .h5ad 파일에서 데이터 로드
        newAD = anndata.read_h5ad(tmp_h5ad_path)

        # 필요한 정보 추출
        phenotype_1 = config_data["Nearest Neighbor"]["Phenotype 1"]
        phenotype_2 = config_data["Nearest Neighbor"]["Phenotype 2"]
        custom_cMap = config_data["cMap"]
        output_dir = tempfile.mkdtemp()

        # 주석 처리 코드
        annot_df, detect_df, annotations = read_data(tmp_annotation_path, tmp_detection_path)

        distance_matrix_dict = {}

        # newAD.obs의 각 고유한 대상에 대해 반복
        for target in newAD.obs['Target'].unique():
            subsetAD = newAD[newAD.obs['Target'] == target]
            subset_distance_matrix, subset_pair_matrix = calculate_NearestNeighbor(subsetAD, phenotype_1, phenotype_2)
            
            if subset_pair_matrix.empty:
                st.warning(f"No valid data found for target: {target}")
                continue

            draw_phenotype_nearest_neighbor(subsetAD, subset_pair_matrix, f"{phenotype_1} -> {phenotype_2}", custom_cMap, output_dir, f"{target}_to_{phenotype_1}")
            draw_phenotype_nearest_neighbor(subsetAD, subset_pair_matrix, f"{phenotype_2} -> {phenotype_1}", custom_cMap, output_dir, f"{target}_to_{phenotype_2}")
            
            distance_matrix_dict[f"{target}"] = subset_distance_matrix

        st.success("Analysis completed!")

        for target in newAD.obs['Target'].unique():
            if os.path.exists(os.path.join(output_dir, f"{target}_to_{phenotype_1}.png")):
                st.image(os.path.join(output_dir, f"{target}_to_{phenotype_1}.png"), caption=f"{target} to {phenotype_1}")
            if os.path.exists(os.path.join(output_dir, f"{target}_to_{phenotype_2}.png")):
                st.image(os.path.join(output_dir, f"{target}_to_{phenotype_2}.png"), caption=f"{target} to {phenotype_2}")

    else:
        st.warning("Please upload all required files.")

if __name__ == "__main__":
    main()
