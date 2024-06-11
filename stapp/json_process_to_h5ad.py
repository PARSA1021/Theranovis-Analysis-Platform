import pandas as pd
import os
import json
import anndata as ad
import scimap as sm
import scanpy as sc
import streamlit as st

def load_data(uploaded_file):
    try:
        data = sc.read_h5ad(uploaded_file)
        return data
    except Exception as e:
        st.error(f"Error occurred while loading data: {e}")
        return None

def find_leaf(annot_df: pd.DataFrame, level_list: list) -> list:
    children = []
    for i in level_list:
        children = children + list(annot_df.loc[annot_df["Parent"] == i]["Name"])
    if len(children) == 0:
        return level_list
    return find_leaf(annot_df, children)

def read_data(qupath_annotation_path: str, qupath_detection_path: str) -> tuple:
    annot_df = pd.read_csv(qupath_annotation_path, sep='\t')
    annot_df.reset_index(drop=True, inplace=True)

    detect_df = pd.read_csv(qupath_detection_path, sep='\t')
    detect_df.reset_index(drop=True, inplace=True)

    sections = annot_df.loc[annot_df["Parent"] == "Root object (Image)"]["Name"]
    annotations = find_leaf(annot_df, sections)
    
    if set(detect_df.Parent.unique()) != set(annotations):
        raise ValueError("Annotation and Detection files do not match.")

    return annot_df, detect_df, annotations

def preprocess_detection_data(detect_df: pd.DataFrame, new_names: dict, output_path: str):
    protein_df = detect_df[[col for col in detect_df.columns if 'mean' in col and 'Cell:' in col]]
    if not protein_df.empty:
        protein_df.rename(columns=new_names, inplace=True)

        spatial_columns = ['Centroid X µm', 'Centroid Y µm', 'Cell: Area','Image', 'Parent']
        spatial_names = {'Centroid X µm': 'spatial_X', 'Centroid Y µm': 'spatial_Y', 'Cell: Area': 'Area', 'Image': 'ImageID', 'Parent': 'Target'}
        spatial_df = detect_df[spatial_columns]
        spatial_df.rename(columns=spatial_names, inplace=True)
        spatial_df['ImageID'] = spatial_df['ImageID'].apply(lambda x: x.split('.')[0])

        adata = ad.AnnData(protein_df)
        adata.obs = spatial_df

        adata.write_h5ad(output_path)
    else:
        raise ValueError("No valid protein columns found.")

def visualize_anndata(adata: ad.AnnData):
    st.subheader("Variables data:")
    st.dataframe(adata.var)

    st.subheader("Observations data:")
    st.dataframe(adata.obs)

def save_anndata(adata: ad.AnnData, output_path: str, filename: str):
    adata.write(output_path)
    with open(output_path, "rb") as f:
        data = f.read()
    st.download_button(label="Download Anndata Object", data=data, file_name=filename, mime="application/octet-stream")

def add_phenotype_info(phenotype_file_name: str, adata: ad.AnnData) -> ad.AnnData:
    try:
        adata = sm.pp.rescale(adata, imageid='ImageID', method='by_image')
        phenoDF = pd.read_csv(phenotype_file_name)
        sm.tl.phenotype_cells(adata, phenoDF, label='phenotype')
        adata.obsm['X_spatial'] = adata.obs[['spatial_X', 'spatial_Y']].values
        adata.uns['gates'] = adata.obs.columns.tolist()
        adata.obsm_keys = ['X_spatial']
        return adata
    except Exception as e:
        st.error(f"Error occurred: {e}")
        return None

def main():
    st.title("Data Processing and Visualization")

    custom_data_file_name = st.file_uploader("Upload JSON file.", type=['json'])
    if custom_data_file_name is not None:
        custom_content = custom_data_file_name.getvalue()
        custom = json.loads(custom_content)

        data_dir = "Bio_Data"
        annotations_file_name = os.path.join(data_dir, custom["input"]["Annotations"])
        detection_file_name = os.path.join(data_dir, custom["input"]["Detection"])
        phenotype_file_name = os.path.join(data_dir, custom["input"]["Phenotype"])

        try:
            annotDF, detectDF, annotations = read_data(annotations_file_name, detection_file_name)
            preprocessed_ann_data = os.path.join(data_dir, "pre_ann_data.h5ad")
            preprocess_detection_data(detectDF, custom["newNames"], preprocessed_ann_data)

            st.success("Data processing completed successfully.")

            adata = ad.read_h5ad(preprocessed_ann_data)
            visualize_anndata(adata)

            filename = st.text_input("Enter file name for download:", "pre_ann_data.h5ad")
            if st.button("Download Anndata Object"):
                save_anndata(adata, preprocessed_ann_data, filename)

        except Exception as e:
            st.error(f"Error occurred: {e}")

    uploaded_h5ad = st.file_uploader("Upload h5ad file.", type=["h5ad"])
    if uploaded_h5ad is not None:
        adata = load_data(uploaded_h5ad)
        if adata is not None:
            uploaded_csv = st.file_uploader("Upload CSV file.", type=["csv"])
            if st.button("Add Phenotype Information"):
                if uploaded_csv is not None:
                    adata = add_phenotype_info(uploaded_csv, adata)
                    if adata is not None:
                        phenotype_filename = st.text_input("Enter file name for phenotype data:", "pre_phenotype_data.h5ad")
                        save_anndata(adata, "phenotype_ann_data.h5ad", phenotype_filename)
                        st.success("Phenotype information added successfully.")
                    else:
                        st.error("An error occurred while adding phenotype information.")
                else:
                    st.warning("Please upload a CSV file to add phenotype information.")
            st.write(adata)
            st.subheader("Observations information")
            st.write(adata.obs)
            st.subheader("Variables information")
            st.write(adata.var)
            st.subheader("obsm information")
            for key, value in adata.obsm.items():
                st.write(f"{key}: {value}")
            st.subheader("uns information")
            for key, value in adata.uns.items():
                st.write(f"{key}: {value}")

if __name__ == "__main__":
    main()