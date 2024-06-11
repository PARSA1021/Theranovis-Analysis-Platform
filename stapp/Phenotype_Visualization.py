import streamlit as st
import scanpy as sc
import tempfile
import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.cm import get_cmap

def load_data(uploaded_file):
    """Load AnnData object from an uploaded file."""
    adata = sc.read_h5ad(uploaded_file)
    return adata

def draw_cell_phenotype(newAD, target, cMap_file, show=False, transparent=True):
    """
    Plot Cell Phenotype using matplotlib.
    """
    # Load the color map from the JSON file
    with open(cMap_file, 'r') as f:
        cMap = json.load(f)
    
    # Check if the target exists in the 'Target' column
    if target not in newAD.obs['Target'].unique():
        print(f"The target '{target}' does not exist in the 'Target' column of the AnnData object.")
        return False
    
    # Filter the AnnData object to include only the target phenotype and not 'Unknown'
    ss = newAD[(newAD.obs['Target'] == target) & (newAD.obs['phenotype'] != 'Unknown')]
    
    # Check if there are any cells after filtering
    if ss.n_obs == 0:
        print(f"No cells with the phenotype '{target}' and not marked as 'Unknown' were found.")
        return False
    
    # Get unique phenotypes and sort them
    phenotype = ss.obs['phenotype'].unique().tolist()
    phenotype.sort()
    
    # Create a palette for the unique phenotypes that are in the color map
    unique_phenotypes_in_cMap = [pheno for pheno in phenotype if pheno in cMap]
    palette = [cMap[i] for i in unique_phenotypes_in_cMap]
    
    # Plot the cells using scanpy's embedding function
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))  # Increase the width to make the aspect ratio rectangular
    ax.invert_yaxis()
    ax.axis('off')
    sc.pl.embedding(ss, basis='X_spatial', color='phenotype', ax=ax, show=show, palette=palette, frameon=False)
    
    # Adjust the subplot layout
    fig.subplots_adjust(left=0, top=1, bottom=0, right=1)
    
    # Close the figure to free up memory
    plt.close(fig)
    
    return fig

def main():
    st.set_option('deprecation.showPyplotGlobalUse', False)  # Disable PyplotGlobalUseWarning

    st.title("Cell Phenotype visualization")

    uploaded_file = st.file_uploader("Upload h5ad file", type="h5ad")
    uploaded_cmap = st.file_uploader("Upload color map json file", type="json")

    if uploaded_file and uploaded_cmap:
        try:
            # Load data
            adata = load_data(uploaded_file)
            if adata is None:
                raise ValueError("Failed to load AnnData.")
            
            # Save the uploaded color map JSON file to a temporary location
            with tempfile.NamedTemporaryFile(delete=False, suffix='.json') as temp_cmap_file:
                temp_cmap_file.write(uploaded_cmap.read())
                cmap_file_path = temp_cmap_file.name
            
            # Load color map from the temporary file
            with open(cmap_file_path, 'r') as f:
                cmap_data = json.load(f)
            c_map = cmap_data.get('cMap', {})
            if not c_map:
                raise ValueError("Color map is empty. Please check the uploaded file.")

            # Display cell phenotype for each target
            for target in adata.obs['Target'].unique():
                st.subheader(f"Visualization for Target: {target}")
                fig_phenotype = draw_cell_phenotype(adata, target, cmap_file_path)
                if fig_phenotype is not None:
                    st.pyplot(fig_phenotype)
        
        except Exception as e:
            st.error(f"An error occurred: {e}")

# Run the Streamlit app
if __name__ == "__main__":
    main()