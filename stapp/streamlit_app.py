# Import necessary libraries
import streamlit as st
from PIL import Image

# Import modules
import make_json
import json_process_to_h5ad
import Phenotype_Analysis
import Phenotype_Visualization
import Nearest_Neighbor  # Import the new Nearest_Neighbor module

# Define a function to customize the sidebar and set the main content title
def main():
    # Set page configuration
    st.set_page_config(page_title="Theranovis", page_icon="stapp/logo.png")

    # Customize sidebar
    customize_sidebar()

    # Customize main content
    customize_main()

# Define a function to customize the sidebar
def customize_sidebar():
    # Create a container for the sidebar content
    sidebar_container = st.sidebar.container()

    # Load custom logo image
    logo_image = Image.open("logo.png")
    
    # Display logo image with custom styling
    sidebar_container.image(logo_image, use_column_width=True)

    # Customize the sidebar title
    sidebar_container.title("🔬 Theranovis Analysis Platform")
    
    # Company introduction
    sidebar_container.subheader("Welcome to Theranovis!")
    sidebar_container.markdown(
        """
        <div style='text-align: justify;'>
        Theranovis envisions 'Diagnosis as a New Drug,' providing accurate and cost-effective diagnostics and monitoring for various diseases.
        With highly sensitive diagnostic devices using fluorescence imaging and disease-specific tailored diagnostic methods, we can quantify protein biomarkers with high sensitivity.
        In addition, we have developed on-site diagnostic platforms, automation platforms for use in diagnostic medicine, and a digital Multiplex IHC pathology platform, enabling highly sensitive comprehensive analysis and digital transformation of chronic diseases, cancer patients, dementia, myocardial infarction, pulmonary edema, tuberculosis, latent tuberculosis, infectious diseases, and rare diseases.
        Theranovis will continue to develop precision medical in vitro diagnostic devices and solutions for early diagnosis, treatment direction, and prognosis management to evolve into a global diagnostic company.
        </div>
        """,
        unsafe_allow_html=True
    )
    sidebar_container.markdown("---")  # Add a separator

    # Radio buttons for page selection
    page_options = {
        "Make JSON": {"icon": "📄", "description": "Convert data to JSON format."},
        "JSON Process to H5AD": {"icon": "🔄", "description": "Process JSON to H5AD format."},
        "Phenotype Analysis": {"icon": "🔬", "description": "Perform phenotype analysis."},
        "Phenotype_Visualization": {"icon": "📊", "description": "Visualize phenotype data."},
        "Nearest Neighbor Analysis": {"icon": "🔍", "description": "Analyze nearest neighbors between phenotypes."}  # Added Nearest Neighbor Analysis
    }
    
    page = sidebar_container.radio("Select a page", list(page_options.keys()), format_func=lambda x: f"{page_options[x]['icon']} {x}")

    # Navigate to the selected page
    if page in page_options:
        sidebar_container.write(page_options[page]['description'])
        sidebar_container.markdown("---")  # Add a separator

        if page == "Make JSON":
            make_json.main()
        elif page == "JSON Process to H5AD":
            json_process_to_h5ad.main()
        elif page == "Phenotype Analysis":
            Phenotype_Analysis.main()
        elif page == "Phenotype_Visualization":
            Phenotype_Visualization.main()
        elif page == "Nearest Neighbor Analysis":
            Nearest_Neighbor.main()  # Call the Nearest Neighbor Analysis page

# Define a function to customize the main content
def customize_main():
    # Set main content title
    st.title("🧬 Multiplex IHC Analysis Platform")

    # Customize sidebar and main content style
    st.markdown(
        """
        <style>
        /* Sidebar style */
        .sidebar .sidebar-content {
            background-color: #F5F5F5;
            color: #333333;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            border-radius: 10px;
            box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1);
            padding: 10px;
        }
        .sidebar .sidebar-content .sidebar-section .sidebar-title {
            color: #008CBA;
            font-weight: bold;
        }
        .sidebar .sidebar-content .sidebar-section .sidebar-item {
            color: #666666;
        }
        .sidebar .sidebar-content .sidebar-section .sidebar-item.selected {
            background-color: #E6E6E6;
        }
        .sidebar .sidebar-content::-webkit-scrollbar {
            width: 8px;
        }
        .sidebar .sidebar-content::-webkit-scrollbar-thumb {
            background-color: #CCCCCC;
            border-radius: 10px;
        }

        /* Main content style */
        .reportview-container {
            background-color: #FFFFFF;
            color: #333333;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            border-radius: 10px;
            box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1);
            padding: 10px;
        }
        .reportview-container .main .block-container h1 {
            color: #008CBA;
        }
        .reportview-container .main .block-container p {
            color: #666666;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

# Entry point of the program
if __name__ == "__main__":
    main()
