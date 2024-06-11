import streamlit as st
import pandas as pd
import json
import base64

# 파일을 읽어서 DataFrame으로 변환하는 함수를 정의합니다.
def load_data(file):
    if file.name.endswith('.csv'):  # 파일 확장자에 따라 구분자를 결정합니다.
        separator = ","
    else:  
        separator = "\t"  # 탭으로 구분된 파일을 가정합니다.
    data = pd.read_csv(file, sep=separator)  
    return data

# Phenotype 정보가 포함된 열을 결정하는 함수를 정의합니다.
def determine_phenotype_column(data_phenotype):
    potential_columns = ['Phenotype', 'Cell Type', 'Cell_Type']  # Phenotype 정보가 있을 것으로 예상되는 열 이름입니다.
    for column in potential_columns:
        if column in data_phenotype.columns:
            return column
    return None

# 메인 함수를 정의합니다.
def main():
    st.title(' Text, CSV File Transform to Json ')  # 웹 애플리케이션의 타이틀을 설정합니다.

    # Annotation, Detection, Phenotype 파일을 업로드 받습니다.
    uploaded_file_annotation = st.file_uploader("Annotation 파일을 선택하세요.", type=['txt'])
    uploaded_file_detection = st.file_uploader("Detection 파일을 선택하세요.", type=['txt'])
    uploaded_file_phenotype = st.file_uploader("Phenotype.csv 파일을 선택하세요.", type=['csv'])

    # 각 파일이 업로드 되었는지 확인하고, 데이터를 처리합니다.
    data = {}

    if uploaded_file_annotation is not None:
        data_annotation = load_data(uploaded_file_annotation)
        data['input'] = {"Annotations": uploaded_file_annotation.name}
        st.write("### Annotation Table", data_annotation)

    if uploaded_file_detection is not None:
        data_detection = load_data(uploaded_file_detection)
        data['input']['Detection'] = uploaded_file_detection.name
        st.write("### Detection Table", data_detection)
        
        # 사용자가 원하는 열을 선택할 수 있도록 합니다.
        if not data_detection.empty:
            selected_columns = st.multiselect("Detection 데이터에서 확인하고 싶은 열을 선택하세요.", data_detection.columns)
            st.write("### 선택한 열의 데이터", data_detection[selected_columns])

            custom_names = {}
            for column in selected_columns:
                custom_names[column] = st.text_input(f"Enter custom name for '{column}'", key=f'custom_name_{column}')
            data['newNames'] = custom_names

    if uploaded_file_phenotype is not None:
        data_phenotype = load_data(uploaded_file_phenotype)
        data['input']['Phenotype'] = uploaded_file_phenotype.name
        st.write("### Phenotype Table", data_phenotype)

        phenotype_column = determine_phenotype_column(data_phenotype)

        if phenotype_column is not None:
            phenotype_options = data_phenotype[phenotype_column].tolist()
            selected_phenotypes_for_cMap = st.multiselect("Select phenotypes for cMap", options=phenotype_options + ['Unknown'], key='cMap_selection')
            selected_phenotypes_for_nn = st.multiselect("Select phenotypes for Nearest Neighbor", options=phenotype_options, key='nn_selection')

            cMap_info = {}
            if selected_phenotypes_for_cMap:
                for phenotype in selected_phenotypes_for_cMap:
                    if phenotype == 'Unknown':
                        info = st.color_picker(f"Select color for 'Unknown' in cMap", key='cMap_color_unknown')
                        cMap_info['Unknown'] = info
                        continue

                    row = data_phenotype[data_phenotype[phenotype_column] == phenotype].index  
                    if not row.empty:
                        row = row[0]
                        info = st.color_picker(f"Select color for '{phenotype}' in cMap", key=f'cMap_color_{phenotype}')
                        cMap_info[phenotype] = info
                    else:
                        cMap_info[phenotype] = ""

                data['cMap'] = cMap_info

            nn_info = {}
            if selected_phenotypes_for_nn:
                for i, phenotype in enumerate(selected_phenotypes_for_nn, start=1):
                    row = data_phenotype[data_phenotype[phenotype_column] == phenotype].index  
                    if not row.empty:
                        row = row[0]
                        row_data = data_phenotype.loc[row].to_dict()
                        nn_info[f"Phenotype {i}"] = phenotype

                data['Nearest Neighbor'] = nn_info

            if data:
                auto = st.checkbox("Auto", key='auto_checkbox')
                data['Auto'] = str(auto)

                with open('output.json', 'w') as json_file:
                    json.dump(data, json_file, indent=4)

                with open('output.json', 'r') as file:
                    json_data = file.read()
                    b64 = base64.b64encode(json_data.encode()).decode()
                    filename = "output.json"
                    href = f'<a href="data:file/json;base64,{b64}" download="{filename}" style="text-decoration:none; color:white; background-color:#008CBA; padding:10px 20px; border-radius:5px;">Download JSON File</a>'
                    st.markdown(href, unsafe_allow_html=True)

                st.write("### JSON File Contents:")
                with open('output.json', 'r') as json_file:
                    try:
                        json_content = json.load(json_file)
                        st.json(json_content)
                    except json.JSONDecodeError as e:
                        st.error(f"Error decoding JSON file: {e}")

        else:
            st.error("Phenotype 정보가 포함된 열을 찾을 수 없습니다.")

if __name__ == "__main__":
    main()
