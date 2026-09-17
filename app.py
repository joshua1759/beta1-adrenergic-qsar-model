import streamlit as st
import pandas as pd
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import DataStructs
from chembl_webresource_client.new_client import new_client
import os
from sklearn.ensemble import RandomForestRegressor                                         

@st.cache_resource
def load_and_train_model():
    with st.status("Initializing model...", expanded = True) as status:
        try:
            st.write("Acquiring Data from ChEMBL...")

            activities = new_client.activity.filter(                 
                target_chembl_id = "CHEMBL213",                       
                target_type = "SINGLE PROTEIN",                        
                organism = "Homo sapiens",
                standard_type = "IC50"
            )
            df = pd.DataFrame(activities)                             
            df.to_csv('data/beta1_raw.csv', index = False)
            st.write("✓ ChEMBL data acquired sucessfully")
        except:
            st.write("✗ Failed to acquire ChEMBL data")
            if os.path.exists('data/beta1_raw.csv'):
                st.write("Accessing local backup file...")
                df = pd.read_csv('data/beta1_raw.csv')
                st.write("✓ Local file loaded sucessfully")
            else:
                #TODO: add GitHubraw URL after first commit
                st.write("✗ Local file not found")
                st.write("Accessing github backup")
                df = pd.read_csv('https://raw.githubusercontent.com/joshua1759/beta1-adrenergic-qsar-model/refs/heads/master/data/beta1_raw.csv')
                st.write("✓ GitHub file loaded successfully")
        
        st.write("Preprocessing data...")
#=================================================================================================================
        df = df[df['assay_type'] == 'B']    
        selected_columns = ['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'standard_units']
        df2 = df[selected_columns]
        df3 = df2.dropna()
        conversion = pd.to_numeric(df3['standard_value'], errors='coerce')
        df3['standard_value'] = conversion
        df3 = df3.dropna()
        df3 = df3[df3['standard_value'] < 10000]
        df3 = df3[df3['standard_value'] > 0]
        df3 = df3.copy()
        df3['pIC50'] = -np.log10(df3['standard_value'] * 1e-9)
        df_clean = df3.drop(columns=['standard_value','standard_units'])
#=================================================================================================================
        st.write("✓ Data preprocessed successfully")

        st.write("Generating molecular fingerprints...")
#=================================================================================================================
        def smiles_to_fingerprint(canonical_smiles):
            molecule = Chem.MolFromSmiles(canonical_smiles)                                  
            if molecule is None:
                return None
            generator = rdFingerprintGenerator.GetMorganGenerator(radius = 2, fpSize = 2048) 
            fingerprint = generator.GetFingerprint(molecule)                                 
            return fingerprint

        df_clean['fingerprint'] = df_clean['canonical_smiles'].apply(smiles_to_fingerprint)
        df_clean = df_clean[df_clean['fingerprint'].notna()]
        fingerprint_list = []
        for fp in df_clean['fingerprint']:
            arr = np.zeros(2048) 
            DataStructs.ConvertToNumpyArray(fp, arr) 
            fingerprint_list.append(arr)

        X = np.array(fingerprint_list)
        y = df_clean['pIC50'].values
#=================================================================================================================
        st.write("✓ Molecular fingerprints generated successfully")

        st.write("Training model...")
#=================================================================================================================
        model = RandomForestRegressor(n_estimators = 200, max_depth = 7, random_state = 123)
        model.fit(X, y)
        with open('data/beta1_model.pkl', 'wb') as f:
            pickle.dump(model, f)
#=================================================================================================================
        st.write("✓ Model trained successfully")

        status.update(label = "Model ready!", state = "complete")

    return model

model = load_and_train_model()

st.title('Beta-1 Adrenergic Receptor pIC50 Predictor')
st.write('Enter a SMILES string to predict its pIC50 binding affinity against the human beta-1 adrenergic receptor.')

smiles_input = st.text_area('Enter smiles strings (one per line)')

st.caption('Note: Higher pIC50 indicates greater predicted binding potency against beta-1 adrenergic receptor.')

if st.button('Predict'):
    if smiles_input:
        smiles_list = smiles_input.splitlines()
        valid_smiles = []
        fingerprints = []
        for smiles in smiles_list:
            fp = Chem.MolFromSmiles(smiles)
            if fp is not None:
                generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
                fingerprint = generator.GetFingerprint(fp)
                arr = np.zeros(2048)
                DataStructs.ConvertToNumpyArray(fingerprint, arr)
                fingerprints.append(arr)
                valid_smiles.append(smiles)
            else:
                st.write(f"Invalid SMILES: {smiles}")

        if fingerprints:
            X_new = np.array(fingerprints)
            predictions = model.predict(X_new)
            results_df = pd.DataFrame({'SMILES': valid_smiles, 'Predicted pIC50': predictions})
            results_df = results_df.sort_values('Predicted pIC50', ascending=False)
            results_df['Rank'] = range(1, len(results_df) + 1)
            st.write(results_df)
        else:
            st.write("No valid SMILES strings provided.")
    else:
        st.write("Please enter at least one SMILES string.")