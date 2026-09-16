import streamlit as st
import pandas as pd
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import DataStructs

# Load the saved model
with open('beta1_model.pkl', 'rb') as f:
    model = pickle.load(f)

# App title
st.title('Beta-1 Adrenergic Receptor pIC50 Predictor')
st.write('Enter a SMILES string to predict its pIC50 binding affinity against the human beta-1 adrenergic receptor.')

# Text area for multiple SMILES (one per line)
smiles_input = st.text_area('Enter SMILES strings (one per line):')

if st.button('Predict and Rank'):
    # Split input into list of SMILES
    smiles_list = smiles_input.strip().split('\n')
    
    valid_smiles = []
    predictions = []

    for smiles in smiles_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            continue
        generator = rdFingerprintGenerator.GetMorganGenerator(radius = 2, fpSize = 2048) 
        fingerprint = generator.GetFingerprint(molecule)
        arr = np.zeros(2048)
        DataStructs.ConvertToNumpyArray(fingerprint, arr)
        prediction = model.predict(arr.reshape(1, -1))[0]
        valid_smiles.append(smiles)
        predictions.append(prediction)
    # After loop - create ranked table
    results = pd.DataFrame({
        'SMILES': valid_smiles,
        'Predicted pIC50': predictions
    })
    results = results.sort_values('Predicted pIC50', ascending=False)
    results['Rank'] = range(1, len(results) + 1)
    st.dataframe(results)