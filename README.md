# Beta-1 Adrenergic Receptor QSAR Model
This project aims to implement a QSAR pipeline where drugs pIC50 scores are predicted against the human beta-1 adrenergic receptor, a G protein-coupled receptor located on the surface of the heart, with a primary purpose of increasing cardiac activity when activated. The data used in this project was acquired from the ChEMBL database, where the molecular structure of the compounds were encoded as Morgan fingerprints which were used to train a Random Forest regression model to predict the pIC50 scores of untested antagonist compounds. This project is useful for saving time and resources when deciding what drug is worth testing in a wet lab.

## Background
The beta-1 adrenergic receptor is one of four adrenergic receptors with a purpose of increasing cardiac output when activated. This GPCR is essential as it serves as a critical mediator of the body's autonomic regulation, particularly under stress or changing physiological demands. Prolonged activation of this receptor can lead to desensitization and downregulation as well as heart failure. Unlike other adrenergic receptors, Beta-1 ARs makes up 75-80% of cardiac beta-AR density as its purpose is to control cardiac contractility. Clinically, beta-1 adrenergic receptors can be targeted with beta blockers to reduce excessive sympathetic stimulation. Beta blockers may be used to treat hypertension, by helping lower blood pressure, as well as assist with performance anxiety as stress-related increases in heart rate becomes controlled. Therefore, computational modeling provides a cost effective, efficient method for researching beta-1 adrenergic receptor pharmacology and aiding in the development of improved therapeutic agents.

## Objectives
-	The data used in this project was acquired from ChEMBL database, where it was processed and filtered for the pIC50 scores of various drugs and their molecular SMILES string acting on the human beta-1 adrenergic receptor. 
-	This project is a QSAR pipeline demonstrating the process through which I acquired the data, cleaned it, processed it, and ran it through a Random Forest Regression model to correlate molecular fingerprints to their pIC50 score.
-	The model was validated using an 80/20 train/test split and 5-fold cross validation in order to test its accuracy in predicting pIC50 scores.

## Pipeline Overview
```
ChEMBL Database
┌───────────────────────────────────────────────────────────────────────────────────────────┐
|- Bioactivity data for human beta-1 adrenergic receptor (CHEMBL213) was retrieved using the| ChEMBL web resource client.                                                                 |
└───────────────────────────────────────────────────────────────────────────────────────────┘
      |
      ▼
Data Preprocessing
┌───────────────────────────────────────────────────────────────────────────────────────────────┐
|- The data acquired was filtered to binding assays only, removing compounds with missing values| and IC50 ceiling values >= 10,000 nM.                                                           |
|- IC50 values were converted to pIC50.                                                         |  
|- Final data set contained 450 compounds.                                                      |  
└───────────────────────────────────────────────────────────────────────────────────────────────┘
      |
      ▼
Morgan Fingerprint Generation
┌───────────────────────────────────────────────────────────────────────────────────────────────┐
|- SMILES strings for the filtered compounds were converted into RDKit molecule objects.        |
|- Morgan fingerprints were generated using a radius of 2 and 2048 bits.                        | 
|- Converted chemical structures into a fixed length of binary vectors which acted as the independent variable for the model.                                                             |
└─────────────────────────────────────────────────────────────────────────────────────────────────┘
      |
      ▼
Random Forest Model
┌─────────────────────────────────────────────────────────────────────────────────────────────────┐
|- The forest consisted of 200 trees with a max depth of 7 and a seed of 123 for reproducibility. |
|- The random forest model was used to learn relationships from the Morgan fingerprints to predict pIC50 binding values.                                                                             |
└─────────────────────────────────────────────────────────────────────────────────────────────────┘
      |
      ▼
pIC50 Prediction
┌─────────────────────────────────────────────────────────────────────────────────────────────────┐
|- Model output predicted pIC50 values of untested data.                                          |
|- Used to quantify binding affinity of compounds against the human beta-1 adrenergic receptor for practical prioritization.                                                                         |
└─────────────────────────────────────────────────────────────────────────────────────────────────┘
```
## Results

|              Metric                         | Value |
|---------------------------------------------|-------|
| Training R^2                                | 0.900 |
| Test R^2                                    | 0.622 |
| Cross-Validated Mean R^2 (shuffled, 5-fold) | 0.715 |
| Cross-Validated Std R^2                     | 0.060 |

For full results, discussion, and visualizations, see ['beta1_qsar_exploratory.ipynb']

## Requirements

Install the required libraries using the following commands:
```bash
pip install pandas numpy scikit-learn matplotlib chembl-webresource-client
conda install -c conda-forge rdkit
```

|        Library            |         Purpose                        |
|---------------------------|----------------------------------------|
| pandas                    | Data manipulation and filtering        |
| numpy                     | Numerical computation                  |
| scikit-learn              | Random Forest model and evaluation     |
| matplotlib                | Data visualization                     |
| rdkit                     | SMILES to Morgan fingerprint conversion|
| chembl-webresource-client | ChEMBL database querying               |

## How to Run

1. Clone the repository:
  ```bash
   git clone https://github.com/joshua1759/beta1-adrenergic-qsar-model.git
  ```
2. Install requirements (See Requirements section)
3. Open Jupyter Notebook (beta1_qsar_exploratory.ipynb
4. Run the Notebook from top to bottom. (If ChEMBL is unavailable, the pipeline will automatically fall back to the locally saved dataset or the GitHub Hosted CSV.)

## References

1. Alhayek S, Preuss CV. Beta 1 Receptors. [Updated 2023 Aug 14]. In: StatPearls [Internet]. Treasure Island (FL): StatPearls Publishing; 2026 Jan-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK532904/
2. Srimanee, N., Chinda, K., & Chatsudthipong, V. (2024). Regulation of beta-adrenergic receptors in the heart: A review on emerging therapeutic strategies for heart failure. Cells, 13(20), 1674. https://doi.org/10.3390/cells13201674
3. Yu BH, Kang EH, Ziegler MG, Mills PJ, Dimsdale JE. Mood states, sympathetic activity, and in vivo beta-adrenergic receptor function in a normal population. Depress Anxiety. 2008;25(7):559-64. doi: 10.1002/da.20338. PMID: 17583588; PMCID: PMC2680308.

## Author 
Joshua Pologne
https://github.com/joshua1759
