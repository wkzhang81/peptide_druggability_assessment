# peptide_druggability_assessment
### Just a little tool for peptides screening!!!

### Copy the 'peptide_druggability_assessment.py' and 'batch.py' files to your folder and run the following command:
````
python3 batch.py
````

### Sample file was provided as 'test_peptides.csv' and the result was provided as 'test_peptides_results.csv'
<img width="3287" height="408" alt="image" src="https://github.com/user-attachments/assets/2a64a941-13b4-41e4-9afb-839e2a38dd17" />

### Basic notions:
- Solubility Assessment (Based on the Principles of the CamSol Method)
- Charge/pKa Analysis (Based on the Henderson-Hasselbalch Equation)
- Protein Binding Affinity Prediction (Based on Machine Learning Models)
- Molecular Descriptor Calculation (Using RDKit)
- Overall = (solubility * 0.3) + (permeability * 0.3) + (stability * 0.4)

### Requirements
- numpy>=1.21.0
- pandas>=1.3.0
- matplotlib>=3.4.0
- seaborn>=0.11.0
