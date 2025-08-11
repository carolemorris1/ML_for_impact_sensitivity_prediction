# ML_for_impact_sensitivity_prediction
Repository containing energetic materials for impact sensitivity database and ML models

Information on included files:
 
- Feature extraction:
	- '00-featureextractor.py' - this takes in two files: 'Dataset' and 'Descriptor List'
	- 'Dataset' contains SMILES strings for input into model
	- 'Descriptor List' contains all of the model descriptors to be called by the script to be extracted for each molecule
	
Output is 'DataOutput'. NOTE: if adding in new data (which should go into TestDataset.xlsx, and the file renamed in the featureextractor script), this should be checked by the user before moving on to the next step, particularly for the Class II trigger bonds (see paper for description of these). The script is not perfect at picking these out.
 
- Models:
	- '01-DataPreprocessing' Jupyter notebook - this takes in the output of the feature extraction script ('DataOutput') and scales all data to be of similar magnitudes. Output is 'df_processed'
	- '02-ModelFitting' Jupyter notebook - takes in 'df_processed' and 'targets' (a file containing IS classification values) and allows for all models to be run. Also takes in parameters in the form of .pkl files for faster notebook running if desired.
 	- '03-Predict_IS' Jupyter notebook - takes in DataOutput (output of 00-featureextractor.py for a TestDataset file; preprocessing is applied directly in this notebook) and makes predictions on new SMILES strings.

Instructions for use:

1. Run 00-featureextractor.py. Ensure that the code is running in the environment in which RDKit is installed. Ensure that the input file name is as appropriate for intended use.

2. Run the Jupyter Notebook '01-DataPreprocessing'. 

3. Train and test any models in the Jupyter Notebook '02-ModelFitting'.

4. Optionally, run '03-Predict_IS' notebook for any new molecules to be tested.
