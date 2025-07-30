# ML_for_impact_sensitivity_prediction
Repository containing energetic materials for impact sensitivity database and ML models

Information on included files:

- 'ESI Data Set' - all compounds (including SMILES strings, structures & IS data with reference) for all 485 compounds used to train and test the models.
 
- Feature extraction folder: 
	- 'FeatureExtractor.py' - this takes in two files: 'Dataset' and 'Descriptor List'
	- 'Dataset' contains SMILES strings for input into model
    - 'additionaltestset' is a spreadsheet which should be added to (i.e. add compound 31 onwards) in the event a user wishes to add their own compounds to be tested only. Do not delete the included compounds 1-30, add below. These are included in order for features of any added compounds to be correctly scaled.
	- 'Descriptor List' contains all of the model descriptors to be called by the script to be extracted for each molecule
	Output is 'DataOutput'. NOTE: if adding in new data, this should be checked by the user before moving on to the next step, particularly for the Class II trigger bonds (see paper for description of these). The script is not 100% accurate at picking these out as sometimes SMILES strings aren't structured in the expected way.
 
- Notebooks folder:
	- 'DataPreprocessing' Jupyter notebook - this takes in the output of the feature extraction script ('DataOutput') and scales all data to be of similar magnitudes. Output is 'df_minmax'
	- 'ModelsRun' Jupyter notebook - takes in 'df_minmax' and 'targets' (a file containing IS classification values) and allows for all models to be run.

** If the user wishes to check the classification of their own compound (with known SMILES), use FeatureExtractor.py first as usual, but change the filenames on lines 148 (to additionaltestset) and 161 (to DataOutputTest).

Instructions for use:

1. Run FeatureExtractor.py. Ensure that the code is running in the environment in which RDKit is installed. Ensure that the input file name is as appropriate for intended use.

2. Run the Jupyter Notebook 'DataPreprocessing'. 

3. Train and test any models in the Jupyter Notebook 'ModelsRun'.
