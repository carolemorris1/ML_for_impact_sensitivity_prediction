# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:55:11 2024

@author: heatherquayle

Adapted from codes by Jack M. Hemingway (University of Edinburgh, In-person communication) 
and J. Rein et al. (ref: doi.org/10.1002/anie.202218213)

"""
### script to extract molecular descriptors from sheet containing SMILES ###

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import rdkit
from rdkit import Chem
import math
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import GraphDescriptors as gdesc
from rdkit.Chem import Lipinski
from rdkit.Chem import EState
from rdkit.Chem import MolSurf
from collections import OrderedDict



def counter(molecule, fg):
    molecule_mol = Chem.MolFromSmiles(molecule)
    funcgrp = Chem.MolFromSmarts(fg)
    output = (len(molecule_mol.GetSubstructMatches(funcgrp)))
    return output

def chemdescriptor(molecule):
    # calculate molecular formula weight from SMILES string
    MolWt = Descriptors.MolWt(Chem.MolFromSmiles(molecule))
    HeavyMolWt = Descriptors.HeavyAtomMolWt(Chem.MolFromSmiles(molecule))
    numH = math.floor(MolWt-HeavyMolWt)
    # deduce atom count for the molecule
    numC = counter(molecule, '[#6]')
    numN = counter(molecule, '[#7]')
    numO = counter(molecule, '[#8]')
    # calculate oxygen balance for the molecule
    OxBal = (-1600/MolWt)*(2*numC+(numH/2)-numO)
    # count smarts needed for h-bond ratio calculation
    no2 = counter(molecule, '[$([NX3](=O)=O),$([NX3+](=O)[O-])]')
    nh2 = counter(molecule, '[NH2]')
    co = counter(molecule, '[#6X3]=[OX1]')
    coc = counter(molecule, '[$([#8])D2]([#6])[#6]')
    con = counter(molecule, '[$([#8])D2]([#7])[#6]')
    nhr2 = counter(molecule, '[#7X3H1]')
    ringn = counter(molecule, '[#7D2H0]')
    oh = counter(molecule, '[$(O)H]')
    nr3 = counter(molecule, '[$([#7]);X3;!$(N[O]);H0]')
    # calculate h-bond ratio
    donors = (2*nh2) + nhr2 + oh
    acceptors = (2*no2) + nh2 + co + coc + con + nhr2 + ringn + oh + nr3
    hbondratio = (donors**2)/acceptors
    # define other features not from smarts queries
    rotbonds = Lipinski.NumRotatableBonds(Chem.MolFromSmiles(molecule))
    heteroatoms = Lipinski.NumHeteroatoms(Chem.MolFromSmiles(molecule))
    totrings = Lipinski.RingCount(Chem.MolFromSmiles(molecule))
    aromrings = Lipinski.NumAromaticRings(Chem.MolFromSmiles(molecule))
    alirings = Lipinski.NumAliphaticRings(Chem.MolFromSmiles(molecule))
    VSA8 = EState.EState_VSA.VSA_EState8(Chem.MolFromSmiles(molecule))
    VSA5 = MolSurf.SMR_VSA5(Chem.MolFromSmiles(molecule))
    TPSA = MolSurf.TPSA(Chem.MolFromSmiles(molecule))
    # extract values needed for KMF calculation & calculate
    K1 = gdesc.Kappa1(Chem.MolFromSmiles(molecule))
    K2 = gdesc.Kappa2(Chem.MolFromSmiles(molecule))
    heavy_atom_count=Descriptors.HeavyAtomCount(Chem.MolFromSmiles(molecule))
    kmf=K1*K2/heavy_atom_count
    # count smarts for trigger bonds
    class1 = no2
    cnali = counter(molecule, '[$([#7]-[#6]);!$(N[O])]') 
    coali = counter(molecule, '[C;!R]-[OH0]')
    calinaro = counter(molecule, '[C;!R]-[#7;R]')
    nnh2 = counter(molecule, '[#7]~[NH2]')
    cp = counter(molecule, '[#6]-[#15]')
    class2 = coali + calinaro + nnh2 + cp + cnali
    # function can return values listed below
    return (rotbonds, heteroatoms, totrings, aromrings, alirings, VSA8, VSA5, TPSA, MolWt, numC, numN, numO, OxBal, hbondratio, kmf, class1, class2)


  # this function runs each molecule in molecules_inp through the entire list of descriptors descriptors_df.
def analyse(molecules_inp, descriptors_inp):
    # define empty list to store outputs.
    allmolecules = []
    # iterate through each molecule in molecules_inp dataframe.
    for molecule in range(len(molecules_inp)):
        # define temp list to store descriptor values during iteration through molecules.
        temp = []
        currentmolecule = molecules_inp.iloc[molecule]['SMILES']
        # calculate values for descriptors that cannot be calculated from SMARTS queries.
        (rotbonds, heteroatoms, totrings, aromrings, alirings, VSA8, VSA5, TPSA, MolWt, numC, numN, numO, OxBal, hbondratio, kmf, class1, class2) = chemdescriptor(currentmolecule)
        # iterate through each descriptor in descriptors_inp dataframe.
        for descriptor in range(len(descriptors_inp)):
            currentdescriptor = descriptors_inp.iloc[descriptor]['SMARTS']
            # these conditional branches deal with descriptors that cannot be calculated from SMARTS queries.
            if pd.isnull(descriptors_inp.loc[descriptor]['SMARTS']):
                # deal with 'special case' descriptors: 
                if descriptors_inp.loc[descriptor]['Feature'] == 'Rot Bonds':
                    temp.append(rotbonds)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Num Heteroatoms':
                    temp.append(heteroatoms)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Total Rings':
                    temp.append(totrings)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Aromatic Rings':
                    temp.append(aromrings)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Aliphatic Rings':
                    temp.append(alirings)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'VSAEState8':
                    temp.append(VSA8)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'SMRVSA5':
                    temp.append(VSA5)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'TPSA':
                    temp.append(TPSA)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'MW':
                    temp.append(MolWt)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Number of C':
                    temp.append(numC)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Number of N':
                    temp.append(numN)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Number of O':
                    temp.append(numO)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Oxygen Balance':
                    temp.append(OxBal)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'H Bond Ratio':
                    temp.append(hbondratio)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'KMF':
                    temp.append(kmf)
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Class I TBs':
                    temp.append(class1)  
                elif descriptors_inp.loc[descriptor]['Feature'] == 'Class II TBs':
                    temp.append(class2)
            # call counter to deal with 'normal' descriptors using SMARTS queries.
            else:
                smartscount = counter(currentmolecule,currentdescriptor)
                temp.append(smartscount)
        # append values from molecule to list containing all molecules.
        allmolecules.append(temp)
    return allmolecules

# # # MAIN JOB # # #

# load 'Dataset.xlsx' as a pandas dataframe molecules_inp
molecules_inp = pd.read_excel('c:/Projects/Chemistry/Files for GitHub/PredictTest/TestDataset.xlsx')
# load 'Descriptor List.xlsx' as a pandas dataframe descriptors_inp
descriptors_inp = pd.read_excel('c:/Projects/Chemistry/Files for GitHub/PredictTest/Descriptor List.xlsx')
# analyse molecules and obtain counting descriptor list of lists.
listofmolecules = analyse(molecules_inp, descriptors_inp)


# create column titles for output workbook.
columns = list((descriptors_inp.T).iloc[0,:])
#concatenate column titles, molecules_df, and moleculesWithDescriptors_df.
moleculeswithdescriptors_out = pd.DataFrame(listofmolecules, columns=columns)

# write to output Excel workbook (in row order of input SMILES)
moleculeswithdescriptors_out.to_excel('c:/Projects/Chemistry/Files for GitHub/PredictTest/DataOutput.xlsx')