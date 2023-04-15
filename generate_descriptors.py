import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

def generate_descriptors(file_path,
                             smiles_column_name, disease_column_name):
        """
        :param file_path: path to the .csv table.
        :param smiles_column_name: string reflecting the SMILES column name.
        :return: a DataFrame object containing N columns with chemical descriptors from Rdkit
        + 1 column with disease categories.
        """
        Smiles_and_Disease = pd.read_csv(file_path)[[smiles_column_name, disease_column_name]]
        Smiles_and_Disease = Smiles_and_Disease.drop_duplicates().dropna().reset_index(drop=True)
        
        descriptors = [x[0] for x in Chem.Descriptors._descList]
        calculator = MolecularDescriptorCalculator(descriptors)

        dataset_with_descriptors = []
        
        for smile in Smiles_and_Disease[smiles_column_name]:
          mol = Chem.MolFromSmiles(smile)                             
          dataset_with_descriptors.append(calculator.CalcDescriptors(mol))
          
        dataset_with_descriptors = pd.DataFrame(dataset_with_descriptors)
        dataset_with_descriptors.columns = descriptors
        
        dataset_with_descriptors["Disease"] = Smiles_and_Disease[disease_column_name]
        dataset_with_descriptors = dataset_with_descriptors.dropna().reset_index(drop=True)                                              

        #dataset_with_descriptors = None
        
        return dataset_with_descriptors