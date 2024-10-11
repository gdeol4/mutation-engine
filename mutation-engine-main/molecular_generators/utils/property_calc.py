import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

class MolPropertyCalc:
    """
    A class for calculating molecular properties given a DataFrame with a 'smiles' column.
    """

    def __init__(self, df):
        self.df = df

    def calculate_properties(self):
        properties = []
        existing_columns = set(self.df.columns)

        for _, row in self.df.iterrows():
            smiles = row['smiles']
            mol = Chem.MolFromSmiles(smiles)
            prop_dict = {'smiles': smiles, 'inchikey': None}  # Initialize inchikey as None

            if mol is not None:
                prop_dict['inchikey'] = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
                if 'molecular_weight' not in existing_columns:
                    prop_dict['molecular_weight'] = Descriptors.ExactMolWt(mol)
                if 'nhet' not in existing_columns:
                    prop_dict['nhet'] = Lipinski.NumHeteroatoms(mol)
                if 'nrot' not in existing_columns:
                    prop_dict['nrot'] = Lipinski.NumRotatableBonds(mol)
                if 'nring' not in existing_columns:
                    prop_dict['nring'] = Chem.rdMolDescriptors.CalcNumRings(mol)
                if 'nha' not in existing_columns:
                    prop_dict['nha'] = Lipinski.NumHAcceptors(mol)
                if 'nhd' not in existing_columns:
                    prop_dict['nhd'] = Lipinski.NumHDonors(mol)
                if 'logp' not in existing_columns:
                    prop_dict['logp'] = Descriptors.MolLogP(mol)
            else:
                for col in ['inchikey', 'molecular_weight', 'nhet', 'nrot', 'nring', 'nha', 'nhd', 'logp']:
                    if col not in existing_columns:
                        prop_dict[col] = None

            for col in self.df.columns:
                if col != 'smiles':
                    prop_dict[col] = row[col]

            properties.append(prop_dict)

        return pd.DataFrame(properties)