# Python built-in modules

# Third-party modules
import re
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolFromSmiles, MolToSmiles
import pandas as pd
import datamol as dm
from rdkit import Chem
from rdkit import DataStructs

from mutation_engine_3.settings import BASE_DIR

from shared.constants import SUBSTRUCTURES
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

# Local imports
import os

def is_valid_inchikey(inchikey):
    """
    Checks if the given string is a valid InChIKey.
    InChIKeys have the following properties:
    - Length of 27 characters
    - Consists of uppercase letters and dashes (-)
    - Has two dashes separating it into three parts
    """
    pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-NP-Z]$'
    return bool(re.match(pattern, inchikey))


def sanitize_smiles(smi):
    '''
    Returns a canonical smile representation of smi.
    Parameters:
    smi (string): smile string to be canonicalized

    Returns:
    mol (rdkit.Chem.rdchem.Mol): RDKit mol object (None if invalid smile string smi)
    smi_canon (string): Canonicalized smile representation of smi (None if invalid smile string smi)
    conversion_successful (bool): True/False to indicate if conversion was successful
    '''
    try:
        mol = MolFromSmiles(smi, sanitize=True)
        smi_canon = MolToSmiles(mol, isomericSmiles=False, canonical=True)
        return (mol, smi_canon, True)
    except Exception as e:
        print(f"Error in sanitize_smiles: {e}")
        return (None, None, False)
    

def filter_by_substructure(df):
    # Use the first element in SUBSTRUCTURES as the reference SMILES for flavone
    reference_smiles = SUBSTRUCTURES[0]
    
    # Convert reference SMILES to a Mol object
    reference_mol = Chem.MolFromSmiles(reference_smiles)
    if not reference_mol:
        print(f"Could not parse the reference SMILES: {reference_smiles}")
        return None

    try:
        # Convert the Mol object to a SMARTS string and create a pattern for substructure matching
        reference_pattern = Chem.MolFromSmarts(Chem.MolToSmarts(reference_mol))
    except:
        print("Failed to convert to SMARTS. The molecule may have issues with kekulization.")
        return None

    # Use a lambda function to check for substructure and filter the DataFrame
    has_substructure = lambda x: Chem.MolFromSmiles(x) and Chem.MolFromSmiles(x).HasSubstructMatch(reference_pattern)
    
    return df[df['smiles'].apply(has_substructure)]



def filter_unstable_molecules(df):
    """Filter out unstable molecules based on specific criteria."""
    stable_df = df[df['smiles'].apply(is_stable)]
    return stable_df

def is_stable(smiles):
    """Determine if a molecule is stable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    # Example criteria: check for the absence of radicals
    num_radical_electrons = Descriptors.NumRadicalElectrons(mol)
    if num_radical_electrons > 0:
        return False
    return True

def filter_unusual_oxygen_bindings(df):
    """Filter out molecules with unusual oxygen bindings, including three consecutive oxygens."""
    def has_unusual_oxygen_bindings(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        # Check for consecutive oxygen atoms
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'O':
                return True
        # Check for three consecutive oxygens in the SMILES string
        if "OOO" in smiles:
            return True
        return False

    # Apply the filter to the DataFrame
    filtered_df = df[~df['smiles'].apply(has_unusual_oxygen_bindings)]
    return filtered_df


def process_generated_smiles(unique_mutants_df):
    canonical_smiles_set = set()

    # Loop through each SMILES string and apply transformations
    for smiles in unique_mutants_df['smiles']:
        mol = dm.to_mol(smiles, ordered=True)

        if not mol:
            continue

        mol = dm.fix_mol(mol)
        sanitized_mol = dm.sanitize_mol(mol, sanifix=True, charge_neutral=False)
        standardized_mol = dm.standardize_mol(sanitized_mol, disconnect_metals=False, reionize=True, uncharge=True, stereo=True)

        canonical_smiles = dm.to_smiles(standardized_mol, canonical=True)
        canonical_smiles_set.add(canonical_smiles)

    # Convert the set to a DataFrame
    canonical_df = pd.DataFrame(list(canonical_smiles_set), columns=['smiles'])

    # Apply the function to get only those molecules with a flavone backbone
    substructure_filtered_df = filter_by_substructure(canonical_df)

    # Apply the function to remove unstable molecules
    stability_filtered_df = filter_unstable_molecules(substructure_filtered_df)

    # Apply the function to remove molecules with unusual oxygen bindings
    filtered_df = filter_unusual_oxygen_bindings(stability_filtered_df)

    print(f'Count of filtered molecules: {len(filtered_df)}')

    return filtered_df