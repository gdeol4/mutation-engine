# Local imports
import pandas as pd
from mutation_engine_3.settings import BASE_DIR
from shared.constants import BACKBONE, SUBSTRUCTURES
from .smiles_tools import sanitize_smiles

# Python built-in modules
import time

# Third-party libraries
import numpy as np
import rdkit 
from rdkit import Chem
from rdkit.Chem import MolToSmiles
import selfies
from selfies import encoder, decoder

import os

def generate_mutations(iterations, input_molecule, sample_size):
    unique_mutants_set = set()
    
    for i in range(iterations):
        start_time = time.time()

        # Randomize molecules
        rand_mols = randomize_molecules(input_molecule, sample_size=sample_size)
        if time.time() - start_time > 60:
            print(f"Iteration {i+1} took too long after randomization. Stopping the process.")
            break

        # Obtain mutations
        mutants = obtain_mutations(rand_mols)
        if time.time() - start_time > 60:
            print(f"Iteration {i+1} took too long after mutation. Stopping the process.")
            break

        # Update data
        new_mutants_count = len(unique_mutants_set)
        unique_mutants_set.update(mutants)
        new_mutants_count = len(unique_mutants_set) - new_mutants_count
        print(f"Iteration {i+1}: New mutants this iteration: {new_mutants_count}")

        if time.time() - start_time > 60:
            print(f"Iteration {i+1} took too long after fingerprinting. Stopping the process.")
            break
    
    # Convert the set of unique mutants to a DataFrame with a 'smiles' column
    unique_mutants_df = pd.DataFrame(list(unique_mutants_set), columns=['smiles'])
    
    return unique_mutants_df

def randomize_molecules(input_molecule, sample_size):
    '''Randomize the molecules represented by the input molecule using SELFIES.
    
    Parameters:
    - input_molecule (string): The input molecule as a SELFIE string
    - sample_size (int): Number of random samples to generate
    
    Returns:
    - selfies_ls (list): A list of SELFIES representing the randomized molecules
    '''
    starting_smi = BACKBONE[input_molecule]
    num_random_samples = sample_size
    
    mol = Chem.MolFromSmiles(starting_smi)
    if mol is None:
        raise Exception('Invalid starting structure encountered')
    
    randomized_smile_orderings = [randomize_smiles(mol) for _ in range(num_random_samples)]
    
    # Convert all the molecules to SELFIES
    selfies_ls = [encoder(x) for x in randomized_smile_orderings]
    
    return selfies_ls

def randomize_smiles(mol):
    '''
    Returns a random (dearomatized) SMILES given an RDKit mol object of a molecule.

    Parameters:
    mol (rdkit.Chem.rdchem.Mol): RDKit mol object (None if invalid smile string smi)

    Returns:
    mol (rdkit.Chem.rdchem.Mol): RDKit mol object (None if invalid smile string smi)
    '''
    if not mol:
        return None

    Chem.Kekulize(mol)
    return MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False, kekuleSmiles=True)

def obtain_mutations(selfies_ls):
    '''Obtain mutations of the SELFIES within the given list.
    
    Parameters:
    - selfies_ls (list): A list of SELFIES
    
    Returns:
    - all_smiles_collect (list): A list of mutated SMILES strings
    '''
    num_mutation_ls = [1, 2, 3, 4, 5]
    all_smiles_collect = []
    
    for num_mutations in num_mutation_ls:
        # Mutate the SELFIES:
        selfies_mut = get_mutated_SELFIES(selfies_ls.copy(), num_mutations=num_mutations)

        # Convert back to SMILES:
        smiles_back = [decoder(x) for x in selfies_mut]
        all_smiles_collect.extend(smiles_back)
    
    return all_smiles_collect

def get_mutated_SELFIES(selfies_ls, num_mutations): 
    ''' Mutate all the SELFIES in 'selfies_ls' 'num_mutations' number of times. 
    
    Parameters:
    - selfies_ls (list): A list of SELFIES 
    - num_mutations (int): Number of mutations to perform on each SELFIES within 'selfies_ls'
    
    Returns:
    - selfies_ls (list): A list of mutated SELFIES
    '''
    for _ in range(num_mutations): 
        selfie_ls_mut_ls = []
        for str_ in selfies_ls: 
            str_chars = get_selfie_chars(str_)
            max_molecules_len = len(str_chars) + num_mutations
            
            selfie_mutated, _ = mutate_selfie(str_, max_molecules_len, index=None)
            selfie_ls_mut_ls.append(selfie_mutated)
        
        selfies_ls = selfie_ls_mut_ls.copy()
    return selfies_ls

def get_selfie_chars(selfie):
    '''Obtain a list of all selfie characters in the selfie string
    
    Parameters: 
    - selfie (string): A selfie string - representing a molecule 
    
    Example: 
    >>> get_selfie_chars('[C][=C][C][=C][C][=C][Ring1][Branch1_1]')
    ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[Branch1_1]']
    
    Returns:
    - chars_selfie (list): List of selfie characters present in the molecule selfie
    '''
    chars_selfie = []  # A list of all SELFIE symbols from the selfie string
    while selfie != '':
        chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
        selfie = selfie[selfie.find(']')+1:]
    return chars_selfie

def mutate_selfie(selfie, max_molecules_len, write_fail_cases=True, index=None):
    '''Return a mutated selfie string (only one mutation on selfie is performed)
    
    Parameters:
    - selfie (string): SELFIE string to be mutated 
    - max_molecules_len (int): Mutations of SELFIE string are allowed up to this length
    - write_fail_cases (bool): If true, failed mutations are recorded in "selfie_failure_cases.txt"
    - index (int or None): Index value of the character to mutate. If None, a random index is selected.
    
    Returns:
    - selfie_mutated (string): Mutated SELFIE string or original SELFIE if valid mutation isn't found.
    - smiles_canon (string): Canonical SMILES of mutated SELFIE string or original SELFIE's SMILES.
    '''
    valid = False
    fail_counter = 0
    MAX_ATTEMPTS = 1000  # Maximum number of mutation attempts before giving up
    
    chars_selfie = get_selfie_chars(selfie)
    
    while not valid and fail_counter < MAX_ATTEMPTS:
        fail_counter += 1
                
        alphabet = ['[OH1-1]', '[C][OH0]'] # optional prenol group '[O][C][\\C][=C][Branch1][C][/C][C]'
        
        random_choice = 2  # simplified to account for substitution only  
        
        # Insert a character at a random location or choose an index value 
        if random_choice == 2:
            if index is None:
                random_index = np.random.randint(len(chars_selfie))
            else:
                random_index = index
                
            random_character = np.random.choice(alphabet, size=1)[0]
            
            if random_index == 0:
                selfie_mutated_chars = [random_character] + chars_selfie[random_index+1:]
            else:
                selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index+1:]              
        else: 
            raise Exception('Invalid operation trying to be performed')

        selfie_mutated = "".join(x for x in selfie_mutated_chars)
        sf = "".join(x for x in chars_selfie)
        
        try:
            smiles = selfies.decoder(selfie_mutated)
            mol, smiles_canon, done = sanitize_smiles(smiles)
            if len(selfie_mutated_chars) > max_molecules_len or smiles_canon == "" or substructure_preserver(mol) == False:
                done = False
            if done:
                valid = True
            else:
                valid = False
        except:
            valid = False
            if fail_counter > 1 and write_fail_cases == True:
                with open("selfie_failure_cases.txt", "a+") as f:
                    f.write('Tried to mutate SELFIE: ' + str(sf) + ' To Obtain: ' + str(selfie_mutated) + '\n')
    
    # If a valid mutation is not found after MAX_ATTEMPTS, return the original selfie and its corresponding smiles.
    if not valid:
        print(f"Warning: Could not find a valid mutation for SELFIE: {selfie} after {MAX_ATTEMPTS} attempts.")
        return (selfie, selfies.decoder(selfie))
    
    return (selfie_mutated, smiles_canon)

def substructure_preserver(mol):
    '''
    Check for substructure violations
    Return True: contains a substructure violation
    Return False: No substructure violation
    '''        
    if mol.HasSubstructMatch(rdkit.Chem.MolFromSmarts(SUBSTRUCTURES[0])) is True:
        return True  # The molecule has a substructure violation
    elif mol.HasSubstructMatch(rdkit.Chem.MolFromSmarts(SUBSTRUCTURES[1])) is True:
        return True  # The molecule has a substructure violation
    else: 
        return False  # Molecule does not have a substructure violation