import os
import numpy as np
import pandas as pd
import xgboost as xgb
from rdkit import Chem
import deepchem as dc

from molecular_generators.models import GeneratedFlavonoid
from predictors.models import AdmetProperties
from shared.constants import ADMET_PROPERTIES

ADMET_MODELS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'admet_models')

def admet_featurizer(mol_list):
    """
    Featurizes the given molecules and returns the features as a numpy array.
    Parameters:
    - mol_list (list): List of RDKit Mol objects
    Returns:
    - fp_val (numpy.ndarray): Featurized representation of the molecules
    """
    maccskeys = dc.feat.MACCSKeysFingerprint()
    circular = dc.feat.CircularFingerprint()
    mol2vec = dc.feat.Mol2VecFingerprint()
    mordred = dc.feat.MordredDescriptors(ignore_3D=True)
    rdkit = dc.feat.RDKitDescriptors()

    maccskeys_mol = maccskeys.featurize(mol_list)
    circular_mol = circular.featurize(mol_list)
    mol2vec_mol = mol2vec.featurize(mol_list)
    mordred_mol = mordred.featurize(mol_list)
    rdkit_mol = rdkit.featurize(mol_list)

    fp_val = np.column_stack((maccskeys_mol,
                              circular_mol,
                              mol2vec_mol,
                              mordred_mol,
                              rdkit_mol))

    fp_val = np.nan_to_num(fp_val, nan=0, posinf=0)

    return fp_val


def pad_features(features, expected_length):
    if features.shape[1] < expected_length:
        placeholder_features = np.zeros((features.shape[0], expected_length - features.shape[1]))
        features = np.hstack((features, placeholder_features))
    else:
        features = features[:, :expected_length]
    return features


def predict_admet_properties(features, inchikey_list, smiles_list):
    """
    Predict ADMET properties for the molecules using XGBoost models.
    
    Parameters:
    - features (numpy.ndarray): Numpy array containing the featurized representation of molecules.
    - inchikey_list (list): List of InChIKeys corresponding to the molecules.
    - smiles_list (list): List of SMILES strings corresponding to the molecules.
    
    Returns:
    - pd.DataFrame: A DataFrame containing predicted ADMET properties for each molecule.
    """
    admet_df = pd.DataFrame()
    admet_df['inchikey'] = inchikey_list
    admet_df['smiles'] = smiles_list

    features = pad_features(features, expected_length=4336)

    for model in ADMET_PROPERTIES:
        model_path = os.path.join(ADMET_MODELS_DIR, f"{model}.json")
        if os.path.exists(model_path):
            x = xgb.Booster()
            x.load_model(model_path)
            x.set_param({"predictor": "cpu_predictor"})
            array = x.inplace_predict(features)
            admet_df[model] = array
            print(f"Finished predicting for {model} model and added to dataframe.")
        else:
            print(f"Model {model} not found in {ADMET_MODELS_DIR}. Skipping predictions for this model.")

    return admet_df


def featurize_and_predict_admet(input_df):
    print("Featurizing molecules...")
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in input_df['smiles']]
    featurized_data = admet_featurizer(mol_list)
    print("Featurization completed.")
    
    inchikey_list = input_df['inchikey'].tolist()
    smiles_list = input_df['smiles'].tolist()
    
    print("Predicting ADMET properties...")
    admet_df = predict_admet_properties(featurized_data, inchikey_list, smiles_list)
    print("ADMET prediction completed.")
    
    return admet_df


# def get_inchikeys_to_predict():
#     # Get all InChIKeys from the generated_flavonoids table
#     generated_inchikeys = GeneratedFlavonoid.objects.values_list('inchikey', flat=True)

#     # Get all InChIKeys from the admet_properties table
#     existing_inchikeys = AdmetProperties.objects.values_list('inchikey', flat=True)

#     # Find the InChIKeys that are in generated_flavonoids but not in admet_properties
#     inchikeys_to_predict = generated_inchikeys.exclude(inchikey__in=existing_inchikeys)

#     return list(inchikeys_to_predict)

def get_inchikeys_to_predict(table1, table2, column):
    # Get all values from the specified column in table1
    table1_inchikeys = table1.objects.values_list(column, flat=True)
    
    # Get all values from the specified column in table2
    table2_inchikeys = table2.objects.values_list(column, flat=True)
    
    # Find the values that are in table1 but not in table2
    inchikeys_to_predict = table1_inchikeys.exclude(**{f'{column}__in': table2_inchikeys})
    
    return list(inchikeys_to_predict)


def upload_admet_properties(admet_df):
    # Upload the predicted ADMET properties to the database
    for _, row in admet_df.iterrows():
        admet_properties = AdmetProperties(
            inchikey_id=row['inchikey'],
            ames=row['ames'],
            bbb=row['bbb'],
            bioavailability=row['bioavailability'],
            caco2=row['caco2'],
            clearance_hepatocyte=row['clearance_hepatocyte'],
            clearance_microsome=row['clearance_microsome'],
            cyp2c9=row['cyp2c9'],
            cyp2c9_substrate=row['cyp2c9_substrate'],
            cyp2d6=row['cyp2d6'],
            cyp2d6_substrate=row['cyp2d6_substrate'],
            cyp3a4=row['cyp3a4'],
            cyp3a4_substrate=row['cyp3a4_substrate'],
            dili=row['dili'],
            half_life=row['half_life'],
            herg=row['herg'],
            hia=row['hia'],
            ld50=row['ld50'],
            lipophilicity=row['lipophilicity'],
            pgp=row['pgp'],
            ppbr=row['ppbr'],
            solubility=row['solubility'],
            vdss=row['vdss']
        )
        admet_properties.save()
    print("ADMET properties uploaded to the database.")


def predict_and_upload_admet_properties(request):
    # Get the InChIKeys to predict
    inchikeys_to_predict = get_inchikeys_to_predict(GeneratedFlavonoid, AdmetProperties, 'inchikey')
    print(f"Number of InChIKeys to predict: {len(inchikeys_to_predict)}")

    # Get the molecules corresponding to the InChIKeys
    molecules_to_predict = GeneratedFlavonoid.objects.filter(inchikey__in=inchikeys_to_predict)
    print(f"Number of molecules to predict: {molecules_to_predict.count()}")

    # Create a DataFrame from the molecules
    input_df = pd.DataFrame(list(molecules_to_predict.values('inchikey', 'smiles')))
    print("Input DataFrame created.")

    if input_df.empty:
        # If the input DataFrame is empty, create an empty DataFrame with the expected columns
        admet_df = pd.DataFrame(columns=['inchikey', 'smiles', 'caco2_permeability', 'hia_absorption', 'pgp_inhibition', 'cyp_inhibition', 'cyp_substrate', 'half_life', 'clearance', 'bbb_permeability', 'renal_ocr', 'ld50'])
    else:
        # Predict ADMET properties
        admet_df = featurize_and_predict_admet(input_df)

    # Upload the predicted ADMET properties to the database
    upload_admet_properties(admet_df)

    # Fetch the first 10 rows from the admet_properties table
    rows = AdmetProperties.objects.all()[:20]

    # Get the number of rows for each table
    generated_flavonoids_count = GeneratedFlavonoid.objects.count()
    admet_properties_count = AdmetProperties.objects.count()

    print(f"Number of rows in GeneratedFlavonoid table: {generated_flavonoids_count}")
    print(f"Number of rows in AdmetProperties table: {admet_properties_count}")

    #return render(request, 'admet_results.html', {'rows': rows})
    
    # return render(request, 'admet_results.html', {'rows': rows})