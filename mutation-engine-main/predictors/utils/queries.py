# # from django.db.models import OuterRef, Subquery
# # from molecular_generators.models import GeneratedFlavonoid
# # from predictors.models import AdmetProperties
# # from predictors.utils.admet_tools import featurize_and_predict_admet
# # import pandas as pd

# # def get_inchikeys_to_predict():
# #     # Get all InChIKeys from the generated_flavonoids table
# #     generated_inchikeys = GeneratedFlavonoid.objects.values_list('smiles', flat=True)

# #     # Get all InChIKeys from the admet_properties table
# #     existing_inchikeys = AdmetProperties.objects.values_list('inchikey', flat=True)

# #     # Find the InChIKeys that are in generated_flavonoids but not in admet_properties
# #     inchikeys_to_predict = generated_inchikeys.exclude(inchikey__in=existing_inchikeys)

# #     return list(inchikeys_to_predict)

# # def predict_and_upload_admet_properties():
# #     # Get the InChIKeys to predict
# #     inchikeys_to_predict = get_inchikeys_to_predict()

# #     # Get the molecules corresponding to the InChIKeys
# #     molecules_to_predict = GeneratedFlavonoid.objects.filter(inchikey__in=inchikeys_to_predict)

# #     # Create a DataFrame from the molecules
# #     input_df = pd.DataFrame(list(molecules_to_predict.values('inchikey', 'smiles')))

# #     # Featurize and predict ADMET properties
# #     admet_df = featurize_and_predict_admet(input_df)

# #     # Upload the predicted ADMET properties to the database
# #     for _, row in admet_df.iterrows():
# #         admet_properties = AdmetProperties(
# #             inchikey_id=row['inchikey'],
# #             ames=row['ames'],
# #             bbb=row['bbb'],
# #             bioavailability=row['bioavailability'],
# #             caco2=row['caco2'],
# #             clearance_hepatocyte=row['clearance_hepatocyte'],
# #             clearance_microsome=row['clearance_microsome'],
# #             cyp2c9=row['cyp2c9'],
# #             cyp2c9_substrate=row['cyp2c9_substrate'],
# #             cyp2d6=row['cyp2d6'],
# #             cyp2d6_substrate=row['cyp2d6_substrate'],
# #             cyp3a4=row['cyp3a4'],
# #             cyp3a4_substrate=row['cyp3a4_substrate'],
# #             dili=row['dili'],
# #             half_life=row['half_life'],
# #             herg=row['herg'],
# #             hia=row['hia'],
# #             ld50=row['ld50'],
# #             lipophilicity=row['lipophilicity'],
# #             pgp=row['pgp'],
# #             ppbr=row['ppbr'],
# #             solubility=row['solubility'],
# #             vdss=row['vdss']
# #         )
# #         admet_properties.save()
# #------------------------------------------------------------------------------------------
# # def compare_table_row_counts():
# #     generated_flavonoids_count = GeneratedFlavonoid.objects.count()
# #     admet_properties_count = AdmetProperties.objects.count()

# #     if generated_flavonoids_count == admet_properties_count:
# #         return True
# #     else:
# #         return False

# # from django.db.models import OuterRef, Subquery
# # from molecular_generators.models import GeneratedFlavonoid
# # from predictors.models import AdmetProperties
# # from predictors.utils.admet_tools import featurize_and_predict_admet
# # import pandas as pd



# # def predict_and_upload_admet_properties():
# #     # Get the InChIKeys to predict
# #     inchikeys_to_predict = get_inchikeys_to_predict()
# #     print(f"Number of InChIKeys to predict: {len(inchikeys_to_predict)}")

# #     # Get the molecules corresponding to the InChIKeys
# #     molecules_to_predict = GeneratedFlavonoid.objects.filter(inchikey__in=inchikeys_to_predict)
# #     print(f"Number of molecules to predict: {molecules_to_predict.count()}")

# #     # Create a DataFrame from the molecules
# #     input_df = pd.DataFrame(list(molecules_to_predict.values('inchikey', 'smiles')))
# #     print("Input DataFrame created.")

# #     # Featurize and predict ADMET properties
# #     admet_df = featurize_and_predict_admet(input_df)
# #     print("ADMET properties predicted.")

# #     # Upload the predicted ADMET properties to the database
# #     for _, row in admet_df.iterrows():
# #         admet_properties = AdmetProperties(
# #             inchikey_id=row['inchikey'],
# #             ames=row['ames'],
# #             bbb=row['bbb'],
# #             bioavailability=row['bioavailability'],
# #             caco2=row['caco2'],
# #             clearance_hepatocyte=row['clearance_hepatocyte'],
# #             clearance_microsome=row['clearance_microsome'],
# #             cyp2c9=row['cyp2c9'],
# #             cyp2c9_substrate=row['cyp2c9_substrate'],
# #             cyp2d6=row['cyp2d6'],
# #             cyp2d6_substrate=row['cyp2d6_substrate'],
# #             cyp3a4=row['cyp3a4'],
# #             cyp3a4_substrate=row['cyp3a4_substrate'],
# #             dili=row['dili'],
# #             half_life=row['half_life'],
# #             herg=row['herg'],
# #             hia=row['hia'],
# #             ld50=row['ld50'],
# #             lipophilicity=row['lipophilicity'],
# #             pgp=row['pgp'],
# #             ppbr=row['ppbr'],
# #             solubility=row['solubility'],
# #             vdss=row['vdss']
# #         )
# #         admet_properties.save()
# #     print("ADMET properties uploaded to the database.")

# # def compare_table_row_counts():
# #     generated_flavonoids_count = GeneratedFlavonoid.objects.count()
# #     admet_properties_count = AdmetProperties.objects.count()

# #     if generated_flavonoids_count == admet_properties_count:
# #         return True
# #     else:
# #         return False

# #----------------------------------------------------------------------------------
# from predictors.models import AdmetProperties
# from predictors.utils.admet_tools import featurize_and_predict_admet

# from molecular_generators.models import GeneratedFlavonoid


# def get_inchikeys_to_predict():
#     # Get all InChIKeys from the generated_flavonoids table
#     generated_inchikeys = GeneratedFlavonoid.objects.values_list('inchikey', flat=True)

#     # Get all InChIKeys from the admet_properties table
#     existing_inchikeys = AdmetProperties.objects.values_list('inchikey', flat=True)

#     # Find the InChIKeys that are in generated_flavonoids but not in admet_properties
#     inchikeys_to_predict = generated_inchikeys.exclude(inchikey__in=existing_inchikeys)

#     return list(inchikeys_to_predict)

# def predict_admet_properties(input_df):
#     # Featurize and predict ADMET properties
#     admet_df = featurize_and_predict_admet(input_df)
#     print("ADMET properties predicted.")
#     return admet_df

# def upload_admet_properties(admet_df):
#     # Upload the predicted ADMET properties to the database
#     for _, row in admet_df.iterrows():
#         admet_properties = AdmetProperties(
#             inchikey_id=row['inchikey'],
#             ames=row['ames'],
#             bbb=row['bbb'],
#             bioavailability=row['bioavailability'],
#             caco2=row['caco2'],
#             clearance_hepatocyte=row['clearance_hepatocyte'],
#             clearance_microsome=row['clearance_microsome'],
#             cyp2c9=row['cyp2c9'],
#             cyp2c9_substrate=row['cyp2c9_substrate'],
#             cyp2d6=row['cyp2d6'],
#             cyp2d6_substrate=row['cyp2d6_substrate'],
#             cyp3a4=row['cyp3a4'],
#             cyp3a4_substrate=row['cyp3a4_substrate'],
#             dili=row['dili'],
#             half_life=row['half_life'],
#             herg=row['herg'],
#             hia=row['hia'],
#             ld50=row['ld50'],
#             lipophilicity=row['lipophilicity'],
#             pgp=row['pgp'],
#             ppbr=row['ppbr'],
#             solubility=row['solubility'],
#             vdss=row['vdss']
#         )
#         admet_properties.save()
#     print("ADMET properties uploaded to the database.")

from predictors.models import ProteinTargetPrediction


def upload_binding_affinities(df):
    # Upload the predicted binding affinities to the database
    for _, row in df.iterrows():
        prediction = ProteinTargetPrediction(
            inchikey_id=row['inchikey'],
            smiles=row['smiles'],
            alox5_neg_log10_affinity_M=row['alox5_neg_log10_affinity_M'],
            alox5_affinity_uM=row['alox5_affinity_uM'],
            casp1_neg_log10_affinity_M=row['casp1_neg_log10_affinity_M'],
            casp1_affinity_uM=row['casp1_affinity_uM'],
            cox1_neg_log10_affinity_M=row['cox1_neg_log10_affinity_M'],
            cox1_affinity_uM=row['cox1_affinity_uM'],
            flap_neg_log10_affinity_M=row['flap_neg_log10_affinity_M'],
            flap_affinity_uM=row['flap_affinity_uM'],
            jak1_neg_log10_affinity_M=row['jak1_neg_log10_affinity_M'],
            jak1_affinity_uM=row['jak1_affinity_uM'],
            jak2_neg_log10_affinity_M=row['jak2_neg_log10_affinity_M'],
            jak2_affinity_uM=row['jak2_affinity_uM'],
            lck_neg_log10_affinity_M=row['lck_neg_log10_affinity_M'],
            lck_affinity_uM=row['lck_affinity_uM'],
            magl_neg_log10_affinity_M=row['magl_neg_log10_affinity_M'],
            magl_affinity_uM=row['magl_affinity_uM'],
            mpges1_neg_log10_affinity_M=row['mpges1_neg_log10_affinity_M'],
            mpges1_affinity_uM=row['mpges1_affinity_uM'],
            pdl1_neg_log10_affinity_M=row['pdl1_neg_log10_affinity_M'],
            pdl1_affinity_uM=row['pdl1_affinity_uM'],
            trka_neg_log10_affinity_M=row['trka_neg_log10_affinity_M'],
            trka_affinity_uM=row['trka_affinity_uM'],
            trkb_neg_log10_affinity_M=row['trkb_neg_log10_affinity_M'],
            trkb_affinity_uM=row['trkb_affinity_uM'],
            tyk2_neg_log10_affinity_M=row['tyk2_neg_log10_affinity_M'],
            tyk2_affinity_uM=row['tyk2_affinity_uM']
        )
        prediction.save()