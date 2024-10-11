from django.http import HttpResponse
from django.shortcuts import render
from django.db.models import Q
import pandas as pd

from predictors.models import AdmetProperties, ProteinTargetPrediction
from predictors.utils.admet_tools import get_inchikeys_to_predict, predict_and_upload_admet_properties
from predictors.utils.plapt import predict_binding_affinities
from predictors.utils.queries import upload_binding_affinities
from .models import GeneratedFlavonoid
from .utils.mutation_tools import generate_mutations
from .utils.smiles_tools import process_generated_smiles
from .utils.property_calc import MolPropertyCalc
from .utils.queries import upload_df_to_sqlite


from django.db.models import Count

def debug_admet_prediction(request):
    output = ""

    # Check the InChIKeys to predict
    inchikeys_to_predict = get_inchikeys_to_predict()
    output += f"<p>Number of InChIKeys to predict: {len(inchikeys_to_predict)}</p>"
    output += f"<p>InChIKeys to predict: {inchikeys_to_predict}</p>"

    # Get the number of rows for each table
    generated_flavonoids_count = GeneratedFlavonoid.objects.count()
    admet_properties_count = AdmetProperties.objects.count()
    output += f"<p>Number of rows in GeneratedFlavonoid table: {generated_flavonoids_count}</p>"
    output += f"<p>Number of rows in AdmetProperties table: {admet_properties_count}</p>"

    # Print the last 10 rows from the GeneratedFlavonoid table
    output += "<p>Last 10 rows from GeneratedFlavonoid table:</p>"
    rows = GeneratedFlavonoid.objects.all().order_by('-pk')[:10]
    for row in rows:
        output += f"<p>{row}</p>"

    # Print the last 10 rows from the AdmetProperties table
    output += "<p>Last 10 rows from AdmetProperties table:</p>"
    rows = AdmetProperties.objects.all().order_by('-pk')[:10]
    for row in rows:
        output += f"<p>{row}</p>"



    return HttpResponse(output)

def generate_view(request):
    return render(request, 'generate.html')


def generate_and_process_mutants(request):
    if request.method == 'GET':
        # Generate mutants
        iterations = 5
        input_molecule = 'luteolin'
        sample_size = 10
        unique_mutants_df = generate_mutations(iterations, input_molecule, sample_size)

        # Process the DataFrame
        processed_df = process_generated_smiles(unique_mutants_df)

        # Calculate molecular properties
        mol_property_calc = MolPropertyCalc(processed_df)
        processed_df = mol_property_calc.calculate_properties()

        # Upload the data to the database
        upload_df_to_sqlite(processed_df)

        # Fetch the first 10 rows from the database
        rows = GeneratedFlavonoid.objects.all()[:10]

        # Get the number of rows for each table
        generated_flavonoids_count = GeneratedFlavonoid.objects.count()
        admet_properties_count = AdmetProperties.objects.count()

        print(f"Number of rows in GeneratedFlavonoid table: {generated_flavonoids_count}")
        print(f"Number of rows in AdmetProperties table: {admet_properties_count}")

        return render(request, 'result.html', {'rows': rows})
    else:
        return HttpResponse("Invalid request method.")
    
def predict_admet(request):
    # Run the prediction and upload process
    predict_and_upload_admet_properties(request)

    # Fetch the first 10 rows from the admet_properties table
    rows = AdmetProperties.objects.all()[:10]

    # Get the number of rows for each table
    generated_flavonoids_count = GeneratedFlavonoid.objects.count()
    admet_properties_count = AdmetProperties.objects.count()

    print(f"Number of rows in GeneratedFlavonoid table: {generated_flavonoids_count}")
    print(f"Number of rows in AdmetProperties table: {admet_properties_count}")

    return render(request, 'admet_results.html', {'rows': rows})

def predict_and_upload_binding_affinities(request):
    # Get the InChIKeys to predict 
    inchikeys_to_predict = get_inchikeys_to_predict(GeneratedFlavonoid, ProteinTargetPrediction, 'inchikey')
    print(f"Number of InChIKeys to predict: {len(inchikeys_to_predict)}")
    
    # Predict binding affinities
    predictions_df = predict_binding_affinities(inchikeys_to_predict)
    
    # Upload the predicted binding affinities to the database
    upload_binding_affinities(predictions_df)

    # Fetch the first 10 rows from the binding_affinities table
    rows = ProteinTargetPrediction.objects.all()[:10]

    # Get the number of rows for each table
    generated_flavonoids_count = GeneratedFlavonoid.objects.count()
    admet_properties_count = AdmetProperties.objects.count()
    binding_affinities_count = ProteinTargetPrediction.objects.count()

    print(f"Number of rows in GeneratedFlavonoid table: {generated_flavonoids_count}")
    print(f"Number of rows in AdmetProperties table: {admet_properties_count}") 
    print(f"Number of rows in BindingAffinities table: {binding_affinities_count}")

    return render(request, 'plapt_results.html', {'rows': rows})
