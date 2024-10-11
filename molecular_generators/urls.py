# urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('generate/', views.generate_view, name='generate_view'),
    path('process/', views.generate_and_process_mutants, name='generate_and_process_mutants'),
    path('generate/predict/', views.predict_admet, name='predict_admet'),
    path('debug/', views.debug_admet_prediction, name='debug_admet_prediction'),
    path('predict-binding-affinities/', views.predict_and_upload_binding_affinities, name='predict_and_upload_binding_affinities'),
]