# mutation_engine_3/urls.py

from django.contrib import admin
from django.urls import include, path
from dashboard.views import index

urlpatterns = [
    path('', index, name='home'),  # Add this line
    path('generators/', include('molecular_generators.urls')),
    path('dashboard/', include('dashboard.urls')),
    path('admin/', admin.site.urls),
    path('__reload__/', include('django_browser_reload.urls')),
]