from django.db import models

# Create your models here.

from django.db import models

# Create your models here.

class GeneratedFlavonoid(models.Model):
    inchikey = models.CharField(max_length=100, primary_key=True)
    smiles = models.CharField(max_length=255, default='')
    molecular_weight = models.FloatField(default=0.0)
    nhet = models.IntegerField(default=0)
    nrot = models.IntegerField(default=0)
    nring = models.IntegerField(default=0)
    nha = models.IntegerField(default=0)
    nhd = models.IntegerField(default=0)
    logp = models.FloatField(default=0.0)

    class Meta:
        db_table = 'generated_flavonoids'