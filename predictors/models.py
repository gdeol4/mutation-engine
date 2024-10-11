from django.db import models

from molecular_generators.models import GeneratedFlavonoid

# Create your models here.

class AdmetProperties(models.Model):
    inchikey = models.OneToOneField(GeneratedFlavonoid, on_delete=models.CASCADE, primary_key=True)
    ames = models.FloatField(default=0.0)
    bbb = models.FloatField(default=0.0)
    bioavailability = models.FloatField(default=0.0)
    caco2 = models.FloatField(default=0.0)
    clearance_hepatocyte = models.FloatField(default=0.0)
    clearance_microsome = models.FloatField(default=0.0)
    cyp2c9 = models.FloatField(default=0.0)
    cyp2c9_substrate = models.FloatField(default=0.0)
    cyp2d6 = models.FloatField(default=0.0)
    cyp2d6_substrate = models.FloatField(default=0.0)
    cyp3a4 = models.FloatField(default=0.0)
    cyp3a4_substrate = models.FloatField(default=0.0)
    dili = models.FloatField(default=0.0)
    half_life = models.FloatField(default=0.0)
    herg = models.FloatField(default=0.0)
    hia = models.FloatField(default=0.0)
    ld50 = models.FloatField(default=0.0)
    lipophilicity = models.FloatField(default=0.0)
    pgp = models.FloatField(default=0.0)
    ppbr = models.FloatField(default=0.0)
    solubility = models.FloatField(default=0.0)
    vdss = models.FloatField(default=0.0)

    class Meta:
        db_table = 'admet_properties'

class ProteinTargetPrediction(models.Model):
    inchikey = models.ForeignKey(GeneratedFlavonoid, on_delete=models.CASCADE)
    smiles = models.CharField(max_length=255)
    alox5_neg_log10_affinity_M = models.FloatField()
    alox5_affinity_uM = models.FloatField()
    casp1_neg_log10_affinity_M = models.FloatField()
    casp1_affinity_uM = models.FloatField()
    cox1_neg_log10_affinity_M = models.FloatField()
    cox1_affinity_uM = models.FloatField()
    flap_neg_log10_affinity_M = models.FloatField()
    flap_affinity_uM = models.FloatField()
    jak1_neg_log10_affinity_M = models.FloatField()
    jak1_affinity_uM = models.FloatField()
    jak2_neg_log10_affinity_M = models.FloatField()
    jak2_affinity_uM = models.FloatField()
    lck_neg_log10_affinity_M = models.FloatField()
    lck_affinity_uM = models.FloatField()
    magl_neg_log10_affinity_M = models.FloatField()
    magl_affinity_uM = models.FloatField()
    mpges1_neg_log10_affinity_M = models.FloatField()
    mpges1_affinity_uM = models.FloatField()
    pdl1_neg_log10_affinity_M = models.FloatField()
    pdl1_affinity_uM = models.FloatField()
    trka_neg_log10_affinity_M = models.FloatField()
    trka_affinity_uM = models.FloatField()
    trkb_neg_log10_affinity_M = models.FloatField()
    trkb_affinity_uM = models.FloatField()
    tyk2_neg_log10_affinity_M = models.FloatField()
    tyk2_affinity_uM = models.FloatField()

    class Meta:
        db_table = 'protein_target_predictions'