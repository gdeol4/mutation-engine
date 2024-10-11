import sqlite3
from django.conf import settings
from .smiles_tools import is_valid_inchikey

def upload_df_to_sqlite(df):
    # Get the database path from Django settings
    db_path = settings.DATABASES['default']['NAME']

    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Prepare the SQL query
    query = '''INSERT OR REPLACE INTO generated_flavonoids
                (inchikey, smiles, molecular_weight, nhet, nrot, nring, nha, nhd, logp)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)'''

    # Iterate over the DataFrame rows and insert/replace data
    for _, row in df.iterrows():
        inchikey = row['inchikey']
        smiles = row['smiles']
        molecular_weight = row['molecular_weight']
        nhet = row['nhet']
        nrot = row['nrot']
        nring = row['nring']
        nha = row['nha']
        nhd = row['nhd']
        logp = row['logp']

        if is_valid_inchikey(inchikey):
            c.execute(query, (inchikey, smiles, molecular_weight, nhet, nrot, nring, nha, nhd, logp))
        else:
            print(f"Invalid InChIKey: {inchikey}")

    # Commit the changes and close the connection
    conn.commit()
    conn.close()