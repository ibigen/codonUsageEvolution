import gzip
import os
import pandas as pd
from Bio import SeqIO
from utils.utils import Utils
import random
utils = Utils()

base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"
name = "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"
file = os.path.join(base_path, name)

with (gzip.open(file, mode='rt') if utils.is_gzip(file) else open(file, mode='r')) as handle_read:
    genes = []
    record_dict = SeqIO.to_dict(SeqIO.parse(handle_read, "fasta"))
    for record in record_dict:
        genes.append(record)


data = []
for n in genes:
    x = random.uniform(0, 200)
    y = random.uniform(0, 200)
    z = random.uniform(0, 200)
    w = random.uniform(0, 200)
    data.append([x, y, z, w])

col_names = ['A9_384Bulk_Plate1_S9', 'E20_384Bulk_Plate1_S116',  'F11_384Bulk_Plate2_S131', 'A18_384Bulk_Plate1_S18']

dataframe = pd.DataFrame(data, columns=col_names, index=genes)


dataframe.to_csv(r"C:\Users\Francisca\Desktop\TeseDeMestrado\E.coli_expression_values.csv")

print(dataframe)