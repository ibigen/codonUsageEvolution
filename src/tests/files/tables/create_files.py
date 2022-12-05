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
    for key in record_dict:
        gene_name = record_dict[key].description.\
                            split("[gene=")[1].split(" ")[0].replace("]", "")
        genes.append(gene_name)


data = []
for n in genes:
    x = random.uniform(0, 200)
    y = random.uniform(0, 200)
    z = random.uniform(0, 200)
    w = random.uniform(0, 200)
    v = random.uniform(0, 200)
    u = random.uniform(0, 200)
    data.append([x, y, z, w, v, u])

col_names = ['A9_384Bulk_Plate1_S9','A20_384Bulk_Plate2_S20', 'E20_384Bulk_Plate1_S116',  'F11_384Bulk_Plate2_S131',
             'L19_384Bulk_Plate2_S283', 'A18_384Bulk_Plate1_S18']

dataframe = pd.DataFrame(data, columns=col_names, index=genes)


dataframe.to_csv(r"C:\Users\Francisca\Desktop\TeseDeMestrado\E.coli_expression_values.csv")

print(dataframe)