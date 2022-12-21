import gzip
import os, socket
import pandas as pd
from Bio import SeqIO
from utils.utils import Utils
import random
utils = Utils()

if socket.gethostname() == "cs-nb0008":  # test computer name
    base_path = "/home/projects/ua/master/codon_usage"
else:
    base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"
    
name = "ecoli.fasta"
file = os.path.join(base_path, name)

with (gzip.open(file, mode='rt') if utils.is_gzip(file) else open(file, mode='r')) as handle_read:
    genes = []
    record_dict = SeqIO.to_dict(SeqIO.parse(handle_read, "fasta"))
    for key in record_dict:
        gene_name = record_dict[key].description.\
                            split("[gene=")[1].split(" ")[0].replace("]", "")
        if gene_name in genes: continue
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

#col_names = ['A9_384Bulk_Plate1_S9','A20_384Bulk_Plate2_S20', 'E20_384Bulk_Plate1_S116',  'F11_384Bulk_Plate2_S131',
             #'L19_384Bulk_Plate2_S283', 'A18_384Bulk_Plate1_S18']
col_names = ['A9_384Bulk_Plate1_S9', 'A20_384Bulk_Plate2_S20', 'E20_384Bulk_Plate1_S116','F11_384Bulk_Plate2_S131', 'L19_384Bulk_Plate2_S283', 'A18_384Bulk_Plate1_S18']
dataframe = pd.DataFrame(data, columns=col_names, index=genes)


dataframe.to_csv(os.path.join(base_path, "E.coli_expression_values_test.txt"), sep ='\t')
print("Saved file: " + os.path.join(base_path, "E.coli_expression_values_test.txt"))
