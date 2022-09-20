'''
Created on 09/09/2022

@author: mmp
'''
import gzip
import os
import pandas as pd
import socket
from constants.constants import Constants
from utils.utils import Utils
from Bio import SeqIO
from CAI import RSCU, relative_adaptiveness, CAI

### instanciate two objects
utils = Utils()
constants = Constants()

def read_genome( file_name):
    """ read genome """
    with (gzip.open(file_name, mode='rt') if utils.is_gzip(file_name) else open(file_name, mode='r')) as handle_read:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle_read, "fasta"))
        dic_genes={key:record_dict[key].seq for key in record_dict}

    #create a dictionary with the counts for each codon per gene
    codons={}

    for key,value in dic_genes.items():
        codons[key]=[value[i:i + 3] for i in range(0, len(value), 3)]

    #obtain couts 
    counts={}
    for i in constants.TOTAL_CODONS:
        for key,value in codons.items():
            if i not in value: 
                x=0
                if key not in counts:
                    counts[key]=[x]
                else:
                    counts[key].append(x)
    
            else:
                x=value.count(i)
                if key not in counts:
                    counts[key]=[x]
                else:
                    counts[key].append(x)
    #print(counts)

    
    #creat a dataframe with counts
    data = [[key]+ i for key, i in counts.items()]
    column_labels = ["Gene\Codon"] + constants.TOTAL_CODONS
    dataframe_counts = pd.DataFrame(data, columns=column_labels)

    #RSCU
    sequences = [value for value in dic_genes.values()]
    dic_RSCU = RSCU(sequences)
    dataframe_RSCU = pd.DataFrame([dic_RSCU])

    #CAI
    weights = relative_adaptiveness(sequences)
    list_CAI = [CAI(sequence, weights=weights) for sequence in sequences]
    dic_CAI={gene:cai for gene in dic_genes.keys() for cai in list_CAI}
    dataframe_CAI = pd.DataFrame([dic_CAI])

    print(dataframe_counts)
    return dataframe_counts, dataframe_RSCU, dataframe_CAI
    

def save_table( dataframe_genome, file_out):
    """ save a data frame to a CSV file """
    dataframe_genome.to_csv(file_out)
   

if __name__ == '__main__':
    
    ### several utilities
    utils = Utils()
    
    ## set file name in and out
    if (socket.gethostname() == "cs-nb0008"): ## test computer name
        base_path = "/home/projects/master/codon_usage"
    else:
        base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"
    
    file_name_in = os.path.join(base_path, "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz")
    file_name_out_counts = "table_counts.csv"
    file_name_out_RSCU = "table_RSCU.csv"
    file_name_out_CAI ="table_CAI.csv"

    ## testing existing files
    utils.test_exist_file(file_name_in)
    
    # get dataframes
    dataframe_count_codons_in_genes, dataframe_RSCU, dataframe_CAI = read_genome(file_name_in)
    

    ## save
    save_table(dataframe_count_codons_in_genes, os.path.join(base_path, file_name_out_counts))
    save_table(dataframe_RSCU, os.path.join(base_path, file_name_out_RSCU))
    save_table(dataframe_CAI, os.path.join(base_path, file_name_out_CAI))
    
    ## make expression in genes 
    print("finished")