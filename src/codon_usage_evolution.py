'''
Created on 09/09/2022

@author: mmp
'''
import os
import gzip
import pandas as pd
from Bio import SeqIO
from CAI import RSCU, relative_adaptiveness, CAI

def read_genome(file_name):
    with gzip.open(file_name, 'rt') as file:
        gene_records = [] 
        for record in SeqIO.parse(file, "fasta"):
            #print(record.description)
            #print( (record.seq))
            #print(len(record))
            gene_records.append(record)
    return gene_records

def count_codons(gene_records):
    genes=[]
    for record in gene_records:
        element = record.description.split(" ") 
        gene = element[1][6:-1]
        genes.append(gene)
    #print (genes)
    

    dic_genes={}
    for gene in genes:
        for record in gene_records:
            if gene in record.description:
                dic_genes[gene]=record.seq
    #print(dic_genes)

    #create a dictionary with the counts for each codon per gene
    codons={}

    for key,value in dic_genes.items():
        codons[key]=[value[i:i + 3] for i in range(0, len(value), 3)]

    #print(codons)

    total_codons=["UUU", "UUC", "UUA", "UUG", "UCU", "UCC", "UCA", "UCG", "UAU", "UAC","UAA", 
           "UAG", "UGU", "UGC", "UGA", "UGG", "CUU", "CUC", "CUA", "CUG", "CCU", "CCC", "CCA",
          "CCG", "CAU", "CAC", "CAA", "CAG", "CGU", "CGC", "CGA", "CGG", "AUU", "AUC", "AUA",
         "AUG", "ACU", "ACC", "ACA", "ACG", "AAU", "AAC", "AAA", "AAG", "AGU", "AGC", "AGA", 
         "AGG", "GUU", "GUC", "GUA", "GUG", "GCU", "GCC", "GCA", "GCG", "GAU", "GAC", "GAA", 
         "GAG", "GGU", "GGC", "GGA", "GGG"]

    #obtain couts 
    counts={}
    for i in total_codons:
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
    column_labels = ["Gene\Codon"] + [i for i in total_codons]
    dataframe_counts = pd.DataFrame(data, columns=column_labels)
    #print(dataframe_counts)
    return dataframe_counts
    

def calculate_RSCU(records): 
    sequences =[] #To calculate RSCU and CAI
    for record in records:
            sequences.append(record.seq)
    #calculate RSCU for each codon 
    dic_RSCU = RSCU(sequences)
    dataframe_RSCU = pd.DataFrame([dic_RSCU])
    return dataframe_RSCU


def calculate_CAI(records):
    genes=[]
    for record in records:
        element = record.description.split(" ") 
        gene = element[1][6:-1]
        genes.append(gene)
    sequences =[] #To calculate RSCU and CAI
    for record in records:
            sequences.append(record.seq)
    weights = relative_adaptiveness(sequences)
    list_CAI = [CAI(sequence, weights=weights) for sequence in sequences]
    dic_CAI={gene:cai for gene in genes for cai in list_CAI}
    dataframe_CAI = pd.DataFrame([dic_CAI])
    return dataframe_CAI
            

def save_table(dataframe_genome, file_out):
    dataframe_genome.to_csv(file_out)
   

if __name__ == '__main__':

    ## set file name in and out
    #base_path = "/home/project/master/codon_usage"
    base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"
    
    file_name_in = "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"
    file_name_out_counts = "table_counts.csv"
    file_name_out_RSCU = "table_RSCU.csv"
    file_name_out_CAI ="table_CAI.csv"

    # get dataframes 
    dataframe_genome = count_codons(read_genome(os.path.join(base_path, file_name_in)))
    dataframe_RSCU = calculate_RSCU(read_genome(os.path.join(base_path, file_name_in)))
    dataframe_CAI = calculate_CAI(read_genome( os.path.join(base_path, file_name_in)))
    
    ## save
    save_table(dataframe_genome, os.path.join(base_path, file_name_out_counts))
    save_table(dataframe_RSCU, os.path.join(base_path, file_name_out_RSCU))
    save_table(dataframe_CAI, os.path.join(base_path, file_name_out_CAI))
    print("finished")