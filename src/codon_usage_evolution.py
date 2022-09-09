'''
Created on 09/09/2022

@author: mmp
'''
import os
import gzip
import pandas as pd
from Bio import SeqIO

def read_genome_and_count_codons(file_name):
    with gzip.open(file_name, 'rt') as file:
        sequences = []  # Setup an empty list
        for record in SeqIO.parse(file, "fasta"):
            #print(record.description)
            #print( (record.seq))
            #print(len(record))
            sequences.append(record)
        genes=[]
        for record in sequences:
            element = record.description.split(" ") 
            gene = element[1][6:-1]
            genes.append(gene)
        #print (genes)

    dic_genes={}
    for gene in genes:
        for record in sequences:
            if gene in record.description:
                dic_genes[gene]=record.seq
    #print(dic_genes)

#fazer a contagem de cada codão para cada gene e criar novo dicionário 
    codoes={}

    for key,value in dic_genes.items():
        codoes[key]=[value[i:i + 3] for i in range(0, len(value), 3)]

    #print(codoes)

    total_codoes=["UUU", "UUC", "UUA", "UUG", "UCU", "UCC", "UCA", "UCG", "UAU", "UAC","UAA", 
           "UAG", "UGU", "UGC", "UGA", "UGG", "CUU", "CUC", "CUA", "CUG", "CCU", "CCC", "CCA",
          "CCG", "CAU", "CAC", "CAA", "CAG", "CGU", "CGC", "CGA", "CGG", "AUU", "AUC", "AUA",
         "AUG", "ACU", "ACC", "ACA", "ACG", "AAU", "AAC", "AAA", "AAG", "AGU", "AGC", "AGA", 
         "AGG", "GUU", "GUC", "GUA", "GUG", "GCU", "GCC", "GCA", "GCG", "GAU", "GAC", "GAA", 
         "GAG", "GGU", "GGC", "GGA", "GGG"]

##obter contagens 
    contagens={}
    counts=[]
    for i in total_codoes:
        for key,value in codoes.items():
            if i not in value: 
                x=0
                if key not in contagens:
                    contagens[key]=[x]
                else:
                    contagens[key].append(x)
    
            else:
                x=value.count(i)
                if key not in contagens:
                    contagens[key]=[x]
                else:
                    contagens[key].append(x)
    #print(contagens)

#criar dataframe com as contagens
    data = [[key]+ i for key, i in contagens.items()]
    column_labels = ["Gene\Codon"] + [i for i in total_codoes]
    dataframe_genome = pd.DataFrame(data, columns=column_labels)
    return dataframe_genome
    #print(dataframe_genome)

def save_table(dataframe_genome, file_out):
    dataframe_genome.to_csv(file_out)
   
if __name__ == '__main__':

    ## set file name in and out
    #base_path = "/home/project/master/codon_usage"
    base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"
    
    file_name_in = "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"
    file_name_out = "table_out.csv"
    
    ## get dataframe for genes
    dataframe_genome = read_genome_and_count_codons( os.path.join(base_path, file_name_in) )
    
    ## save
    save_table(dataframe_genome, os.path.join(base_path, file_name_out))
    
    print("finished")