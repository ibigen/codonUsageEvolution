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
    print("Genome reading: " + file_name)
    with (gzip.open(file_name, mode='rt') if utils.is_gzip(file_name) else open(file_name, mode='r')) as handle_read:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle_read, "fasta"))
        ### record_dict == dic_genes
        ## dic_genes={key:record_dict[key].seq for key in record_dict}

        #create a dictionary with the counts for each codon per gene
        #codons={}
    
        ## can do this when 
        #for key,value in record_dict.items():
        #    codons[key]=[value[i:i + 3] for i in range(0, len(value), 3)]
    
        ## i can interact directly with record_dict
        data = {} ## { gene : { condon1: 2, codon2 : 4, codon3 : 5 } }
        for key in record_dict:
            counts_gene = {}
            for i in range(0, len(record_dict[key].seq), 3):
                codon = str(record_dict[key].seq)[i:i + 3].upper().replace('T', 'U')
                if codon in counts_gene: counts_gene[codon] += 1
                else: counts_gene[codon] = 1
            ### add gene counts
            codon_count =  [0] * len(constants.TOTAL_CODONS)
            for index, codon in enumerate(constants.TOTAL_CODONS):
                if codon in counts_gene: codon_count[index] = counts_gene[codon]
            
            ## add gene with count codon
            data[key] = codon_count
               
        #obtain counts 
        # counts={}
        # for i in constants.TOTAL_CODONS:
        #     for key,value in codons.items():
        #         if i not in value:
        #             x=0
        #             if key not in counts:
        #                 counts[key]=[x]
        #             else:
        #                 counts[key].append(x)
        #
        #         else:
        #             x=value.count(i)
        #             if key not in counts:
        #                 counts[key]=[x]
        #             else:
        #                 counts[key].append(x)
        #print(counts)
    
        
        #creat a dataframe with counts
#        data = [[key]+ i for key, i in counts.items()]
        column_labels = ["Gene\Codon"] + constants.TOTAL_CODONS
        dataframe_counts = pd.DataFrame(data, columns=column_labels)
    
        #RSCU
        dic_RSCU = {key:RSCU([record_dict[key].seq]) for key in record_dict.keys()}
        dataframe_RSCU = pd.DataFrame([dic_RSCU])
        #print(dic_RSCU)

        #CAI
        sequences = [seq_record.seq for seq_record in record_dict.values()]
        weights = relative_adaptiveness(sequences)
        list_CAI = [CAI(sequence, weights=weights) for sequence in sequences]
        dic_CAI={gene:cai for gene in record_dict.keys() for cai in list_CAI}
        dataframe_CAI = pd.DataFrame([dic_CAI])

    #print(dataframe_counts)
    return dataframe_counts, dataframe_RSCU, dataframe_CAI
    

def save_table( dataframe_genome, file_out):
    """ save a data frame to a CSV file """
    dataframe_genome.to_csv(file_out)
   

if __name__ == '__main__':
    
    ### several utilities
    utils = Utils()
    
    ## set file name in and out
    if (socket.gethostname() == "cs-nb0008"): ## test computer name
        base_path = "/home/projects/ua/master/codon_usage"
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