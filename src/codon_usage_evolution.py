'''
Created on 09/09/2022

@author: mmp
'''
import os



def read_genome_and_count_codons(file_name):
    pass

def save_table(dataframe_genome, file_out):
    pass

if __name__ == '__main__':

    ## set file name in and out
    base_path = "/home/project/master/codon_usage"
    
    file_name_in = "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"
    file_name_out = "table_out.csv"
    
    ## get dataframe for genes
    dataframe_genome = read_genome_and_count_codons( os.path.join(base_path, file_name_in) )
    
    ## save
    save_table(dataframe_genome, os.path.join(base_path, file_name_out))
    
    print("finished")