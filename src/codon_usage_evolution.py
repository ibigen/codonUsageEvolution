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

### instantiate two objects
utils = Utils()
constants = Constants()


def read_genome(file_name):
    """ read genome """
    print("Genome reading: " + file_name)
    with (gzip.open(file_name, mode='rt') if utils.is_gzip(file_name) else open(file_name, mode='r')) as handle_read:
        record_dict_1 = SeqIO.to_dict(SeqIO.parse(handle_read, "fasta"))
        record_dict = {}
        for key in record_dict_1:
            if len(record_dict_1[key].seq) % 3 == 0:
                record_dict[key] = record_dict_1[key]
        #print(record_dict)

        ## i can interact directly with record_dict
        data = {}  ## { gene : { condon1: 2, codon2 : 4, codon3 : 5 } }
        for key in record_dict:
            counts_gene = {}
            for i in range(0, len(record_dict[key].seq), 3):
                codon = str(record_dict[key].seq)[i:i + 3].upper().replace('T', 'U')
                if codon in counts_gene:
                    counts_gene[codon] += 1
                else:
                    counts_gene[codon] = 1

            ### add gene counts
            codon_count = [0] * len(constants.TOTAL_CODONS)
            for indexes, codon in enumerate(constants.TOTAL_CODONS):
                if codon in counts_gene: codon_count[indexes] = counts_gene[codon]
            ## add gene with count codon
            data[key] = codon_count

        # creat a dataframe with counts
        rows = [[key] + i for key, i in data.items()]
        print("Create codon data frame counts")
        column_labels = ["Gene\Codon"] + constants.TOTAL_CODONS
        dataframe_counts = pd.DataFrame(rows, columns=column_labels)

        # RSCU
        print("Create RSCUs")
        # Try to eliminate bases that do not form a codon
        initial_dic_RSCU = {key: RSCU([record_dict[key].seq]) for key in record_dict.keys()}

        # Create table sorted by amino acid.

        dic_values = list(initial_dic_RSCU.values())
        sorted_by_aminoacid = {}
        for codon in constants.TOTAL_CODONS:
            for value in range(0, len(dic_values)):
                if str(codon).upper() in dic_values[value]:
                    if codon in sorted_by_aminoacid:
                        sorted_by_aminoacid[codon].append(dic_values[value][codon])
                    else:
                        sorted_by_aminoacid[codon] = [dic_values[value][codon]]

        data_RSCU = [n for key, n in sorted_by_aminoacid.items()]
        columns_RSCU = [list(value) for value in data_RSCU]
        dataframe_RSCU = pd.DataFrame(columns_RSCU,
                                      index=[key for key in sorted_by_aminoacid.keys()],
                                      columns=[key for key in initial_dic_RSCU.keys()])
        dataframe_RSCU = dataframe_RSCU.T

        # Add CAI value at the end.
        # CAI

        sequences = [record_dict[key].seq for key in record_dict]
        weights = relative_adaptiveness(sequences)
        list_CAI = [CAI(sequence, weights=weights) for sequence in sequences]
        dic_CAI = {gene: cai for gene in record_dict.keys() for cai in list_CAI}
        dataframe_CAI = pd.DataFrame([dic_CAI])

        #dataframe_RSCU.join(list_CAI)
        #cai_column = pd.Series(dic_CAI)
        print("Creating dataframe ")
        dataframe_RSCU_and_CAI = pd.concat([dataframe_RSCU, dataframe_CAI.T])
        dataframe_RSCU_and_CAI.rename(columns={0: 'CAI'}, inplace=True)


    return dataframe_counts, dataframe_RSCU_and_CAI, dataframe_CAI


def save_table(dataframe_genome, file_out):
    """ save a data frame to a CSV file """
    dataframe_genome.to_csv(file_out)


if __name__ == '__main__':

    ### several utilities
    utils = Utils()

    ## set file name in and out
    if (socket.gethostname() == "cs-nb0008"):  ## test computer name
        base_path = "/home/projects/ua/master/codon_usage"
    else:
        base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"

    file_name_in = os.path.join(base_path, "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz")
    file_name_out_counts = "table_counts.csv"
    file_name_out_RSCU_CAI = "table_RSCU_CAI.csv"
    file_name_out_CAI = "table_CAI.csv"

    # testing existing files
    utils.test_exist_file(file_name_in)

    # get dataframes
    dataframe_count_codons_in_genes, dataframe_RSCU_CAI, dataframe_CAI = read_genome(file_name_in)

    # save
    save_table(dataframe_count_codons_in_genes, os.path.join(base_path, file_name_out_counts))
    save_table(dataframe_RSCU_CAI, os.path.join(base_path, file_name_out_RSCU_CAI))
    save_table(dataframe_CAI, os.path.join(base_path, file_name_out_CAI))

    # make expression in genes 
    print("finished")
