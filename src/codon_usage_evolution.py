'''
Created on 09/09/2022

@author: mmp
'''
import gzip
import os
import glob
import pandas as pd
import socket
from constants.constants import Constants
from Bio import SeqIO
from CAI import RSCU, CAI
from utils.utils import Utils
from utils.count_sequences import CountSequences
from utils.expression import Expression

# instantiate two objects
utils = Utils()
constants = Constants()


def read_genome(file_name):
    """ read genome """

    counts_stats = CountSequences()
    print("Genome reading: " + file_name)
    with (gzip.open(file_name, mode='rt') if utils.is_gzip(file_name) else open(file_name,
                                                                                mode='r')) as handle_read:

        record_dict = SeqIO.to_dict(SeqIO.parse(handle_read, "fasta"))

        dt_gene_name = {}  ### {gene_name : [id, len(seq)], gene_name1 : [id1, len(seq)], ....]
        ### normalize gene IDs
        for key in record_dict:
            if record_dict[key].description.find("[gene=") > 0:
                ### get gene name
                gene_name = record_dict[key].description. \
                    split("[gene=")[1].split(" ")[0].replace("]", "")
                if gene_name in dt_gene_name:

                    if len(str(record_dict[key].seq)) > dt_gene_name[gene_name][1]:
                        dt_gene_name[gene_name] = [key, len(str(record_dict[key].seq))]
                else:
                    dt_gene_name[gene_name] = [key, len(str(record_dict[key].seq))]
            else:
                dt_gene_name[key] = [key, len(str(record_dict[key].seq))]

        # i can interact directly with record_dict
        data = {}  # { gene : { condon1: 2, codon2 : 4, codon3 : 5 } }
        initial_dic_RSCU = {}
        dic_CAI, dic_genome_CAI = {}, {}
        codon_count_total = [0] * len(constants.TOTAL_CODONS)
        print("Calculating RSCU and CAI")
        for gene_name in dt_gene_name:
            key = dt_gene_name[gene_name][0]
            if len(record_dict[key].seq) % 3 != 0:
                counts_stats.add_divisible_3()
                continue

            counts_stats.add_pass()
            counts_gene = {}
            for i in range(0, len(record_dict[key].seq) + 1, 3):
                codon = str(record_dict[key].seq)[i:i + 3].upper().replace('U', 'T')

                if codon in counts_gene:
                    counts_gene[codon] += 1
                else:
                    counts_gene[codon] = 1

            # add gene counts
            codon_count = [0] * len(constants.TOTAL_CODONS)
            for indexes, codon in enumerate(constants.TOTAL_CODONS):
                if codon in counts_gene:
                    codon_count[indexes] = counts_gene[codon]
                    codon_count_total[indexes] += counts_gene[codon]

            # add gene with count codon
            data[gene_name] = codon_count

            # RSCU
            initial_dic_RSCU[gene_name] = RSCU([record_dict[key].seq])  # {gene: {codon1: RSCU1, codon2: RSCU2}}

            # CAI
            dic_CAI[gene_name] = float(
                CAI(record_dict[key].seq, RSCUs=initial_dic_RSCU[gene_name]))  # {gene1: CAI1} {GENE2: CAI2}
        print("Calculating counts of all codons")
        # count of all codons
        data[Constants.GENOME_KEY] = codon_count_total
        # Global RSCU
        initial_dic_RSCU[Constants.GENOME_KEY] = RSCU(
            [record_dict[dt_gene_name[key][0]].seq for key in
             initial_dic_RSCU])  # {gene: {codon1: RSCU1, codon2: RSCU2}}

        # Global CAI, from RSCU from all genome
        for key in dic_CAI:
            dic_genome_CAI[dt_gene_name[key][0]] = float(
                CAI(record_dict[dt_gene_name[key][0]].seq,
                    RSCUs=initial_dic_RSCU[Constants.GENOME_KEY]))  # {gene1: CAI1} {GENE2: CAI2}
        dic_genome_CAI[Constants.GENOME_KEY] = 1

        # Global CAI, doesn't matter in this case
        dic_CAI[Constants.GENOME_KEY] = 1

        # create a dataframe with counts
        rows = [i for key, i in data.items()]
        print("Create codon counts data frame")
        column_labels = constants.TOTAL_CODONS
        dataframe_counts = pd.DataFrame(rows, columns=column_labels, index=[key for key in data.keys()])

        # Create table sorted by amino acid.
        dic_values = list(initial_dic_RSCU.values())  # [{codon1: 1, codon2:0, codon3: 1 ...}, {codon1: 1, codon2: 0,
        # codon3: 1...}, {...}, {...}]
        sorted_by_aminoacid = {}
        print("Creating dataframes to RSCU and CAI values")
        for codon in constants.TOTAL_CODONS:
            for value in range(0, len(dic_values)):
                if str(codon).upper() in dic_values[value]:
                    if codon in sorted_by_aminoacid:
                        sorted_by_aminoacid[codon].append(float(dic_values[value][codon]))
                    else:
                        sorted_by_aminoacid[codon] = [
                            float(dic_values[value][codon])]  # {codon1: [1,2,0,...], codon2: [1,0,0,2,0,...]}
        data_RSCU = [n for key, n in sorted_by_aminoacid.items()]
        columns_RSCU = [list(value) for value in data_RSCU]
        dataframe_RSCU = pd.DataFrame(columns_RSCU,
                                      index=[key for key in sorted_by_aminoacid.keys()],
                                      columns=[key for key in initial_dic_RSCU.keys()])
        dataframe_RSCU = dataframe_RSCU.T

        # CAI
        dataframe_RSCU[Constants.GENE_CAI] = dic_CAI.values()
        dataframe_RSCU[Constants.GENOME_CAI] = dic_genome_CAI.values()

    return dataframe_counts, dataframe_RSCU, counts_stats


def save_table(dataframe_genome, file_out):
    """ save a data frame to a CSV file """
    dataframe_genome.to_csv(file_out)


if __name__ == '__main__':

    # several utilities
    utils = Utils()
    b_ecoli = True
    test = True
    # set file name in and out
    if socket.gethostname() == "cs-nb0008":  # test computer name
        name = "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"  # ecoli genome
        base_path = "/home/projects/ua/master/codon_usage"
    else:
        base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"
        #name = "GCF_000001635.27_GRCm39_cds_from_genomic.fna.gz"  # mouse genome
        #name = "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"  # ecoli genome
        name = "ecoli.fasta"  # to create tables for tes

    # expression file
    if b_ecoli:
        if test:
            information_file = os.path.join(base_path, "E.coli_information_file.txt")
            expression_file = os.path.join(base_path, "E.coli_expression_values_test.txt")
        else:
            information_file = os.path.join(base_path, "E.coli_information_file.txt")
            expression_file = os.path.join(base_path, "E.coli_expression_values.txt")
    else:
        information_file = os.path.join(base_path, "coldata_brain.txt")
        expression_file = os.path.join(base_path, "norm_brain_counts.txt")

    if name == "GCF_000001635.27_GRCm39_cds_from_genomic.fna.gz":
        animal = "mouse"
    elif name == "GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz":
        animal = "ecoli"
    elif name == "ecoli.fasta":
        animal = "test"

    file_name_in = os.path.join(base_path, name)
    file_name_out_counts = f"{animal}/table_counts_{animal}.csv"
    file_name_out_RSCU_CAI = f"{animal}/table_RSCU_CAI_{animal}.csv"
    file_name_out_CAI = f"{animal}/table_CAI_{animal}.csv"

    # testing existing files
    utils.test_exist_file(file_name_in)
    utils.test_exist_file(information_file)
    utils.test_exist_file(expression_file)

    # get dataframes
    dataframe_count_codons_in_genes, dataframe_RSCU_CAI, counts_stats = read_genome(file_name_in)

    # show stats
    # print(counts_stats)

    # save
    save_table(dataframe_count_codons_in_genes, os.path.join(base_path, file_name_out_counts))
    save_table(dataframe_RSCU_CAI, os.path.join(base_path, file_name_out_RSCU_CAI))

    # make expression in genes
    print("Loading expression and samples")
    expression = Expression(information_file, expression_file)

    # analyse the different samples
    if b_ecoli:
        if test:
            temp_samples = ['A9_384Bulk_Plate1_S9', 'A20_384Bulk_Plate2_S20', 'E20_384Bulk_Plate1_S116',
                            'F11_384Bulk_Plate2_S131', 'L19_384Bulk_Plate2_S283', 'A18_384Bulk_Plate1_S18']
        else:
            temp_samples = ['A9_384Bulk_Plate1_S9', 'F11_384Bulk_Plate2_S131', 'A20_384Bulk_Plate2_S20',
                            'L19_384Bulk_Plate2_S283', 'E20_384Bulk_Plate1_S116']
    else:
        temp_samples = ['A9_384Bulk_Plate1_S9', 'B12_384Bulk_Plate3_S36',
                        'J19_384Bulk_Plate2_S235', 'A9_384Bulk_Plate2_S9', 'A20_384Bulk_Plate2_S20',
                         'E20_384Bulk_Plate1_S116', 'I20_384Bulk_Plate1_S212',
                        'P15_384Bulk_Plate2_S375', 'B12_384Bulk_Plate2_S36', 'A10_384Bulk_Plate1_S10']
    # get the list of the one hundred most differentially expressed genes between sample A9_384Bulk_Plate1_S9 and
    # E20_384Bulk_Plate1_S116
    # dt_genes_diff_expressed = expression.most_differentially_expressed_genes(sample_1, sample_2)

    print("Calculating counts with expression values")
    unor_samples = temp_samples
    same_age = {}
    media = {}
    time_points = [expression.sample.dt_sample[sample].age for sample in unor_samples]

    indexes = []
    for n, time in enumerate(time_points):
        if time_points.count(time) > 1:
            multi = True
            indexes.append(n)

    def time_point(sample):
        return int(expression.sample.dt_sample[sample].age)
    samples = sorted(unor_samples, key=time_point, reverse=False)
    print(samples)
    print(time_points)
    for n, sample in enumerate(samples):
        for i in indexes:
            if n == i:
                if len(indexes) - 1 != indexes.index(i):
                    samples.remove(sample)
                for key, value in expression.sample.dt_sample[sample].dt_gene.items():

                    if key not in same_age:
                        same_age[key] = [value, ]
                    else:
                        same_age[key].append(value)
    print(samples)
    time_points = [expression.sample.dt_sample[sample].age for sample in samples]
    print(time_points)

    for key, value in same_age.items():
        exp = [float(n) for n in value]
        if key not in media:
            media[key] = sum(exp) / len(exp)
    counts = []
    for n, sample in enumerate(samples):
        multi = False
        if temp_samples.index(sample) in indexes:
            multi = True
            counts.append(
                expression.counts_with_expression(sample, dataframe_count_codons_in_genes.to_dict(orient='index'),
                                                  multi=multi, media=media))
        else:
            counts.append(
                expression.counts_with_expression(sample, dataframe_count_codons_in_genes.to_dict(orient='index'), multi=multi))
        save_table(counts[n], os.path.join(base_path, f'{animal}/Counts-with-expression-{sample}.csv'))

    print("Comparing different time points")
    for n, dataframe in enumerate(counts):
        dif = expression.compare_timepoints(dataframe, counts[n - 1])
        # print(dif)
        save_table(dif.T, os.path.join(base_path, f'{animal}/Differences_{samples[n - 1]}_{samples[n]}.csv'))

    print('Searching for patterns')
    folder = os.path.join(base_path, f'{animal}')
    patterns = expression.compare_counts(folder, samples)
    save_table(patterns, os.path.join(base_path, f'{animal}/Patterns_between_{samples}.csv'))

    print("Illustrating patterns")
    table_direction = expression.ilustrate_patterns(patterns)
    save_table(table_direction, os.path.join(base_path, f'{animal}/Table_directions_from_{samples}.csv'))
    hist = expression.plot_counts(counts, samples, b_ecoli, test)

    print("Finished")
