"""Open files with information of samples and expression values"""
import constants.constants
from utils.utils import Utils
import sys
import itertools
import pandas as pd
from constants.constants import Constants
import csv

class Tissue(object):

    def __init__(self, *args):
        self.tissue, self.age, self.sex = args
        self.dt_gene = {}
        self.count_codon_expression = {}

    def add_value(self, gene, value):

        if gene not in self.dt_gene:
            self.dt_gene[gene] = float(value)
        ##compor repetição dos genes
        else:
            sys.exit(f"Error gene {gene} already exists.")
        return self.dt_gene

    def get_number_gene(self):
        return len(self.dt_gene)


class Sample(object):

    def __init__(self):
        """
		:parm sample_name
		"""
        self.dt_sample = {}  ## { sample_name : tissue, sample_name1 : tissue, sample_name2 : tissue, ...

    def add_sample(self, sample_name, tissue, age, sex):

        if sample_name not in self.dt_sample:
            self.dt_sample[sample_name] = Tissue(tissue, age, sex)
        else:
            sys.exit("Sample already exist: " + sample_name)

        return self.dt_sample

    def get_number_gene(self, sample_name):
        return self.dt_sample[sample_name].get_number_gene()

    def add_gene(self, sample_name, gene, value):
        if sample_name not in self.dt_sample: sys.exit("Sample does not exist: " + sample_name)
        self.dt_sample[sample_name].add_value(gene, value)

    def get_number_sample(self):
        return len(self.dt_sample)


class Expression(object):
    utils = Utils()
    constants = Constants()

    def __init__(self, sample_info, sample_expression):
        """
		:param sample_info			info to the samples
		:param sample_expression   file with expression values
		"""
        self.most_dif_expressed = {}
        self.sample = Sample()

        # set the names of the files
        self.file_information = sample_info
        self.file_expression = sample_expression

        # test if files exists
        self.utils.test_exist_file(self.file_information)
        self.utils.test_exist_file(self.file_expression)

        # read files
        self.__samples_information()
        self.__expression_values()

    def __str__(self):
        return f"samples: {self.sample.get_number_sample()} genes: {self.sample.get_number_gene()}"

    def get_number_sample(self):
        return self.sample.get_number_sample()

    def get_number_gene(self, sample_name):
        return self.sample.get_number_gene(sample_name)

    def most_differentially_expressed_genes(self, sample_name1, sample_name2):
        dif_expression = sorted(
            [abs(self.sample.dt_sample[sample_name2].dt_gene[key] - self.sample.dt_sample[sample_name1].dt_gene[key])
            for key in self.sample.dt_sample[sample_name1].dt_gene.keys()], reverse=True)
        dif_expression_dict = {}
    
        for dif in dif_expression:
            for key in self.sample.dt_sample[sample_name1].dt_gene.keys():
                if dif == abs(self.sample.dt_sample[sample_name2].dt_gene[key] - \
						self.sample.dt_sample[sample_name1].dt_gene[key]):
                    dif_expression_dict[key] = dif
                self.most_dif_expressed = dict(itertools.islice(dif_expression_dict.items(), 100))
        return self.most_dif_expressed

    def counts_with_expression(self, sample, counts):
        try:
            most_expressed_counts = {gene: {codon: self.sample.dt_sample[sample].dt_gene[gene] * counts[gene][codon]
                                for codon in list(counts[gene].keys())} for gene in counts.keys() if
                                gene != 'genome' and gene in self.sample.dt_sample[sample].dt_gene }
        except KeyError as e:
            print(str(e))
            sys.exit("Error")


        dataframe_counts_expression = pd.DataFrame.from_dict(data=most_expressed_counts, orient='index')
        totals = dataframe_counts_expression.sum(axis=0).T
        dataframe_counts_expression.loc['Total'] = totals
        #print(dataframe_counts_expression)
        return dataframe_counts_expression

    def compare_T0_T1(self, dataframe1, dataframe2):
        dif = []
        for codon in dataframe1:
            dif.append(dataframe2[codon]['Total'] - dataframe1[codon]['Total'])

        dataframe_dif = pd.DataFrame(dif, columns=['Total'], index=Constants.TOTAL_CODONS)

        return dataframe_dif


    def compare_counts(self, file1, file2):
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        not_patterns = dict()
        patterns = dict()
        for codon in df1:

            if df1[codon][0] != 'Total':
                if (float(df1[codon][0]) > 0 and float(df2[codon][0] < 0)):
                    not_patterns[codon] =  ('Decrease', df1[codon][0], df2[codon][0])
                elif(float(df1[codon][0]) < 0 and float(df2[codon][0] > 0)):
                    not_patterns[codon] = ('Increase', df1[codon][0], df2[codon][0])
                else:
                    if float(df1[codon][0]) > 0:
                        patterns[codon] = ('Increase', Constants.codons_per_aminoacid[codon.upper().replace('T', 'U')])
                    else:
                        patterns[codon] = ('Decrease', Constants.codons_per_aminoacid[codon.upper().replace('T', 'U')])
        patterns_df = pd.DataFrame(patterns)
        print(patterns)
        print(patterns_df)
        return patterns_df

    def ilustrate_patterns(self, dataframe):
        direction = dict()
        for codon in dataframe:
            if dataframe[codon][0] == 'Increase':
                direction[codon].append('↗')
            else:
                direction[codon].append('↘')
        print(direction)
        return direction


    def __samples_information(self):
        """Open, read and save information from samples
		File:
			Sample	Tissue	Age	Sex
		A9_384Bulk_Plate1_S9	A9_384Bulk_Plate1_S9	Brain	3	Male
		A20_384Bulk_Plate2_S20	A20_384Bulk_Plate2_S20	Brain	12	Male
		E20_384Bulk_Plate1_S116	E20_384Bulk_Plate1_S116	Brain	21	Male

		"""
        main_header = ["Sample", "Tissue", "Age", "Sex"]
        with open(self.file_information, 'r') as information_file:
            file = information_file.readlines()

            ## test header
            if len(file) > 0 and len(file[0]) > len(main_header):
                lst_data = file[0].strip().split("\t")
                if len(lst_data) != len(main_header):
                    sys.exit("Error: header must be '{}'".format("','".join(main_header)))

                for i in range(len(lst_data)):
                    if lst_data[i] != main_header[i]:
                        sys.exit("Error: header must be '{}'".format("','".join(main_header)))

            ## read body file
            lines = file[1:]
            for line in lines:

                line = line.strip()
                if len(line) == 0: continue

                lst_line = line.split('\t')
                if len(lst_line) != len(lst_data) + 1:
                    sys.exit("Wrong line: " + line)
                self.sample.add_sample(lst_line[0], lst_line[2], lst_line[3], lst_line[4])

    def __expression_values(self):
        """Open, read and save values of expression from the different samples

		File:
		  A9_384Bulk_Plate1_S9  A20_384Bulk_Plate2_S20 E20_384Bulk_Plate1_S116 F11_384Bulk_Plate2_S131 L19_384Bulk_Plate2_S283 A18_384Bulk_Plate1_S18
		lcl|NC_000913.3_cds_NP_414542.1_1 2109.15514707196	2347.99870835941	1745.26561494942	1901.93328159018	2160.2317805036

		"""

        genes = 0
        with open(self.file_expression, 'r') as expression_file:

            dict_header = {}
            for line in expression_file:

                line = line.strip()
                if len(line) == 0: continue

                lst_data = line.split()
                if len(dict_header) == 0 and len(lst_data) > 0:
                    dict_header = dict(zip(range(len(lst_data)), lst_data))
                elif (len(lst_data) > 0):  ## body
                    ### test number of samples
                    if (len(lst_data) - 1 != self.sample.get_number_sample()):
                        print(self.sample.get_number_sample())
                        sys.exit("Error: line doesn't have values for all samples - " + line)

                    gene = None
                    for index, value in enumerate(lst_data):
                        if gene is None:
                            gene = value
                        else:
                            self.sample.add_gene(dict_header[index - 1], gene, value)

                    genes += 1
