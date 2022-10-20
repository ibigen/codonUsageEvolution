
"""Open files with information of samples and expression values"""
from utils.utils import Utils
import sys


class Tissue(object):

    def __init__(self, tissue, age, sex):

        self.tissue = tissue
        self.age = age
        self.sex = sex
        self.dt_gene = {}

    def add_value(self, gene, value):

        if not gene in self.dt_gene:
            self.dt_gene[gene] = value
        else:
            sys.exit("Gene already in dictonary: " + gene)


class Sample(object):

    def __init__(self):
        """
        :parm sample_name
        """
        self.dt_sample = {}  ## { sample_name : tissue, sample_name1 : tissue, sample_name2 : tissue, ...

    def add_sample(self, sample_name, tissue, age, sex):

        if not sample_name in self.dt_sample:
            self.dt_sample[sample_name] = Tissue(tissue, age, sex)
        else:
            sys.exit("Sample already exist: " + sample_name)

    def add_gene(self, sample_name, gene, value):
        if not sample_name in self.dt_sample: sys.exit("Sample does not exist: " + sample_name)
        self.dt_sample[sample_name].add_value(gene, value)

    def get_number_sample(self):
        return len(self.dt_sample)


class Expression(object):
    utils = Utils()

    def __init__(self, sample_info, sample_expression):
        """
        :param sample_info			info to the samples
        :param sample_expression   file with expression values
        """
        self.sample = Sample()

        # set the names of the files
        self.file_information = sample_info
        self.file_expression = sample_expression

        # test if files exists
        self.utils.test_exist_file(self.file_information)
        self.utils.test_exist_file(self.file_expression)

        # read files
        self.samples_information()
        self.expression_values()

    def samples_information(self):
        """Open, read and save information from samples
        File:
            Sample	Tissue	Age	Sex
        A9_384Bulk_Plate1_S9	A9_384Bulk_Plate1_S9	Brain	3	Male
        A20_384Bulk_Plate2_S20	A20_384Bulk_Plate2_S20	Brain	12	Male
        E20_384Bulk_Plate1_S116	E20_384Bulk_Plate1_S116	Brain	21	Male

        """

        with open(self.file_information, 'r') as information_file:
            file = information_file.readlines()
            lines = file[1:]
            for line in lines:

                line = line.strip()
                if len(line) == 0: continue

                lst_line = line.split('\t')
                if len(lst_line) != 5:
                    sys.exit("Wrong line: " + line)
                self.sample.add_sample(lst_line[0], lst_line[2], lst_line[3], lst_line[4])

        print("Number of samples: {}".format(self.sample.get_number_sample()))
        return self.sample.get_number_sample()

    def expression_values(self):
        """Open, read and save values of expression from the different samples

        File:
          A9_384Bulk_Plate1_S9  A20_384Bulk_Plate2_S20 E20_384Bulk_Plate1_S116 F11_384Bulk_Plate2_S131 L19_384Bulk_Plate2_S283 A18_384Bulk_Plate1_S18
        lcl|NC_000913.3_cds_NP_414542.1_1 2109.15514707196	2347.99870835941	1745.26561494942	1901.93328159018	2160.2317805036	2359.80082045078

        """

        with open(self.file_expression, 'r') as expression_file:
            file = expression_file.readlines()
            samples = file[0].strip()
            for line in file[1:]:
                line = line.strip()
                if len(line) == 0: continue
                lst_line = line.split('\t')
                if len(lst_line) < 6:
                    sys.exit("Wrong line: " + line)
            for n in range(len(samples)):
                self.sample.add_sample(samples[n], lst_line[0], lst_line[n + 1])

        return len(samples)
