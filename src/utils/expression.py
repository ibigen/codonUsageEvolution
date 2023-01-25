"""Open files with information of samples and expression values"""
import numpy as np

from utils.utils import Utils
import sys
import itertools
import pandas as pd
from constants.constants import Constants
import glob, os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.cm import ScalarMappable
import seaborn as sb

class Tissue(object):

    def __init__(self, *args):
        self.tissue, self.age, self.sex = args
        self.dt_gene = {}
        self.count_codon_expression = {}

    def add_value(self, gene, value):

        if gene not in self.dt_gene:
            self.dt_gene[gene] = float(value)

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
        :param sample_info            info to the samples
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
                if dif == abs(self.sample.dt_sample[sample_name2].dt_gene[key] -
                              self.sample.dt_sample[sample_name1].dt_gene[key]):
                    dif_expression_dict[key] = dif
                self.most_dif_expressed = dict(itertools.islice(dif_expression_dict.items(), 100))
        return self.most_dif_expressed

    def counts_with_expression(self, sample, counts, **kwargs):

        #print(kwargs)
        multi = kwargs['multi']
        if multi:
            media = kwargs['media']
            try:
                most_expressed_counts = {gene: {codon: media[gene] * counts[gene][codon]
                                            for codon in list(counts[gene].keys())} for gene in counts.keys() if
                                     gene != 'genome' and gene in media.keys()}
            except KeyError as e:
                print(str(e))
                sys.exit("Error")
        else:

            try:
                most_expressed_counts = {gene: {codon: self.sample.dt_sample[sample].dt_gene[gene] * counts[gene][codon]
                                            for codon in list(counts[gene].keys())} for gene in counts.keys() if
                                     gene != 'genome' and gene in self.sample.dt_sample[sample].dt_gene}
            except KeyError as e:
                print(str(e))
                sys.exit("Error")

        dataframe_counts_expression = pd.DataFrame.from_dict(data=most_expressed_counts, orient='index')
        totals = dataframe_counts_expression.sum(axis=0).T
        dataframe_counts_expression.loc['Total'] = totals
        return dataframe_counts_expression

    def compare_timepoints(self, df1, df0):
        dif = {}
        for codon in df1:
            dif[codon] = [df1[codon]['Total'] - df0[codon]['Total']]

        dataframe_dif = pd.DataFrame.from_dict(dif, orient='index')
        return dataframe_dif

    def compare_counts(self, folder, samples):
        files_lst = []
        for n, file in enumerate(glob.glob(os.path.join(folder, "Counts-with-expression-*.csv"))):
            print(samples[n], file)
            """
A9_384Bulk_Plate1_S9 /home/mmp/git/codonUsageEvolution/src/tests/files/result/Counts-with-expression-A9_384Bulk_Plate1_S9.csv
A20_384Bulk_Plate2_S20 /home/mmp/git/codonUsageEvolution/src/tests/files/result/Counts-with-expression-A20_384Bulk_Plate2_S20.csv
E20_384Bulk_Plate1_S116 /home/mmp/git/codonUsageEvolution/src/tests/files/result/Counts-with-expression-E20_384Bulk_Plate1_S116.csv
F11_384Bulk_Plate2_S131 /home/mmp/git/codonUsageEvolution/src/tests/files/result/Counts-with-expression-A18_384Bulk_Plate1_S18.csv
L19_384Bulk_Plate2_S283 /home/mmp/git/codonUsageEvolution/src/tests/files/result/Counts-with-expression-L19_384Bulk_Plate2_S283.csv
A18_384Bulk_Plate1_S18 /home/mmp/git/codonUsageEvolution/src/tests/files/result/Counts-with-expression-F11_384Bulk_Plate2_S131.csv
            """
            file_df = pd.read_csv(file, index_col=0, sep=',')  # , index_col=0
            # file_df = file_df.split('\n')
            files_lst.append(file_df)

        patterns = {}
        for n, dataframe in enumerate(files_lst):

            for value in dataframe:
                print(value, files_lst[n - 1][value]['Total'], dataframe[value]['Total'])
                if files_lst[n - 1][value]['Total'] < dataframe[value]['Total']:
                    if value not in patterns:
                        patterns[value] = ['Increase']
                    else:
                        patterns[value] += ['Increase']
                else:
                    if value not in patterns:
                        patterns[value] = ['Decrease']
                    else:
                        patterns[value] += ['Decrease']

        columns = [f'{samples[n-1]}_{sample}' for n, sample in enumerate(samples)]

        data = [n for key, n in patterns.items()]
        final_dataframe = pd.DataFrame(data, columns=columns, index=[key for key in patterns.keys()])
        #print(final_dataframe)
        return final_dataframe

    def ilustrate_patterns(self, patterns_lst):
        direction = {}
        for sample in patterns_lst:
            for codon in patterns_lst[sample]:
                if codon == 'Increase':
                    if sample not in direction:
                        direction[sample] = ['+']
                    else:
                        direction[sample] += ['+']
                else:
                    if sample not in direction:
                        direction[sample] = ['-']
                    else:
                        direction[sample] += ['-']
        dataframe_direction = pd.DataFrame(direction, columns=[key for key in direction.keys()], index=Constants.TOTAL_CODONS)
        #print(dataframe_direction)
        return dataframe_direction


    def plot_counts(self, lst_counts, samples, b_ecoli, test):
        if b_ecoli:
            if test:
                directory = r'C:\Users\Francisca\Desktop\TeseDeMestrado\test'
            else:
                directory = r'C:\Users\Francisca\Desktop\TeseDeMestrado\ecoli'
        else:
            directory = r'C:\Users\Francisca\Desktop\TeseDeMestrado\mouse'
        dic_codons = {}
        final_dict = {}
        for n, dataframe in enumerate(lst_counts):
            dic_codons[samples[n]] = {}
            for codon in dataframe:
                if codon not in dic_codons[samples[n]]:
                    dic_codons[samples[n]][codon] = dataframe[codon]['Total']
                else:
                    dic_codons[samples[n]][codon].append(dataframe[codon]['Total'])

        data = pd.DataFrame(dic_codons, columns=[sample for sample in samples])
        codons = Constants.TOTAL_CODONS
        data['Codon'] = codons
        df = pd.melt(data, id_vars='Codon', value_vars=samples, value_name='Counts', ignore_index=True)
        max = 0
        min = 100000
        for value in df['Counts']:
            if type(value) == float:
                if value > max:
                    max = value
                elif value < min:
                    min = value


        norm = TwoSlopeNorm(vcenter=(max-min)/2, vmin=min-100, vmax=max+100)
        print(norm)
        #cmap = plt.get_cmap('PuBuGn')
        #cmap = plt.get_cmap('YlGnBu')
        #cmap = plt.get_cmap('brg')
        cmap = plt.get_cmap('brg')


        def my_bar_plot(x, y, **kwargs):
            plt.barh(y=y, width=np.abs(x), color=cmap(norm(x)))


        g = sb.FacetGrid(data=df, col='variable', height=9, aspect=0.2, sharey=True)
        g.map(my_bar_plot, 'Counts', 'Codon')
        g.fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', ax=g.axes, fraction=0.1, shrink=0.2)
        plt.show()


    def __samples_information(self):
        """Open, read and save information from samples
        File:
            Sample    Tissue    Age    Sex
        A9_384Bulk_Plate1_S9    A9_384Bulk_Plate1_S9    Brain    3    Male
        A20_384Bulk_Plate2_S20    A20_384Bulk_Plate2_S20    Brain    12    Male
        E20_384Bulk_Plate1_S116    E20_384Bulk_Plate1_S116    Brain    21    Male

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
        lcl|NC_000913.3_cds_NP_414542.1_1 2109.15514707196    2347.99870835941    1745.26561494942    1901.93328159018    2160.2317805036

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
