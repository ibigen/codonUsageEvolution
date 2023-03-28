"""Open files with information of samples and expression values"""
import numpy as np

from utils.utils import Utils
import sys, os
import itertools
import pandas as pd
from collections import OrderedDict
from constants.constants import Constants
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.cm import ScalarMappable
import seaborn as sb
from sklearn.decomposition import PCA



class Tissue(object):
    GENDER_BOTH = "BOTH"
    GENDER_MALE = "MALE"
    GENDER_FEMALE = "FEMALE"

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

    def get_list_samples(self, gender=Tissue.GENDER_BOTH):
        """ return all sample names ordered by time point """
        sample_list = list(self.dt_sample)
        if gender == Tissue.GENDER_MALE:
            sample_list = [sample for sample in list(self.dt_sample)
                           if self.dt_sample[sample].sex.lower() == Tissue.GENDER_MALE.lower()]
        if gender == Tissue.GENDER_FEMALE:
            sample_list = [sample for sample in list(self.dt_sample)
                           if self.dt_sample[sample].sex.lower() == Tissue.GENDER_FEMALE.lower()]
        return sorted(sample_list,
                      key=lambda sample: int(self.dt_sample[sample].age), reverse=False)

    def get_list_time_points(self, gender=Tissue.GENDER_BOTH):
        """ return all time points ordered and not repeated """
        timepoints = [int(self.dt_sample[sample].age) for sample in self.dt_sample]
        if gender == Tissue.GENDER_MALE:
            timepoints = [int(self.dt_sample[sample].age) for sample in self.dt_sample
                          if self.dt_sample[sample].sex.lower() == Tissue.GENDER_MALE.lower()]
        if gender == Tissue.GENDER_FEMALE:
            timepoints = [int(self.dt_sample[sample].age) for sample in self.dt_sample
                          if self.dt_sample[sample].sex.lower() == Tissue.GENDER_FEMALE.lower()]

        return sorted(list(dict.fromkeys(timepoints)), reverse=False)


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

    def save_table(self, dataframe_genome, file_out):
        """ save a data frame to a CSV file """
        dataframe_genome.to_csv(file_out)

    def get_number_sample(self):
        return self.sample.get_number_sample()

    def get_list_samples(self, gender=Tissue.GENDER_BOTH):
        """ ordered by time point """
        return self.sample.get_list_samples(gender)

    def get_list_time_points(self, gender=Tissue.GENDER_BOTH):
        """ ordered by time point and not repeated"""
        return self.sample.get_list_time_points(gender)

    def get_number_gene(self, sample_name):
        return self.sample.get_number_gene(sample_name)


    #def most_differentially_expressed_genes(self, sample_name1, sample_name2):
        #dif_expression = sorted(
            #[abs(self.sample.dt_sample[sample_name2].dt_gene[key] - self.sample.dt_sample[sample_name1].dt_gene[key])
             #for key in self.sample.dt_sample[sample_name1].dt_gene.keys()], reverse=True)
        #dif_expression_dict = {}

        #for dif in dif_expression:
            #for key in self.sample.dt_sample[sample_name1].dt_gene.keys():
                #if dif == abs(self.sample.dt_sample[sample_name2].dt_gene[key] -
                              #self.sample.dt_sample[sample_name1].dt_gene[key]):
                    #dif_expression_dict[key] = dif
                #self.most_dif_expressed = dict(itertools.islice(dif_expression_dict.items(), 100))
        #return self.most_dif_expressed

    def counts_with_expression(self, counts, average):

        try:
            most_expressed_counts = {gene: {codon: average[gene] * counts[gene][codon]
                                            for codon in list(counts[gene].keys())} for gene in counts.keys() if
                                     gene != 'genome' and gene in average.keys()}
        except KeyError as e:
            print(str(e))
            sys.exit("Error")

        dataframe_counts_expression = pd.DataFrame.from_dict(data=most_expressed_counts, orient='index')
        totals = dataframe_counts_expression.sum(axis=0).T
        dataframe_counts_expression.loc['Total'] = totals
        rscu = self.calculate_RSCU(dataframe_counts_expression)
        rscu_dataframe = pd.DataFrame(rscu, index=['RSCU'])
        final_dataframe = pd.concat([dataframe_counts_expression, rscu_dataframe], axis=0)

        return final_dataframe

    def get_counts(self, gender, codons_in_genes, b_make_averages_for_same_time_points):
        """ 
        :param gender
        :codons_in_genes {}return counts by gender making average with the same time points 
        :out   dict_samples_out { sample_name : [sample_name, sample_name same time point, sample_name same time point 2],
                                sample_name1 : [sample_name1, sample_name same time point, sample_name same time point 2], } """

        list_samples = self.get_list_samples(gender)
        dict_samples_out = OrderedDict()  ## Key, first sample of each time point, []
        return_counts = []
        ## this list
        if b_make_averages_for_same_time_points:  ## make averages...
            list_time_poins = self.get_list_time_points(gender)
            for timepoint in list_time_poins:
                fist_sample_name = None
                same_age = {}
                ## return all the sample for this time point
                for sample in [_ for _ in list_samples if int(self.sample.dt_sample[_].age) == timepoint]:

                    ### control of the samples with same time points
                    if fist_sample_name is None:
                        fist_sample_name = sample
                        dict_samples_out[fist_sample_name] = [fist_sample_name]
                    else:
                        dict_samples_out[fist_sample_name].append(sample)

                    ## create an array for average
                    for key, value in self.sample.dt_sample[sample].dt_gene.items():
                        if key not in same_age:
                            same_age[key] = [value]
                        else:
                            same_age[key].append(value)

                average = {}
                for key, list_values in same_age.items():
                    #    if key not in media:
                    average[key] = sum(list_values) / len(list_values)

                #### append
                return_counts.append(self.counts_with_expression(codons_in_genes, average))
        else:  ## return all data for a specific gender
            for sample in list_samples:
                if sample in dict_samples_out: sys.exit("Error: sample already in dictonary - " + str(sample))
                dict_samples_out[sample] = [sample]
                return_counts.append(
                    self.counts_with_expression(codons_in_genes, self.sample.dt_sample[sample].dt_gene))

        ### dictionary with expression X codons
        return return_counts, dict_samples_out

    def compare_timepoints(self, counts, samples, working_path):
        data = 'RSCU'
        differences_abs, differences = OrderedDict(), OrderedDict()   # need to be ordered
        for n, dataframe in enumerate(counts):
            key_to_process = f'{self.sample.dt_sample[samples[n - 1]].age}_{self.sample.dt_sample[samples[n]].age}'
            repeat = 1
            ### get unique name
            while key_to_process in differences_abs:
                key_to_process = f'{self.sample.dt_sample[samples[n - 1]].age}_{self.sample.dt_sample[samples[n]].age}_{repeat}'
                repeat += 1
            differences_abs[key_to_process] = {}
            differences[key_to_process] = {}
            for value in dataframe:

                if value not in differences[key_to_process]:
                    differences[key_to_process][
                        value] = counts[n - 1][value][data] - dataframe[value][data]
                    differences_abs[key_to_process][
                        value] = abs(counts[n - 1][value][data] - dataframe[value][data])

                else:
                    differences[key_to_process][
                        value] += counts[n - 1][value][data] - dataframe[value][data]
                    differences_abs[key_to_process][
                        value] += abs(counts[n - 1][value][data] - dataframe[value][data])

        ### save differences
        dataframe = pd.DataFrame(differences)
        print(
            "File with differences: " + str(os.path.join(working_path, "Differences_between_time_points.csv")))
        dataframe.to_csv(os.path.join(working_path, "Differences_between_time_points_27vs3.csv"))

        ### start making chart
        dataframe_abs = pd.DataFrame(differences_abs)

        dataframe_abs['Codon'] = [f'{str(key).upper().replace("U", "T")}_{value}' for key, value in
                                  Constants.codons_per_aminoacid.items()]
        columns = [f'{self.sample.dt_sample[samples[n - 1]].age}_{self.sample.dt_sample[sample].age}' for
                   n, sample in enumerate(samples)]
        df = pd.melt(dataframe_abs, id_vars='Codon', value_vars=columns, value_name='Difference')
        df.rename(columns={"variable": "ID"}, inplace=True)
        max_ = 0
        min_ = 100000

        for value in df['Difference']:
            if type(value) == float:
                if value > max_:
                    max_ = value
                elif value < min_:
                    min_ = value

        norm = TwoSlopeNorm(vcenter=max_ - (max_ / 2), vmin=min_, vmax=max_)
        cmap = plt.get_cmap('brg')

        def my_bar_plot(x, y, **kwargs):
            plt.barh(y=y, width=np.abs(x), color=cmap(norm(x)))

        g = sb.FacetGrid(data=df, col='ID', height=9, aspect=0.2,
                         col_order=list(dataframe.columns), sharey=True)
        g.map(my_bar_plot, 'Difference', 'Codon')
        g.fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', ax=g.axes, fraction=0.1,
                       shrink=0.2)

        # plt.title(f'Difference between Time points:{[self.sample.dt_sample[sample].age for sample in samples]} ')
        print("Create image: {}".format(os.path.join(working_path, f'Barplot_to_differences_{data}.png')))
        plt.savefig(os.path.join(working_path, f'Barplot_to_differences_{data}.png'))

        return df


    def compare_counts(self, counts, samples):
        data = 'RSCU'
        patterns = {}
        for n, dataframe in enumerate(counts):
            for value in dataframe:
                if counts[n - 1][value][data] < dataframe[value][data]:
                    if value not in patterns:
                        patterns[value] = ['Increase']
                    else:
                        patterns[value] += ['Increase']
                else:
                    if value not in patterns:
                        patterns[value] = ['Decrease']
                    else:
                        patterns[value] += ['Decrease']

        columns = []
        for n, sample in enumerate(samples):
            key_to_process = f'{self.sample.dt_sample[samples[n - 1]].age}_{self.sample.dt_sample[samples[n]].age}'
            repeat = 1
            ### get unique name
            while key_to_process in columns:
                key_to_process = f'{self.sample.dt_sample[samples[n - 1]].age}_{self.sample.dt_sample[samples[n]].age}_{repeat}'
                repeat += 1
            columns.append(key_to_process)
        data_values = [n for key, n in patterns.items()]
        final_dataframe = pd.DataFrame(data_values, columns=columns, index=[key for key in patterns.keys()])
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

        dataframe_direction = pd.DataFrame(direction, columns=[time for time in direction.keys()],
                                           index=Constants.TOTAL_CODONS)

        return dataframe_direction

    def plot_counts(self, lst_counts, samples, working_path, b_make_averages_for_same_time_points):
        data = 'RSCU'
        time_points = []
        if b_make_averages_for_same_time_points:
            time_points = [f'{self.sample.dt_sample[sample].age}' for sample in samples]
        else:
            time_points = samples

        dic_codons = OrderedDict()
        for n, dataframe in enumerate(lst_counts):
            dic_codons[time_points[n]] = {}
            for codon in dataframe:
                if codon not in dic_codons[time_points[n]]:
                    dic_codons[time_points[n]][codon] = dataframe[codon][data]
                else:
                    dic_codons[time_points[n]][codon].append(dataframe[codon][data])

        data_values = pd.DataFrame(dic_codons, columns=time_points)

        codons = Constants.TOTAL_CODONS
        data_values['Codon'] = codons

        df = pd.melt(data_values, id_vars='Codon', value_vars=time_points, value_name='Counts')

        max = 0
        min = 100000

        for value in df['Counts']:
            if type(value) == float:
                if value > max:
                    max = value
                elif value < min:
                    min = value

        norm = TwoSlopeNorm(vcenter=max - (max / 2), vmin=min, vmax=max)
        cmap = plt.get_cmap('brg')

        def my_bar_plot(x, y, **kwargs):
            plt.barh(y=y, width=np.abs(x), color=cmap(norm(x)))

        g = sb.FacetGrid(data=df, col='variable', height=9, aspect=0.2,
                         col_order=time_points, sharey=True)
        g.map(my_bar_plot, 'Counts', 'Codon')
        g.fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', ax=g.axes, fraction=0.1,
                       shrink=0.2)

        # plt.title(f'Counts to samples from {[samples for samples in samples]}')
        plt.savefig(os.path.join(working_path, f'Barplot_to_counts_{data}.png'))
        # plt.show()

    def calculate_RSCU(self, counts):
        rscu = {}
        res = {}
        for key, val in Constants.codons_per_aminoacid.items():
            if val not in res:
                res[val] = [key]
            else:
                res[val].append(key)

        for amino, codons in res.items():
            codons_T = [str(key).upper().replace('U', 'T') for key in codons]
            total = [counts[codon]['Total'] for codon in codons_T]
            for n, codon in enumerate(codons_T):
                rscu[codon] = len(codons) * ((counts[codon]['Total']) / sum(total))
                # nº de codões*(codão/soma(codões por aminoacido))
        return rscu

    def PCA_analysis(self, counts, samples, working_path):
        data = 'RSCU'
        RSCU_dic = OrderedDict()
        time_points = [f'{self.sample.dt_sample[sample].age}' for sample in samples]
        for n, dataframe in enumerate(counts):
            for codon in dataframe:
                if samples[n] not in RSCU_dic:
                    RSCU_dic[samples[n]] = [dataframe[codon][data]]
                else:
                    RSCU_dic[samples[n]].append(dataframe[codon][data])

        RSCU_dataframe = pd.DataFrame.from_dict(RSCU_dic, orient='columns')

        RSCU_dataframe['Codon'] = [str(key).upper().replace('U', 'T') for key in Constants.TOTAL_CODONS]
        RSCU_dataframe.set_index('Codon', inplace=True)
        sample_time = {}
        times = [self.sample.dt_sample[sample].age for sample in RSCU_dataframe.columns]
        RSCU_dataframe = RSCU_dataframe.transpose()
        RSCU_dataframe["time"] = times
        RSCU_dataframe = RSCU_dataframe.transpose()
        print(RSCU_dataframe)

        # Obter a linha com os time points
        time_points = RSCU_dataframe.iloc[-1, :].values

        # Remover a linha com os time points do dataframe
        RSCU_dataframe.drop(RSCU_dataframe.tail(1).index, inplace=True)

        # Executar a análise PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(RSCU_dataframe.transpose())

        # Plotar o gráfico com os pontos coloridos de acordo com os time points
        fig, ax = plt.subplots()
        for tp in np.unique(time_points):
            samples = np.where(time_points == tp)
            c = plt.cm.Set1(int(tp) / np.max([int(time_point) for time_point in time_points]))
            ax.scatter(pca_result[:, 0][samples], pca_result[:, 1][samples], color=c, label=f'Time {tp}')
        ax.legend()
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
