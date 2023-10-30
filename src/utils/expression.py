"""Open files with information of samples and expression values"""
import numpy as np
from utils.utils import Utils
import sys, os
import pandas as pd
from collections import OrderedDict
from constants.constants import Constants
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable
import seaborn as sb
from scipy.stats import chi2_contingency
#from sklearn.decomposition import PCA
#from sklearn.decomposition import FastICA

utils = Utils()

class Tissue(object):
    GENDER_BOTH = "BOTH"
    GENDER_MALE = "MALE"
    GENDER_FEMALE = "FEMALE"

    def __init__(self, *args):
        self.tissue, self.age, self.sex = args
        self.dt_gene = {}

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

    constants = Constants()

    def __init__(self, file_information, file_expression):
        """
        :param sample_info            info to the samples
        :param sample_expression   file with expression values
        """
        #self.most_dif_expressed = {}
        self.sample = Sample()

        # set the names of the files
        self.file_information = file_information
        self.file_expression = file_expression

        # test if files exists
        utils.test_exist_file(self.file_information)
        utils.test_exist_file(self.file_expression)

        # read files
        self.__samples_information()
        self.__expression_values()

    def __str__(self):
        return f"samples: {self.sample.get_number_sample()} genes: {self.sample.get_number_gene()}"

    def save_table(self, dataframe_genome, file_out):
        """ save a data frame to a CSV file """
        dataframe_genome.to_csv(file_out)

    def get_number_sample(self):
        """return number of samples"""
        return self.sample.get_number_sample()

    def get_list_samples(self, gender=Tissue.GENDER_BOTH):
        """ ordered by time point """
        return self.sample.get_list_samples(gender)

    def get_list_time_points(self, gender=Tissue.GENDER_BOTH):
        """ ordered by time point and not repeated"""
        return self.sample.get_list_time_points(gender)

    def get_number_gene(self, sample_name):
        """return number of genes"""
        return self.sample.get_number_gene(sample_name)

    # def most_differentially_expressed_genes(self, sample_name1, sample_name2):
    #   dif_expression = sorted(
    #   [abs(self.sample.dt_sample[sample_name2].dt_gene[key] - self.sample.dt_sample[sample_name1].dt_gene[key])
    #   for key in self.sample.dt_sample[sample_name1].dt_gene.keys()], reverse=True)
    #       dif_expression_dict = {}

    #       for dif in dif_expression:
    #           for key in self.sample.dt_sample[sample_name1].dt_gene.keys():
    #               if dif == abs(self.sample.dt_sample[sample_name2].dt_gene[key] -
    #                   self.sample.dt_sample[sample_name1].dt_gene[key]):
    #                   dif_expression_dict[key] = dif
    #                   self.most_dif_expressed = dict(itertools.islice(dif_expression_dict.items(), 100))
    #   return self.most_dif_expressed

    def counts_with_expression(self, counts, average):
        """
                :param average: dataframe with the b_make_averages_for_same_time_points of counts
                :codons_in_genes {}return counts by gender making b_make_averages_for_same_time_points with the same time points
                :out   dataframe with counts calculated with expression appended to RSCU values to each codon """

        try:
            counts_with_expression_dict = {gene: {codon: average[gene] * counts[gene][codon]
                                            for codon in list(counts[gene].keys())} for gene in counts.keys() if
                                     gene != 'genome' and gene in average.keys()}
        except KeyError as e:
            print(str(e))
            sys.exit("Error")

        dataframe_counts_expression = pd.DataFrame.from_dict(data=counts_with_expression_dict, orient='index')
        totals = dataframe_counts_expression.sum(axis=0).T
        dataframe_counts_expression.loc['Total'] = totals
        rscu = self.calculate_RSCU(dataframe_counts_expression)
        rscu_dataframe = pd.DataFrame(rscu, index=['RSCU'])
        final_dataframe = pd.concat([dataframe_counts_expression, rscu_dataframe], axis=0)

        return final_dataframe

    def get_counts(self, gender, codons_in_genes, b_make_averages_for_same_time_points):
        """ 
        :param gender
        :codons_in_genes {}return counts by gender making b_make_averages_for_same_time_points with the same time points
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

                    ## create an array for b_make_averages_for_same_time_points
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
        differences_abs, differences = OrderedDict(), OrderedDict()  # need to be ordered

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
                        value] = (counts[n - 1][value][data] - dataframe[value][data])

                else:
                    differences[key_to_process][
                        value] += counts[n - 1][value][data] - dataframe[value][data]
                    differences_abs[key_to_process][
                        value] += (counts[n - 1][value][data] - dataframe[value][data])

        ### save differences
        dataframe = pd.DataFrame(differences)
        print(
            "File with differences: " + str(os.path.join(working_path, "Differences_between_samples.csv")))

        dataframe.to_csv(os.path.join(working_path, f"Differences_between_samples.csv"))

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

        def my_bar_plot(x, y, **kwargs):
            colors = ['red' if val < 0 else 'green' for val in x]
            plt.barh(y=y, width=np.abs(x), color=colors)

        g = sb.FacetGrid(data=df, col='ID', height=9, aspect=0.2,
                         col_order=list(dataframe.columns), sharey=True)
        g.map(my_bar_plot, 'Difference', 'Codon')
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

    def ilustrate_patterns(self, patterns_df):
        direction = {}

        for sample in patterns_df:

            for codon in patterns_df[sample]:
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

    def plot_reference(self, RSCU_reference, working_path):
        DATA = RSCU_reference.drop(columns=['GENE CAI'])
        DATA = DATA.drop(columns=['GENOME CAI'])
        data = DATA.T
        codons = Constants.TOTAL_CODONS[:-3]
        data['Codon'] = codons
        fig, ax = plt.subplots(figsize=(6, 10))
        fig.set_figwidth(10)
        plt.subplots_adjust(left=0.3)

        norm = plt.Normalize(vmin=min(data['genome']), vmax=max(data['genome']))
        cmap = plt.get_cmap('brg')
        color = cmap(norm(data['genome']))
        ax.barh(data['Codon'], data['genome'], color=color)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.xlabel('RSCU')
        plt.ylabel('Codon')
        plt.subplots_adjust(wspace=0.3)
        plt.title(f'RSCU values for each codon in the reference gCUB')
        print("Create image: {}".format(os.path.join(working_path, f'Barplot_to_counts_reference.png')))
        plt.savefig(os.path.join(working_path, f'Barplot_to_counts_reference.png'))
        return DATA

    def plot_counts(self, counts, samples, working_path, b_make_averages_for_same_time_points):
        data = 'RSCU'
        if b_make_averages_for_same_time_points:
            time_points = [f'{self.sample.dt_sample[sample].age}' for sample in samples]
        else:
            time_points = samples

        dic_codons = OrderedDict()
        for n, dataframe in enumerate(counts):
            dic_codons[time_points[n]] = {}
            for codon in dataframe:
                if codon not in dic_codons[time_points[n]]:
                    dic_codons[time_points[n]][codon] = dataframe[codon][data]
                else:
                    dic_codons[time_points[n]][codon].append(dataframe[codon][data])

        data_values = pd.DataFrame(dic_codons, columns=time_points)
        codons = Constants.TOTAL_CODONS
        data_values['Codon'] = codons

        df = pd.melt(data_values, id_vars='Codon', value_vars=time_points, value_name='RSCU values')

        max = 0
        min = 100000

        for value in df['RSCU values']:
            if type(value) == float:
                if value > max:
                    max = value
                elif value < min:
                    min = value

        norm = TwoSlopeNorm(vcenter=max - (max / 2), vmin=min, vmax=max)
        cmap = plt.get_cmap('brg')

        def my_bar_plot(x, y, **kwargs):
            plt.barh(y=y, width=np.abs(x), color=cmap(norm(x)))

        df.rename(columns={"variable": "ID"}, inplace=True)
        g = sb.FacetGrid(data=df, col='ID', height=9, aspect=0.2,
                         col_order=time_points, sharey=True)
        g.map(my_bar_plot, 'RSCU values', 'Codon')
        g.fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', ax=g.axes, fraction=0.1,
                       shrink=0.2)
        print("Create image: {}".format(os.path.join(working_path, f'Barplot_to_counts_{data}.png')))
        plt.savefig(os.path.join(working_path, f'Barplot_to_counts_{data}.png'))

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

    ''' def PCA_analysis(self, counts, samples, working_path, **kwargs):

        ## define working path
        working_path_PCA = os.path.join(working_path, 'DEGs')
        self.utils.make_path(working_path_PCA)

        diff_expressed_genes = kwargs['genes']
        comparisons = kwargs['comparisons']
        consecutive = kwargs['consecutive']
        self.gender = kwargs['gender']
        times = [self.sample.dt_sample[sample].age for sample in samples]
        explained_variances = OrderedDict()

        for n, comparison in enumerate(comparisons):
            comparison_counts = []
            final_samples = []
            final_counts = []
            for m, time in enumerate(times):
                if int(time) in comparisons[n]:
                    comparison_counts.append(counts[m])
                    final_samples.append(samples[m])

            if f'{comparison[0]}vs{comparison[1]}' in diff_expressed_genes.keys():
                wanted_indexes = diff_expressed_genes[f'{comparison[0]}vs{comparison[1]}']
                new_counts = [dataframe.loc[dataframe.index.isin(wanted_indexes)]
                              for dataframe in comparison_counts]

                for dataframe in new_counts:
                    dataframe_copy = dataframe.copy()
                    totals = dataframe_copy.sum(axis=0)
                    dataframe_copy.loc['Total'] = totals
                    rscu = self.calculate_RSCU(dataframe_copy)
                    rscu_dataframe = pd.DataFrame(rscu, index=['RSCU'])
                    final_dataframe = pd.concat([dataframe_copy, rscu_dataframe], axis=0)
                    final_counts.append(final_dataframe)

                RSCU_dic = OrderedDict()
                for x, dataframe in enumerate(final_counts):
                    for codon in dataframe:
                        if final_samples[x] not in RSCU_dic:
                            RSCU_dic[final_samples[x]] = [dataframe[codon]['RSCU']]
                        else:
                            RSCU_dic[final_samples[x]].append(dataframe[codon]['RSCU'])

                RSCU_dataframe = pd.DataFrame.from_dict(RSCU_dic, orient='columns')
                RSCU_dataframe['Codon'] = [str(key).upper().replace('U', 'T') for key in Constants.TOTAL_CODONS]
                RSCU_dataframe.set_index('Codon', inplace=True)
                times_dataframe = [int(self.sample.dt_sample[sample].age) for sample in RSCU_dataframe.columns]

                ## save file, to REMOVE
                if consecutive:
                    print("Create RSCU_dataframe csv: {}".format(
                        os.path.join(working_path_PCA, f'ICA_Clustering_DATA_{comparison[0]}vs{comparison[1]}.csv')))
                    RSCU_dataframe.to_csv(
                        os.path.join(working_path_PCA, f'ICA_Clustering_DATA_{comparison[0]}vs{comparison[1]}.csv'),
                        index=True)
                else:
                    print("Create RSCU_dataframe csv: {}".format(
                        os.path.join(working_path_PCA, f'PCA_Clustering_DATA_{comparison[0]}vs{comparison[1]}.csv')))
                    RSCU_dataframe.to_csv(
                        os.path.join(working_path_PCA, f'PCA_Clustering_DATA_{comparison[0]}vs{comparison[1]}.csv'),
                        index=True)
                ## END save file, to REMOVE

                if consecutive:
                    RSCU_matrix = RSCU_dataframe.to_numpy()
                    ica = FastICA(n_components=2, whiten='unit-variance')
                    ica_result = ica.fit_transform(RSCU_matrix)
                    fig, ax = plt.subplots()
                    colors = plt.cm.Set1(np.linspace(0, 1, len(np.unique(times_dataframe))))
                    for i, tp in enumerate(sorted(np.unique(times_dataframe))):
                        indices = np.where(times_dataframe == tp)
                        c = colors[i]
                        ax.scatter(ica_result[indices, 0], ica_result[indices, 1], color=c, label=f'Time {tp}')
                    lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    print("Create image consecutive: {}".format(os.path.join(working_path_PCA,
                                                                             f'ICA_Clustering_consecutive_{comparison[0]}vs{comparison[1]}.png')))
                    plt.subplots_adjust(right=0.7)
                    plt.title(f'ICA analysis consecutive {comparison[0]}vs{comparison[1]}')
                    plt.savefig(os.path.join(working_path_PCA,
                                             f'ICA_Clustering_consecutive_{comparison[0]}vs{comparison[1]}.png'),
                                bbox_extra_artists=(lgd,), bbox_inches='tight')
                    componentes = ica.components_
                    # explained_variance = ica.explained_variance_
                    num_components = 2
                    num_codons_highest_weight = 10
                    for i in range(num_components):
                        components_weights = componentes[i]
                        indexes_highest_weight = np.argsort(components_weights)[::-1][:num_codons_highest_weight]
                        codons_highest_weight = list(RSCU_dataframe.index[indexes_highest_weight])
                        highest_weights = [float(components_weights[n]) for n in indexes_highest_weight]
                        df = pd.DataFrame(codons_highest_weight)
                        highest_weights_df = pd.DataFrame(highest_weights)
                        final_df = pd.concat([df.T, highest_weights_df.T])

                        final_df.to_excel(os.path.join(working_path_PCA,
                                                       f'Major_influence_codons_consecutive_{comparison[0]}vs{comparison[1]}_ICA.xlsx'),
                                          index=False)
                else:
                    # PCA analysis
                    print("Stat PCA analysis")
                    pca = PCA(n_components=2)
                    pca_result = pca.fit_transform(RSCU_dataframe.transpose())
                    fig, ax = plt.subplots()
                    colors = plt.cm.Set1(np.linspace(0, 1, len(np.unique(times_dataframe))))
                    for i, tp in enumerate(sorted(np.unique(times_dataframe))):
                        a = np.where(times_dataframe == tp)
                        c = colors[i]
                        ax.scatter(pca_result[:, 0][a], pca_result[:, 1][a], color=c, label=f'Time {tp}')
                    lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    plt.subplots_adjust(right=0.7)

                    plt.title(f'PCA analysis non consecutive {comparison[0]}vs{comparison[1]}')
                    print("Create image not consecutive: {}".format(
                        os.path.join(working_path_PCA,
                                     f'PCA_analysis_non_consecutive_{comparison[0]}vs{comparison[1]}.png')))
                    plt.savefig(os.path.join(working_path_PCA,
                                             f'PCA_analysis_non_consecutive_{comparison[0]}vs{comparison[1]}.png'),
                                bbox_extra_artists=(lgd,), bbox_inches='tight')
                    componentes = pca.components_
                    explained_variance = pca.explained_variance_ratio_
                    if f'{comparison[0]}vs{comparison[1]}' not in explained_variances:
                        explained_variances[f'{comparison[0]}vs{comparison[1]}'] = list(explained_variance)
                    print(
                        f'explained variance to comparison {comparison[0]}vs{comparison[1]}: C1-{explained_variance[0]}, C2-{explained_variance[1]}')
                    num_components = 2
                    num_codons_highest_weight = 10
                    for i in range(num_components):
                        components_weights = componentes[i]
                        indexes_highest_weight = np.argsort(components_weights)[::-1][:num_codons_highest_weight]
                        codons_highest_weight = list(RSCU_dataframe.index[indexes_highest_weight])
                        # To present weights in per cent
                        # highest_weights = [float(weight)*100 for n, weight in list(components_weights) if n in indexes_highest_weight]
                        highest_weights = [float(components_weights[n]) for n in indexes_highest_weight]
                        df = pd.DataFrame(codons_highest_weight)
                        highest_weights_df = pd.DataFrame(highest_weights)
                        final_df = pd.concat([df.T, highest_weights_df.T])
                        final_df.to_excel(os.path.join(working_path_PCA,
                                                       f'Major_influence_codons_non_consecutive_{comparison[0]}vs{comparison[1]}_PCA.xlsx'),
                                          index=False)
        if len(explained_variances) != 0:
            dataframe = pd.DataFrame(explained_variances)
            dataframe.to_excel(os.path.join(working_path_PCA,
                                            f'Explained_variances.xlsx'),
                               index=False)'''

    def test_X2(self, counts, samples, comparisons, working_path, liver, consecutive, gender):

        utils.make_path(os.path.join(working_path, gender, 'Chi2'))
        significant_aminoacids = []
        path = os.path.join(working_path, gender, 'Chi2')
        totals = OrderedDict()
        time = [int(self.sample.dt_sample[sample].age) for sample in samples]
        n_comparisons = len(comparisons)
        alpha = 0.05/ (19 * n_comparisons)
        print(n_comparisons)
        print(alpha)
        print("time: ", time)
        for i, dataframe in enumerate(counts):
            line = dataframe.iloc[-2]
            totals[str(time[i])] = line
        totals_dataframe = pd.DataFrame(totals)
        #totals_dataframe.to_csv(os.path.join(working_path, 'Totals_dataframe.csv'))

        aminoacids = []
        codons = []

        for aminoacido, lista_codons in Constants.ordered_codons.items():
            aminoacids.extend([aminoacido] * len(lista_codons))
            codons.extend(lista_codons)

        aminoacid_codons_groups = pd.DataFrame({'Aminoacid': aminoacids, 'Codon': codons})

        merged_dataframe = pd.merge(totals_dataframe, aminoacid_codons_groups, left_index=True, right_on='Codon')
        dataframe = merged_dataframe.groupby('Aminoacid')
        df_final = None
        print('Chi2 test to each comparison')
        for comparison in comparisons:
            time_point_1 = str(comparison[0])
            time_point_2 = str(comparison[1])
            results = []
            for aminoacid, codon_group in dataframe:
                #print("################## aminoacid: ", aminoacid)
                if len(codon_group) > 1:
                    codons_to_compare = codon_group['Codon'].tolist()

                    time_points_to_compare = [str(time_point_1), str(time_point_2)]
                    codon_counts_to_compare = codon_group[['Codon'] + time_points_to_compare]

                    contingency_table = []

                    for codon in codons_to_compare:
                        codon_counts_sample_1 = codon_counts_to_compare[codon_counts_to_compare['Codon'] == codon][
                            str(time_point_1)].values
                        codon_counts_sample_2 = codon_counts_to_compare[codon_counts_to_compare['Codon'] == codon][
                            str(time_point_2)].values

                        # Contigency table
                        if len(contingency_table) == 0:
                            contingency_table = np.array([[codon_counts_sample_1[0], codon_counts_sample_2[0]]])
                        else:
                            contingency_table = np.append(contingency_table,
                                                          [[codon_counts_sample_1[0], codon_counts_sample_2[0]]], 0)
                    # Chi2 test
                    chi2, p_value, dof, expected = chi2_contingency(contingency_table)
                    if p_value > alpha:
                        results.append((aminoacid, chi2, p_value, 'Non significant'))
                    else:
                        results.append((aminoacid, chi2, p_value, 'Significant'))
            results_dataframe = pd.DataFrame(results,
                                             columns=[f'Aminoacid_{comparison[0]}vs{comparison[1]}', 'Chi2', 'p_value',
                                                      f'Significance: {alpha}'])
            significant_aminoacids.append(results_dataframe[
                results_dataframe[f'Significance: {0.05 / (19 * n_comparisons)}'] == 'Significant'])
            df_final = pd.concat([df_final, results_dataframe], axis=1)


        if liver:
            if consecutive:
                for n, dataframe in enumerate(significant_aminoacids):
                    dataframe.to_csv(os.path.join(path,
                                                  f'Significant_aminoacids_to_{comparisons[n]}_consecutive.csv'),
                                     index=False)
                print(
                    "File Chi2 test results: " + str(os.path.join(path, f'Chi2_test_liver_consecutive.csv')))
                df_final.to_csv(os.path.join(path, f'Chi2_test_liver_consecutive.csv'), index=False)

        else:
            if consecutive:
                for n, dataframe in enumerate(significant_aminoacids):
                    dataframe.to_csv(os.path.join(path,
                                                  f'Significant_aminoacids_to_{comparisons[n]}_consecutive.csv'),
                                     index=False)
                print(
                    "File Chi2 test results: " + str(os.path.join(path, f'Chi2_test_brain_consecutive.csv')))
                df_final.to_csv(os.path.join(path, f'Chi2_test_brain_consecutive.csv'), index=False)

    def test_X2_with_diff_expressed_genes(self, counts, samples, comparisons, working_path, liver, consecutive, genes, gender):
        path = os.path.join(working_path, gender, 'Chi2')
        times = [int(self.sample.dt_sample[sample].age) for sample in samples]
        significant_aminoacids = []
        n_comparisons = 0
        print('Chi2 test to each comparison with DEGS')
        for comparison in comparisons:
            if f'{comparison[0]}vs{comparison[1]}' in genes.keys():
                n_comparisons += 1
        df_final = None
        alpha = 0.05 / (19 * n_comparisons)
        print(n_comparisons)
        print(alpha)
        for n, comparison in enumerate(comparisons):
            comparison_counts = []
            for m, time in enumerate(times):
                if time in comparisons[n]:
                    comparison_counts.append(counts[m])

            if f'{comparison[0]}vs{comparison[1]}' in genes.keys():
                wanted_indexes = genes[f'{comparison[0]}vs{comparison[1]}']
                new_counts = [dataframe.loc[dataframe.index.isin(wanted_indexes)] for dataframe in comparison_counts]
                counts_with_totals = []
                for dataframe in new_counts:
                    dataframe_copy = dataframe.copy()
                    totals = dataframe_copy.sum(axis=0)
                    dataframe_copy.loc['Total'] = totals
                    counts_with_totals.append(dataframe_copy)

                totals_dict = OrderedDict()
                for t, time in enumerate(comparison):
                    #print("time: ", time)
                    for i, dataframe in enumerate(counts_with_totals):
                        if i == t:
                            line = dataframe.iloc[-1]
                            totals_dict[str(time)] = line
                totals_dataframe = pd.DataFrame(totals_dict)
                aminoacidos = []
                codons = []
                for aminoacido, lista_codons in Constants.ordered_codons.items():
                    aminoacidos.extend([aminoacido] * len(lista_codons))
                    codons.extend(lista_codons)
                aminoacid_codons_groups = pd.DataFrame({'Aminoacid': aminoacidos, 'Codon': codons})
                merged_dataframe = pd.merge(totals_dataframe, aminoacid_codons_groups, left_index=True,
                                            right_on='Codon')
                dataframe = merged_dataframe.groupby('Aminoacid')
                time_point_1 = str(comparison[0])
                time_point_2 = str(comparison[1])
                results = []
                for aminoacid, codon_group in dataframe:
                    #print("##################aminoacid: ", aminoacid)
                    if len(codon_group) > 1:
                        codons_to_compare = codon_group['Codon'].tolist()
                        time_points_to_compare = [str(time_point_1), str(time_point_2)]
                        codon_counts_to_compare = codon_group[['Codon'] + time_points_to_compare]
                        contingency_table = []
                        for codon in codons_to_compare:
                            codon_counts_sample_1 = \
                                codon_counts_to_compare[codon_counts_to_compare['Codon'] == codon][
                                    str(time_point_1)].values
                            codon_counts_sample_2 = \
                                codon_counts_to_compare[codon_counts_to_compare['Codon'] == codon][
                                    str(time_point_2)].values
                            # Contigency table
                            if len(contingency_table) == 0:
                                contingency_table = np.array(
                                    [[codon_counts_sample_1[0], codon_counts_sample_2[0]]])
                            else:
                                contingency_table = np.append(contingency_table, [
                                    [codon_counts_sample_1[0], codon_counts_sample_2[0]]], 0)
                        # Chi2 test
                        contingency_table += 0.5
                        chi2, p_value, dof, expected = chi2_contingency(contingency_table)
                        if p_value > alpha:
                            results.append((aminoacid, chi2, p_value, 'Non significant'))
                        else:
                            results.append((aminoacid, chi2, p_value, 'Significant'))
                        #print(" chi2: ", chi2, "   p_value: ", p_value)
                results_dataframe = pd.DataFrame(results,
                                                 columns=[f'Aminoacid_{comparison[0]}vs{comparison[1]}', 'Chi2',
                                                          'p_value',
                                                          f'Significance: {alpha}'])
                significant_aminoacids.append(results_dataframe[
                    results_dataframe[f'Significance: {0.05 / (19 * n_comparisons)}'] == 'Significant'])
                df_final = pd.concat([df_final, results_dataframe], axis=1)

                if liver:
                    if consecutive:
                        for n, dataframe in enumerate(significant_aminoacids):

                            dataframe.to_csv(os.path.join(path,
                                                         f'Significant_aminoacids_to_{comparison[0]}_{comparison[1]}_DEGs_consecutive.csv'),
                                            index=False)
                        print(
                        "File Chi2 test results: " + str(os.path.join(path, f'Chi2_test_liver_DEGS_consecutive.csv')))
                        df_final.to_csv(os.path.join(path,
                                                              f'Chi2_test_liver_DEGS_consecutive.csv'),
                                                 index=False)
                    else:
                        for n, dataframe in enumerate(significant_aminoacids):

                            dataframe.to_csv(os.path.join(path,
                                                         f'Significant_aminoacids_to_{comparison[0]}_{comparison[1]}_DEGs_non_consecutive.csv'),
                                            index=False)
                        print(
                        "File Chi2 test results: " + str(os.path.join(path, f'Chi2_test_liver_DEGS_non_consecutive.csv')))
                        df_final.to_csv(os.path.join(path,
                                                              f'Chi2_test_liver_DEGS_non_consecutive.csv'),
                                                 index=False)
                else:
                    if consecutive:
                        for n, dataframe in enumerate(significant_aminoacids):

                            dataframe.to_csv(os.path.join(path,
                                                         f'Significant_aminoacids_to_{comparison[0]}_{comparison[1]}_DEGs_consecutive.csv'),
                                            index=False)
                        print(
                        "File Chi2 test results: " + str(os.path.join(path, f'Chi2_test_brain_DEGS_consecutive.csv')))
                        df_final.to_csv(
                            os.path.join(path,
                                         f'Chi2_test_brain_DEGS_consecutive.csv'),
                            index=False)
                    else:
                        for n, dataframe in enumerate(significant_aminoacids):

                            dataframe.to_csv(os.path.join(path,
                                                         f'Significant_aminoacids_to_{comparison[0]}_{comparison[1]}_DEGs_non_consecutive.csv'),
                                            index=False)
                        print(
                        "File Chi2 test results: " + str(os.path.join(path, f'Chi2_test_brain_DEGS_non_consecutive.csv')))
                        df_final.to_csv(
                            os.path.join(path,
                                         f'Chi2_test_brain_DEGS_non_consecutive.csv'),
                            index=False)

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
