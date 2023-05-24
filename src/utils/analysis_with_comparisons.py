from constants.constants import Constants
import socket
import os
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.cm import ScalarMappable
import seaborn as sb
from sklearn.decomposition import PCA



class Comparison(object):
    def __init__(self, counts, samples, gender, liver):
        self.working_path = None
        self.consecutive = False
        self.comparison = None
        self.counts = counts
        self.samples = samples
        self.gender = gender
        self.liver = liver
        self.final_times = []
        self.differentially_expressed_genes, self.counts_to_degs, self.differences_abs, self.differences = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()  # need to be ordered


        if self.consecutive:
            self.times = ['27vs3', '3vs6', '6vs9', '9vs12', '12vs15', '15vs18', '18vs21', '21vs24', '24vs27']
        else:
            self.times = ['3vs6', '3vs9', '3vs12', '3vs15', '3vs18', '3vs21', '3vs24', '3vs27']
            self.comparisons = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8)]
        if socket.gethostname() == "cs-nb0008":  # test computer name

            self.base_path = "/home/projects/ua/master/codon_usage"
        else:
            self.base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"

        self.get_genes()

    def get_genes(self):
        for time in self.times:
            if self.liver:
                with open(os.path.join(self.base_path, f'genes_liver_sig_{time}.csv')) as genes:
                    print( f'genes_liver_sig_{time}.csv')
                    for line in genes.readlines():
                        if 'genes' not in line:
                            line = line.strip()
                            line = line.replace('"', '')
                            if time not in self.differentially_expressed_genes:
                                self.differentially_expressed_genes[time] = [line]
                            else:
                                self.differentially_expressed_genes[time].append(line)
            else:
                with open(os.path.join(self.base_path, f'genes_brain_sig_{time}.csv')) as genes:
                    for line in genes.readlines():
                        if 'genes' not in line:
                            line = line.strip()
                            line = line.replace('"', '')
                            if time not in self.differentially_expressed_genes:
                                self.differentially_expressed_genes[time] = [line]
                            else:
                                self.differentially_expressed_genes[time].append(line)

        self.compare_timepoints()

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

    def compare_timepoints(self):
        if self.consecutive:
            self.comparison = 'consecutive_comparisons'
            final_counts = []
            totals = []
            rscu = []
            rscu_dataframes = []
            final_dataframes = []

            indexes_to_remove = []
            for n, time in enumerate(self.times):
                if time in self.differentially_expressed_genes:

                    self.final_times.append(time)
                    indices_desejados = self.differentially_expressed_genes[time]
                    final_counts.append((self.counts[n-1].loc[self.counts[n-1].index.isin(indices_desejados)],
                                         self.counts[n].loc[self.counts[n].index.isin(indices_desejados)]))
                else:
                    indexes_to_remove.append(n)

            for m, comparison in enumerate(final_counts):
                if m in indexes_to_remove:
                    continue
                else:
                    totals.append((comparison[0].sum(axis=0).T, comparison[1].sum(axis=0).T))
                    comparison_copy_0 = comparison[0].copy()  # Cópia do DataFrame comparison[0]
                    comparison_copy_1 = comparison[1].copy()  # Cópia do DataFrame comparison[1]
                    comparison_copy_0.loc['Total'] = totals[m][0]
                    comparison_copy_1.loc['Total'] = totals[m][1]
                    rscu.append((self.calculate_RSCU(comparison_copy_0), self.calculate_RSCU(comparison_copy_1)))
                    rscu_dataframes.append(
                        (pd.DataFrame(rscu[m][0], index=['RSCU']), pd.DataFrame(rscu[m][1], index=['RSCU'])))
                    final_dataframes.append((pd.concat([comparison_copy_0, rscu_dataframes[m][0]], axis=0),
                                             pd.concat([comparison_copy_1, rscu_dataframes[m][1]], axis=0)))
                    self.counts_to_degs[self.final_times[m]] = final_dataframes[m]
#print(self.counts_to_degs)

            self.plot_differences()
            #self.PCA_analysis()
        else:
            self.comparison = 'fixed_comparisons'
            final_counts = []
            totals = []
            rscu = []
            rscu_dataframes = []
            final_dataframes = []

            indexes_to_remove = []
            for n in range(1, len(self.times)+1):
                if self.times[n-1] in self.differentially_expressed_genes:
                    self.final_times.append(self.times[n-1])
                    indices_desejados = self.differentially_expressed_genes[self.times[n-1]]
                    final_counts.append((self.counts[0].loc[self.counts[0].index.isin(indices_desejados)],
                                         self.counts[n].loc[self.counts[n].index.isin(indices_desejados)]))
                else:
                    indexes_to_remove.append(n-1)

            for m, comparison in enumerate(final_counts):
                if m in indexes_to_remove:
                    continue
                else:
                    totals.append((comparison[0].sum(axis=0).T, comparison[1].sum(axis=0).T))
                    comparison_copy_0 = comparison[0].copy()  # Cópia do DataFrame comparison[0]
                    comparison_copy_1 = comparison[1].copy()  # Cópia do DataFrame comparison[1]
                    comparison_copy_0.loc['Total'] = totals[m][0]
                    comparison_copy_1.loc['Total'] = totals[m][1]
                    rscu.append((self.calculate_RSCU(comparison_copy_0), self.calculate_RSCU(comparison_copy_1)))
                    rscu_dataframes.append(
                        (pd.DataFrame(rscu[m][0], index=['RSCU']), pd.DataFrame(rscu[m][1], index=['RSCU'])))
                    final_dataframes.append((pd.concat([comparison_copy_0, rscu_dataframes[m][0]], axis=0),
                                             pd.concat([comparison_copy_1, rscu_dataframes[m][1]], axis=0)))
                    self.counts_to_degs[self.final_times[m]] = final_dataframes[m]
            # print(self.counts_to_degs)

            self.plot_differences()
            #self.PCA_analysis()

    def plot_differences(self):
        working_path = os.path.join(self.base_path, 'mouse', 'liver' if self.liver else 'brain', 'DEGs')
        for key in self.counts_to_degs:
            for dataframe in self.counts_to_degs[key]:
                self.differences_abs[key] = {}
                self.differences[key] = {}
                for value in dataframe:

                    if value not in self.differences[key]:

                        self.differences[key][
                            value] = self.counts_to_degs[key][1][value]['RSCU'] - self.counts_to_degs[key][0][value]['RSCU']
                        self.differences_abs[key][
                            value] = self.counts_to_degs[key][1][value]['RSCU'] - self.counts_to_degs[key][0][value]['RSCU']

                    else:

                        self.differences[key][
                            value] += self.counts_to_degs[key][1][value]['RSCU'] - self.counts_to_degs[key][0][value]['RSCU']
                        self.differences_abs[key][
                            value] += self.counts_to_degs[key][1][value]['RSCU'] - self.counts_to_degs[key][0][value]['RSCU']


                ### save differences
        dataframe_dif = pd.DataFrame(self.differences)

        print(
                "File with differences: " + str(os.path.join(working_path, f"Differences_between_time_points_{self.comparison}.csv")))

        dataframe_dif.to_csv(os.path.join(working_path, f"Differences_between_time_points_{self.comparison}.csv"))


            ### start making chart
        dataframe_abs = pd.DataFrame(self.differences_abs)

        dataframe_dif['Codon'] = [f'{str(key).upper().replace("U", "T")}_{value}' for key, value in
                                      Constants.codons_per_aminoacid.items()]
        columns = [time for time in self.final_times]
        #columns = ['27vs3', '3vs6', '6vs9', '9vs12', '12vs15']
        df = pd.melt(dataframe_dif, id_vars='Codon', value_vars=columns, value_name='Difference')
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
                             col_order=list(dataframe_dif.columns), sharey=True)
        g.map(my_bar_plot, 'Difference', 'Codon')

        print("Create image: {}".format(os.path.join(working_path, f'Barplot_to_differences_RSCU_DEGs_{self.comparison}.png')))

        plt.savefig(os.path.join(working_path, f'Barplot_to_differences_RSCU_DEGs_{self.comparison}.png'))

        return df

    def PCA_analysis(self):
        RSCU_dic = OrderedDict()
        time_points = self.final_times
        for tuple in self.counts:
            for n, dataframe in enumerate(tuple):
                for codon in dataframe:
                    if self.final_times[n] not in RSCU_dic:
                        RSCU_dic[self.final_times[n]] = [dataframe[codon]['RSCU']]
                    else:
                        RSCU_dic[self.final_times[n]].append(dataframe[codon]['RSCU'])

        RSCU_dataframe = pd.DataFrame.from_dict(RSCU_dic, orient='columns')

        RSCU_dataframe['Codon'] = [str(key).upper().replace('U', 'T') for key in Constants.TOTAL_CODONS]
        RSCU_dataframe.set_index('Codon', inplace=True)

        times = [time for time in RSCU_dataframe.columns]
        RSCU_dataframe = RSCU_dataframe.transpose()
        RSCU_dataframe["time"] = times
        RSCU_dataframe = RSCU_dataframe.transpose()

        # Obtain time points
        time_points = RSCU_dataframe.iloc[-1, :].values
        RSCU_dataframe.drop(RSCU_dataframe.tail(1).index, inplace=True)

        # PCA analysis
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(RSCU_dataframe.transpose())
        fig, ax = plt.subplots()
        colors = plt.cm.Set1(np.linspace(0, 1, len(np.unique(time_points))))
        for i, tp in enumerate(sorted(np.unique(time_points))):
            samples = np.where(time_points == tp)
            c = colors[i]
            ax.scatter(pca_result[:, 0][samples], pca_result[:, 1][samples], color=c, label=f'Time {tp}')
        lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(right=0.7)

        plt.title(f'PCA to comparisons')
        print("Create image: {}".format(os.path.join(self.working_path, f'PCA_analysis_{self.comparison}.png')))
        plt.savefig(os.path.join(self.working_path, f'PCA_analysis_{self.comparison}.png'), bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
