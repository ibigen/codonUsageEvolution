'''
Created on 17/09/2022

@author: mmp
'''
import unittest
import os
from codon_usage_evolution import read_genome, save_table
import filecmp
from utils.utils import Utils
from utils.expression import Expression, Sample, Tissue


class Test(unittest.TestCase):
    utils = Utils()

    def setUp(self):
        self.baseDirectory = os.path.dirname(os.path.abspath(__file__))
        print(self.baseDirectory)

    def tearDown(self):
        pass

    def sum(self, a, b):
        """ sum two numbers, only for example"""
        return a + b

    def test_sum(self):
        self.assertEqual(5, self.sum(3, 2))
        self.assertNotEqual(6, self.sum(3, 2))

    def test_read_fasta(self):
        ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
        self.assertTrue(os.path.exists(ecoli_fasta))

        # call method
        dataframe_counts, dataframe_RSCU_CAI, stats = read_genome(ecoli_fasta)

        ### stats
        self.assertEqual(1, stats.count_divisible_3)
        self.assertEqual(12, stats.count_pass)

        # test if this gene is inside data frame
        self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", dataframe_counts.index)
        self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", dataframe_RSCU_CAI.index)
        self.assertIn("lcl|NC_000913.3_cds_YP_009518733.1_12", dataframe_counts.index)
        self.assertIn("lcl|NC_000913.3_cds_YP_009518733.1_12", dataframe_RSCU_CAI.index)
        self.assertIn('lcl|NC_000913.3_cds_NP_414550.1_9', dataframe_counts.index)
        self.assertIn('lcl|NC_000913.3_cds_NP_414550.1_9', dataframe_RSCU_CAI.index)
        self.assertIn('lcl|NC_000913.3_cds_NP_414546.1_5', dataframe_counts.index)
        self.assertIn('lcl|NC_000913.3_cds_NP_414546.1_5', dataframe_RSCU_CAI.index)
        self.assertEqual(1, dataframe_counts['TGA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
        self.assertEqual(5, dataframe_counts['GTT']['lcl|NC_000913.3_cds_NP_414545.1_4'])
        self.assertEqual(3, dataframe_counts['ACT']['lcl|NC_000913.3_cds_NP_414550.1_9'])
        self.assertEqual(1, dataframe_counts['GTG']['lcl|NC_000913.3_cds_NP_414546.1_5'])
        self.assertTrue(2.0 == dataframe_RSCU_CAI['GGT']['lcl|NC_000913.3_cds_NP_414542.1_1'])
        self.assertEqual(1.22, round(dataframe_RSCU_CAI['TTC']['lcl|NC_000913.3_cds_NP_414551.1_10'], 2))
        self.assertTrue(1.6 == dataframe_RSCU_CAI['TTT']['lcl|NC_000913.3_cds_NP_414550.1_9'])
        self.assertEqual(0.4, dataframe_RSCU_CAI['TTC']['lcl|NC_000913.3_cds_NP_414550.1_9'])
        self.assertEqual(13, len(dataframe_counts))
        self.assertEqual(13, len(dataframe_RSCU_CAI))

    def test_tables(self):
        ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
        expected_result_RSCU_CAI = os.path.join(self.baseDirectory, "files/tables/table_RSCU_CAI_test.csv")
        self.assertTrue(os.path.exists(expected_result_RSCU_CAI))
        dataframe_counts, dataframe_RSCU_CAI, stats = read_genome(ecoli_fasta)
        csv_result_RSCU_CAI = self.utils.get_temp_file("RSCU_and_CAI_to_test", ".csv")
        save_table(dataframe_RSCU_CAI, csv_result_RSCU_CAI)
        self.assertTrue(filecmp.cmp(expected_result_RSCU_CAI, csv_result_RSCU_CAI))

        expected_result_counts = os.path.join(self.baseDirectory, "files/tables/table_counts_test.csv")
        csv_result_counts = self.utils.get_temp_file("counts_to_test", ".csv")
        save_table(dataframe_counts, csv_result_counts)
        self.assertTrue(filecmp.cmp(expected_result_counts, csv_result_counts))

        # remove temp files
        self.utils.remove_file(csv_result_RSCU_CAI)
        self.utils.remove_file(csv_result_counts)

    def test_expression(self):
        self.expressionDirectory = os.path.dirname(os.path.abspath(__file__))
        file_information = os.path.join(self.expressionDirectory, "files/expression/E.coli_information.txt")
        file_expression = os.path.join(self.expressionDirectory, "files/expression/E.coli_expression.txt")
        self.assertTrue(os.path.exists(file_information))
        self.assertTrue(os.path.exists(file_expression))

        # Call method
        expression = Expression(file_information, file_expression)

        self.assertEqual(6, expression.sample.get_number_sample())
        self.assertEqual(6, expression.get_number_sample())
        self.assertEqual(12, expression.get_number_gene("A18_384Bulk_Plate1_S18"))
        self.assertIn("A18_384Bulk_Plate1_S18", expression.sample.dt_sample)
        self.assertEqual('Male', expression.sample.dt_sample["A9_384Bulk_Plate1_S9"].sex)
        self.assertEqual(3, float(expression.sample.dt_sample["A9_384Bulk_Plate1_S9"].age))
        self.assertEqual('Female', expression.sample.dt_sample["A18_384Bulk_Plate1_S18"].sex)
        self.assertEqual('Brain', expression.sample.dt_sample["A18_384Bulk_Plate1_S18"].tissue)

        self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", expression.sample.dt_sample["A9_384Bulk_Plate1_S9"].dt_gene)
        self.assertIn("lcl|NC_000913.3_cds_NP_414551.1_10", expression.sample.dt_sample["A9_384Bulk_Plate1_S9"].dt_gene)
        self.assertEqual(2109.15514707196, expression.sample.dt_sample["A9_384Bulk_Plate1_S9"].dt_gene[
            'lcl|NC_000913.3_cds_NP_414542.1_1'])
        self.assertEqual(13.8749030600905, expression.sample.dt_sample['A18_384Bulk_Plate1_S18'].dt_gene[
            'lcl|NC_000913.3_cds_NP_414552.1_11'])
        self.assertEqual(0.845691718954275, expression.sample.dt_sample['A9_384Bulk_Plate1_S9'].dt_gene[
            'lcl|NC_000913.3_cds_YP_009518733.1_12'])
        self.assertEqual(8.53840188313262, expression.sample.dt_sample['A18_384Bulk_Plate1_S18'].dt_gene[
            'lcl|NC_000913.3_cds_YP_009518733.1_12'])

        self.assertEqual(12, len(expression.most_differentially_expressed_genes('E20_384Bulk_Plate1_S116',
                                                                                'A9_384Bulk_Plate1_S9')))
        self.assertEqual(363.88953212254, expression.most_dif_expressed['lcl|NC_000913.3_cds_NP_414542.1_1'])
        self.assertEqual(0.115760470543813, expression.most_dif_expressed['lcl|NC_000913.3_cds_NP_414551.1_10'])
        self.assertEqual(12.602181753894525,
                         expression.most_differentially_expressed_genes('A20_384Bulk_Plate2_S20',
                                                                        'A9_384Bulk_Plate1_S9')
                         ['lcl|NC_000913.3_cds_NP_414551.1_10'])


        ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
        dataframe_counts, dataframe_RSCU_CAI, stats = read_genome(ecoli_fasta)
        
        self.assertEqual(2109.15514707196, expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                             dataframe_counts.to_dict(orient='index'))['AAA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
        self.assertEqual(34.15360753253048, expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                              dataframe_counts.to_dict(orient='index'))['TTT']['lcl|NC_000913.3_cds_YP_009518733.1_12'])
        self.assertEqual(3.3827668758171, expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                              dataframe_counts.to_dict(orient='index'))[
            'TTT']['lcl|NC_000913.3_cds_YP_009518733.1_12'])
        self.assertEqual(2359.80082045078, expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                              dataframe_counts.to_dict(orient='index'))[
            'AAA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
        self.assertEqual(13215.624492098465, expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                             dataframe_counts.to_dict(orient='index'))[
            'AAA']['Total'])
        self.assertEqual(11088.182145483095, expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                               dataframe_counts.to_dict(
                                                                                   orient='index'))[
            'AAA']['Total'])
        self.assertEqual(11181.037265962166, expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                              dataframe_counts.to_dict(orient='index'))[
            'TTT']['Total'])
        self.assertEqual(12024.890551810831, expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                               dataframe_counts.to_dict(
                                                                                   orient='index'))[
            'TTT']['Total'])

        self.assertEqual(-843.8532858486651, expression.compare_T0_T1(expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                dataframe_counts.to_dict(orient='index')), expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                dataframe_counts.to_dict(orient='index')))['Total']['TTT'])
        self.assertEqual(843.8532858486651,
                         expression.compare_T0_T1(expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                                    dataframe_counts.to_dict(
                                                                                        orient='index')),
                                                  expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                                    dataframe_counts.to_dict(
                                                                                        orient='index')))['Total']['TTT'])
        self.assertEqual(-2127.4423466153694,
                         expression.compare_T0_T1(expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                                    dataframe_counts.to_dict(
                                                                                        orient='index')),
                                                  expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                                    dataframe_counts.to_dict(
                                                                                        orient='index')))['Total'][
                             'AAA'])
        self.assertEqual(2127.4423466153694,
                         expression.compare_T0_T1(expression.counts_with_expression('A18_384Bulk_Plate1_S18',
                                                                                    dataframe_counts.to_dict(
                                                                                        orient='index')),
                                                  expression.counts_with_expression('A9_384Bulk_Plate1_S9',
                                                                                    dataframe_counts.to_dict(
                                                                                        orient='index')))['Total'][
                             'AAA'])




if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_read_fasta']
    unittest.main()
