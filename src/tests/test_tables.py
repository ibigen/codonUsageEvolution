'''
Created on 17/09/2022

@author: mmp
'''
import unittest
import os
from codon_usage_evolution import read_genome, save_table
import filecmp
from utils.utils import Utils
from count_sequences import CountSequences


class Test(unittest.TestCase):
	utils = Utils()

	def setUp(self):
		self.baseDirectory = os.path.dirname(os.path.abspath(__file__))

	def tearDown(self):
		pass

	def sum(self, a, b):
		""" sum two numbers, only for example"""
		return a + b

	def test_sum(self):
		self.assertEqual(5, self.sum(3, 2))
		self.assertNotEqual(6, self.sum(3, 2))

	def test_read_fasta(self):
		counts = CountSequences()
		ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
		self.assertTrue(os.path.exists(ecoli_fasta))

		# call method
		dataframe_counts, dataframe_RSCU_CAI, dataframe_CAI = read_genome(ecoli_fasta)


		# test if this gene is inside data frame
		self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", dataframe_counts.index)
		self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", dataframe_CAI.index)
		self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", dataframe_RSCU_CAI.index)
		self.assertIn("lcl|NC_000913.3_cds_YP_009518733.1_12", dataframe_counts.index)
		self.assertIn("lcl|NC_000913.3_cds_YP_009518733.1_12", dataframe_CAI.index)
		self.assertIn("lcl|NC_000913.3_cds_YP_009518733.1_12", dataframe_RSCU_CAI.index)
		self.assertIn('lcl|NC_000913.3_cds_NP_414550.1_9', dataframe_counts.index)
		self.assertIn('lcl|NC_000913.3_cds_NP_414550.1_9', dataframe_CAI.index)
		self.assertIn('lcl|NC_000913.3_cds_NP_414550.1_9', dataframe_RSCU_CAI.index)
		self.assertIn('lcl|NC_000913.3_cds_NP_414546.1_5', dataframe_counts.index)
		self.assertIn('lcl|NC_000913.3_cds_NP_414546.1_5', dataframe_CAI.index)
		self.assertIn('lcl|NC_000913.3_cds_NP_414546.1_5', dataframe_RSCU_CAI.index)
		self.assertEqual(1, dataframe_counts['TGA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
		self.assertEqual(5, dataframe_counts['GTT']['lcl|NC_000913.3_cds_NP_414545.1_4'])
		self.assertEqual(3, dataframe_counts['ACT']['lcl|NC_000913.3_cds_NP_414550.1_9'])
		self.assertEqual(1, dataframe_counts['GTG']['lcl|NC_000913.3_cds_NP_414546.1_5'])
		self.assertTrue(1.3333333333333333 == dataframe_RSCU_CAI['AAA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
		self.assertTrue(2.0 == dataframe_RSCU_CAI['GGT']['lcl|NC_000913.3_cds_NP_414542.1_1'])
		self.assertTrue(0.6666666666666666 == dataframe_RSCU_CAI['GAG']['lcl|NC_000913.3_cds_YP_009518733.1_12'])
		self.assertTrue((0.8295403880420882 == dataframe_CAI[0]['lcl|NC_000913.3_cds_NP_414542.1_1']))
		self.assertEqual(12, len(dataframe_counts))
		self.assertEqual(12, len(dataframe_RSCU_CAI))
		self.assertEqual(12, len(dataframe_CAI))



	def test_tables(self):
		ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
		expected_result_RSCU_CAI = os.path.join(self.baseDirectory, "files/tables/table_RSCU_CAI_test.csv")
		self.assertTrue(os.path.exists(expected_result_RSCU_CAI))
		dataframe_counts, dataframe_RSCU_CAI, dataframe_CAI = read_genome(ecoli_fasta)
		csv_result_RSCU_CAI = self.utils.get_temp_file("RSCU_and_CAI_to_test", ".csv")
		save_table(dataframe_RSCU_CAI, csv_result_RSCU_CAI)
		print(expected_result_RSCU_CAI, csv_result_RSCU_CAI)
		self.assertTrue(filecmp.cmp(expected_result_RSCU_CAI, csv_result_RSCU_CAI))

		expected_result_counts = os.path.join(self.baseDirectory, "files/tables/table_counts_test.csv")
		csv_result_counts = self.utils.get_temp_file("counts_to_test", ".csv")
		save_table(dataframe_counts, csv_result_counts)
		self.assertTrue(filecmp.cmp(expected_result_counts, csv_result_counts))

		expected_result_CAI = os.path.join(self.baseDirectory, "files/tables/table_CAI_test.csv")
		csv_result_CAI = self.utils.get_temp_file("CAI_to_test", ".csv")
		save_table(dataframe_CAI, csv_result_CAI)
		self.assertTrue(filecmp.cmp(expected_result_CAI, csv_result_CAI))

		## remove temp files
		self.utils.remove_file(csv_result_RSCU_CAI)
		self.utils.remove_file(csv_result_counts)
		
if __name__ == "__main__":
	# import sys;sys.argv = ['', 'Test.test_read_fasta']
	unittest.main()
