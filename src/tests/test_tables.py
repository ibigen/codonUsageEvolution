'''
Created on 17/09/2022

@author: mmp
'''
import unittest
import os
from codon_usage_evolution import read_genome, save_table
import filecmp
from utils.utils import Utils


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
		ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
		self.assertTrue(os.path.exists(ecoli_fasta))

		# call method
		dataframe_counts, dataframe_RSCU_CAI, dataframe_CAI = read_genome(ecoli_fasta)


		# test if this gene is inside data frame
		self.assertIn("lcl|NC_000913.3_cds_NP_414542.1_1", dataframe_counts.index)
		self.assertEqual(1, dataframe_counts['AAA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
		self.assertTrue(1.3333333333333333 == dataframe_RSCU_CAI['AAA']['lcl|NC_000913.3_cds_NP_414542.1_1'])
		self.assertEqual(0.8295403880420882, dataframe_CAI['lcl|NC_000913.3_cds_NP_414542.1_1'])
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

		## remove temp files
		#self.utils.remove_file(csv_result_RSCU_CAI)
		#self.utils.remove_file(csv_result_counts)
		
if __name__ == "__main__":
	# import sys;sys.argv = ['', 'Test.test_read_fasta']
	unittest.main()
