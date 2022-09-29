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

        ## call method
        dataframe_counts, dataframe_RSCU_CAI, dataframe_CAI = read_genome(ecoli_fasta)

        ## not finished, test if this gene is inside data frame
        self.assertTrue(dataframe_counts['lcl|NC_000913.3_cds_NP_414542.1_1'] in dataframe_counts)
        self.assertEqual(1, dataframe_counts['lcl|NC_000913.3_cds_NP_414542.1_1']['AAA'])

    def test_tables(self):
        ecoli_fasta = os.path.join(self.baseDirectory, "files/references/ecoli.fasta")
        expected_result = os.path.join(self.baseDirectory, "files/tables/result_temp.csv")
        self.assertTrue(os.path.exists(expected_result))
        dataframe_counts, dataframe_RSCU_CAI, dataframe_CAI = read_genome(ecoli_fasta)
        csv_result = self.utils.get_temp_file("abc", ".csv")
        save_table(dataframe_counts, csv_result)
        self.assertTrue(filecmp.cmp(expected_result, csv_result))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_read_fasta']
    unittest.main()
