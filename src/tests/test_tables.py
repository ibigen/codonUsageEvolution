'''
Created on 17/09/2022

@author: mmp
'''
import unittest
import os


from codon_usage_evolution import read_genome
import filecmp
import utils as utils

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sum(self, a, b):
        """ sum two numbers, only for example"""
        return a + b

    def test_sum(self):
        self.assertEqual(5, self.sum(3, 2))
        self.assertNotEqual(6, self.sum(3, 2))

    def test_read_fasta(self):
        ecoli_fasta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/references/ecoli.fasta")
        self.assertTrue(os.path.exists(ecoli_fasta))

        ## call method
        dataframe_counts, dataframe_RSCU_CAI, dataframe_CAI = read_genome(ecoli_fasta)

        ## not finished, test if this gene is inside data frame
        self.assertTrue(dataframe_counts['lcl|NC_000913.3_cds_NP_414542.1_1'] in dataframe_counts)

        self.assertEqual(1, dataframe_counts['lcl|NC_000913.3_cds_NP_414542.1_1']['AAA'])

    def test_tables(self):
        expected_result = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/tables/result_temp.csv")
        self.assertTrue(os.path.exists(expected_result))

    def test_get_first_sequence_fasta(self):
        """ Compare created files with the expected ones """

        fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR,
                                  "ecoli.fasta")
        fasta_file_temp = self.utils.get_temp_file("fasta_lower", FileExtensions.FILE_FASTA)
        self.software.set_first_sequence_fasta(fasta_file_temp)
        self.assertTrue(filecmp.cmp(fasta_file_temp, fasta_file))



if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_read_fasta']
    unittest.main()
