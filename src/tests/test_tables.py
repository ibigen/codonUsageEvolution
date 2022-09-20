'''
Created on 17/09/2022

@author: mmp
'''
import unittest
import os
from codon_usage_evolution import read_genome

genome = read_genome()

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

        vect_genome = read_genome(ecoli_fasta)
        self.assertEqual(12, len(vect_genome))
        
        
    def test_tables(self):
        expected_result = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/tables/result_temp.csv")
        self.assertTrue(os.path.exists(expected_result))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_read_fasta']
    unittest.main()