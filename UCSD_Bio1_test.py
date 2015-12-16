import unittest
from UCSD_Bio1 import *


class SimplisticTest(unittest.TestCase):

    def test(self):
        self.assertTrue(True)

class testRos(unittest.TestCase):

    def test_ros_mprt(self):
        inp = ['A2Z669', 'B5ZC00', 'P07204_TRBM_HUMAN', 'P20840_SAG1_YEAST']
        out = ['B5ZC00', '85 118 142 306 395', 'P07204_TRBM_HUMAN',
               '47 115 116 382 409', 'P20840_SAG1_YEAST',
               '79 109 135 248 306 348 364 402 485 501 614',]
        for (out_el, expected_el) in zip(ros_mprt(imp),out):
            self.assertTrue(out_el==expected_el)

    def test_ros_grph(self):
        adjList = ros_grph('dataset_ros_grph_test.txt')
        self.assertEqual(len(adjList),3)
##Rosalind_0498 Rosalind_2391
##Rosalind_0498 Rosalind_0442
##Rosalind_2391 Rosalind_2323

if __name__ == '__main__':
    unittest.main()

