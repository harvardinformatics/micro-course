'''
Created on March 7, 2017
Copyright (c) 2017
Harvard FAS Research Computing
All rights reserved.

@author: Aaron Kitzmiller
'''
from lookkool import findPalindromes

import unittest

class Test(unittest.TestCase):
    '''
    Test the findPalindromes method
    '''
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testPalindromes(self):

        contig = 'ATCGCGAT'
        locations = findPalindromes(contig)
        self.assertTrue(len(locations) == 1, 'Incorrect number of locations %d' % len(locations))
        self.assertTrue(locations[0][0] == 1, 'Incorrect start %d' % locations[0][0])
        self.assertTrue(locations[0][1] == 8, 'Incorrect end %d' % locations[0][1])
