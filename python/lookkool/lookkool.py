'''

lookkool

Uses a shared library to find palindromes

created on March 7, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

from ctypes import c_char_p, c_int, cdll, pointer
import os 

libpath = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'liblookkool.so')
lookkool = cdll.LoadLibrary(libpath)


def findPalindromes(contig):
    '''
    Uses a c library to find palindromes and return the locations
    '''

    c_contig_len = c_int(len(contig))
    max_locs = 2 * len(contig)
    c_locs_array_len = c_int(max_locs)
    c_locs_array = c_int * (2 * len(contig))
    c_locs_array_ptr = pointer(c_locs_array())


    lookkool.find_palindromes(c_char_p(contig),c_contig_len,c_locs_array_ptr,c_locs_array_len)

    locations = []
    i = 0
    while i < max_locs - 2:
        if c_locs_array_ptr.contents[i] != 0:
            start = c_locs_array_ptr.contents[i]
            end = c_locs_array_ptr.contents[i+1]
            locations.append((start,end))
        i += 2
    return locations
