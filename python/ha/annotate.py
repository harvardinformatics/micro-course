'''

ha.annotate

Annotation module for python workshop

created on March 7, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

from time import sleep

def annotateStartStopCodons(seqid, contig):
    '''
    Returns a list of dictionaries representing annotation data.
    Each dictionary includes a start, end, seqid, and key that indicates 
    the name of the annotation and it's range

    Just does start and stop codons
    '''

    stops = [
        'TAG',
        'TAA',
        'TGA',
    ]

    annotations = []
    if contig is not None and len(contig) > 3:
        for i in xrange(0,len(contig) - 2):
            start = i
            end = i + 3
            if contig[start:end] in stops:
                annotations.append(
                    {
                        'seqid' : seqid,
                        'start' : start + 1,
                        'end'   : end,
                        'key'   : 'stop_codon',
                    }
                )
            else:
                if contig[start:end] == 'ATG':
                    annotations.append( 
                        {
                            'seqid' : seqid,
                            'start' : start + 1,
                            'end'   : end,
                            'key'   : 'start_codon',
                        }
                    )
    sleep(5)

    return annotations

def annotatePalindromes(seqid, contig):
    '''
    Returns a list of dictionaries representing annotation data.
    Each dictionary includes a start, end, seqid, and key that indicates 
    the name of the annotation and it's range

    This uses the lookkool package to find palindromes
    '''
    from lookkool import findPalindromes

    palindromes = findPalindromes(contig)
    annotations = []
    for palindrome in palindromes:
        annotations.append(
            {
                'seqid'     : seqid,
                'start'     : palindrome[0],
                'end'       : palindrome[1],
                'key'       : 'palindrome',
            }
        )
    sleep(2)
    return annotations