#!/usr/local/bin/python

'''
hisnhers.py
Harvard Informatics Script for Nextgen HiSeq Extraction and Reporting of Sequences

'''

import os, traceback, re
import json
import subprocess
import time


def fastqToSequenceList(fileh):
    '''
    Takes a fastq file handle, returns a list tuples including
    seqid, sequence bases, and quality scores
    '''
      seqs = []
      seqid = bases = qscores = None

    if fileh.closed:
        raise Exception('fastq file is closed.')
        
    for line in fileh:
        line = line.strip()
        if line == '':
            continue
        if line.startswith('@'):
            # It's a new record
            if seqid is not None:
                seqs.append((seqid,bases,qscores))
            seqid = bases = qscores = None

            m = re.match(r'^@([^ ]+).*',line)
            if m is None:
                raise Exception('No sequence identifier found after the @ symbol')
            seqid = m.group(1).strip()

        elif seqid is not None and bases is None and line.strip() != '+':
            bases = line.strip()
        elif seqid is not None and bases is not None and line.strip() != '+':
            qscores = line.strip()

    return seqs


def main():

    # Read fastq file and report length, base counts
    seqs = []
    fqfilename = '/n/regal/informatics/aaron/testfile.fq'
    seqs = fastqToSequenceList(fqfilename)

    # Process the data that is in the seqs list
    # seqs[0]
    # ('HWUSI-EAS300R_0005_F2AAXX:8:30:18447:12115#0/1\n', 'CGTAGCTAACGTATTCACCGTG', '')
    # for i,seqdata in enumerate(seqs):
    #    Do some stuff here
    #    print basecountline

    # Write out sequences in fasta format
    (path,ext) = os.path.splitext(fqfilename)
    fafilename = path + '.fa'
    print('Writing to %s' % fafilename)
    with open(fafilename,'w') as f:
        for seqdata in seqs:
            f.write('>%s\n%s\n' % (seqdata[0],seqdata[1]))


    # Run megaAssembler with fasta file input and the specified contig output filename
    contigfilename = '%s.contigs' % fafilename
    
    cmd = 'megaAssembler {inputfile} {outputfile}'.format(inputfile=fafilename,outputfile=contigfilename)
    print('Running %s' % cmd)
    os.system('%s > /dev/null 2> /dev/null'  % cmd)

    # Create an array of contig data 
    contigs = []
    with open(contigfilename,'r') as c:
        seqid = None
        for line in c:
            line = line.strip()
            if line == '':
                continue
            # Check for seq id line
            m = re.match(r'^>\s*([^ ]+).*', line)
            if m is not None:
                seqid = m.group(1)
            else:
                # Otherwise, it's the bases
                contigs.append((seqid, line))
                seqid = None

    # Annotate the contigs, store in a dictionary and then write out JSON
    from ha.annotate import annotateStartStopCodons, annotatePalindromes

    # One at a time
    annotations = []
    starttime = time.time()
    for seqid, contig in contigs:
        annotations += annotateStartStopCodons(seqid, contig)
        annotations += annotatePalindromes(seqid, contig)
    endtime = time.time()
    print('Elapsed serial annotation time %d seconds' % int(endtime - starttime))

    # Make a dictionary keyed by contig name
    annotatedcontigs = {}
    for annotation in annotations:
        annotatedcontigs.setdefault(annotation['seqid'],[]).append(annotation)

    # Dump annotations in JSON form
    with open('%s.annotations' % fafilename, 'w') as f:
        f.write(json.dumps(annotatedcontigs,indent=4))

if __name__ == "__main__":
    sys.exit(main())

