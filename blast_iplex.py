#!/usr/bin/env python
#
# Copyright (c) 2018 Ernest K Lee.  All rights reserved.
#
# BLAST iPLEX candidates against nr.
# Input is output from search_iplex.py.
#

import re
import argparse
import os
import tempfile
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Number of threads for BLAST
nthreads = 16
# E-value cufoff for BLAST
eval_cutoff = 1e-10

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("iplexfile", help="iPLEX candidates output from search_iplex.py")
    args = parser.parse_args()
    iplexfile = args.iplexfile

    # Write all iPLEX sequences into a temporary query file.
    queryfile = tempfile.NamedTemporaryFile()
    with open(iplexfile, 'r') as infile:
        next(infile)  # read header
        for line in infile:
            m = re.match('(.*)\t([ACGTN]*)\[([ACGT])/[^\]]*\]([ACGT]*)\t[ACGTN]+$', line)
            seq_name = m.group(1)
            seq_str = m.group(2) + m.group(3) + m.group(4)
            queryfile.write('>{}\n{}\n'.format(seq_name, seq_str))
    queryfile.flush()

    # Run megablast against nt.
    outtempfile = tempfile.NamedTemporaryFile()
    blastn_cline = NcbiblastnCommandline(query=queryfile.name, out=outtempfile.name, db="nt", evalue=eval_cutoff, outfmt=5, num_threads=nthreads, max_target_seqs=1, task='megablast')
    stdout, stderr = blastn_cline()

    # Write results to stdout, including queries with no hits
    print 'Sequence Name\tBLAST Hit\tE-Value'
    if os.stat(outtempfile.name).st_size > 0:
        blast_records = list(NCBIXML.parse(outtempfile))
        for blast_record in blast_records:
            if (len(blast_record.alignments) > 0):
                alignment = blast_record.alignments[0]
                hsp = alignment.hsps[0]
                print '{}\t{}\t{}'.format(blast_record.query, alignment.title, hsp.expect)
            else:
                print blast_record.query
    outtempfile.close()
    queryfile.close()
    return(0)

if __name__ == '__main__':
    main()

