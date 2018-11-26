#!/usr/bin/env python
#
# Copyright (c) 2018 Ernest K Lee.  All rights reserved.
#
# Search for iplex compatible segments from Mauve genome alignment.
# Samples in alignment, but not listed in config file are ignored.
#
# Input:
# - json config file specifying the reference genomes/scaffolds, Mauve alignment file, and replicate BAMs.
#
# Output:
# - iPLEX candidates are written to stdout, log messages to stderr. 
#

import sys
import argparse
import re
import json
from StringIO import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
import pysam
from pyfaidx import Fasta

# Customize the following settings as needed
#
# Minimum number of replicates with the same variant
MIN_REPS = 3
# Maximum number of snps allowed within iplex segment
MAX_SNPS = 10
# Minumum number of bases a snp must be from another snp for it to be chosen
MIN_DIST = 10

# Mauve coordinates of genome/scaffold.
scaf_coord = {}
# Genomes/Scaffolds
genome = {}
# Complement
dna_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def read_faidx(sample, fasta):
    """Read genome reference fasta index"""
    scaffolds = Fasta(fasta, as_raw=True)
    genome[sample] = scaffolds
    return True

def get_iplex(consensus_seq=''):
    """Return list of potential iplex compatible segments given an annotated consensus sequence from an alignment."""
    # Pattern potentially suitable for iPLEX design.
    iplex_pattern = '[ACGT]{40,70}_[ACGT]{40,70}|[ACGT]{30,40}_[ACGT]{50,100}|[ACGT]{50,100}_[ACGT]{30,40}|[ACGT]{40,70}[ACGT_]{1,40}_[ACGT]{40,70}';
    matches = []
    p = re.compile(iplex_pattern)
    for m in p.finditer(consensus_seq):
        matches.append(m)
    return matches

def filter_iplex(seq_str):
    """Additional heuristics to filter iplex candidates.  Return string with all but one snp
    masked, otherwise return an empty string if sequence is not suitable for iplex.
    """
    snp_count = seq_str.count('_')
    last_snp_pos = seq_str.rindex('_')
    if snp_count == 1:
        return seq_str
    else:
        penultimate_snp_pos = seq_str.rindex('_', 0, last_snp_pos)
        if snp_count > MAX_SNPS or (last_snp_pos-penultimate_snp_pos) < MIN_DIST:
            sys.stderr.write('Additional filtering not passed: snp_count={} snp_dist={}\n'.format(snp_count, last_snp_pos-penultimate_snp_pos))
            return ''
        else:
            return seq_str.replace('_', 'N', snp_count-1)

def scaf_to_mauve_coord(ref): 
    """Map genome/scaffold to Mauve coordinates."""
    for sample, fasta in ref.iteritems():
        scaf_coord[sample] = []
        tot_len = 1  # 1-based
        for record in SeqIO.parse(fasta, "fasta"):
            scaf_coord[sample].append((record.id, tot_len, tot_len + len(record) - 1))
            tot_len += len(record)
    return True

def get_scaffold_coord(sample, pos):
    """Return chromosome/scaffold id and position given sample name and Mauve position."""
    # should do this more efficiently
    x = scaf_coord[sample]
    for n in range(1, len(x)-1):
        if pos >= int(x[n-1][1]) and pos < int(x[n][1]):
            return x[n-1][0], pos - int(x[n-1][1]) + 1
    return x[-1][0], pos - int(x[-1][1]) + 1

def get_direction(sample, pos, iseq):
    """Check if iplex sequence in alignment is forward or reverse complement of the scaffold
    Return values: 1=forward -1=reverse-complement 0=unknown
    """
    iseq_nogap = iseq.ungap("-")
    scaf_id = get_scaffold_coord(sample, pos)[0]
    scaf_seq = genome[sample][scaf_id][:]
    if scaf_seq.find(str(iseq_nogap)) > -1:
        return 1
    elif scaf_seq.find(str(iseq_nogap.reverse_complement())) > -1:
        return -1
    else:
        return 0

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Config file specifying reference, BAMs, and Mauve alignment")
    args = parser.parse_args()
    config_file = args.config

    with open(config_file) as json_file:
        config = json.load(json_file)
        reference = config['reference']
        bam_files = config['bams']
        mauve_alignment = config['mauve_alignment']

    sys.stderr.write("Configuration:\n")
    sys.stderr.write('Reference: {}\nBAMs: {}\nMauve alignment: {}\n'.format(reference, bam_files, mauve_alignment))

    samples = sorted(reference.keys())
    reference_files = set(reference.values())

    # Read genome index
    for sample, fasta in reference.iteritems():
        read_faidx(sample, fasta)

    # Convert genome to Mauve coordinates
    scaf_to_mauve_coord(reference)

    # Output header
    print 'Alignment:Interval\tSequence\tSNP {}'.format(samples)

    # Process Mauve alignments
    n = 0
    for aln in AlignIO.parse(open(mauve_alignment), "mauve"):

        # Only looks into alignments with coverage across all samples of consideration.
        aln_sample_files = set()
        for record in aln:
            m = re.match('(.*)/.*', record.id)
            if m:
                aln_sample_files.add(m.group(1))
        if reference_files <= aln_sample_files:
        #if len(aln) >= len(samples)
            new_aln = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
            if reference_files < aln_sample_files:
                # Create new alignment with only relevant samples in it
                for record in aln:
                    m = re.match('(.*)/.*', record.id)
                    if m and m.group(1) in reference_files:
                        new_aln.add_sequence(record.id, str(record.seq))
            else:
                new_aln = aln
            summary_align = AlignInfo.SummaryInfo(new_aln)
            consensus = summary_align.gap_consensus(1.0, '_', require_multiple=1)
            #sys.stderr.write('Consensus: {}\n'.format(str(consensus)))
            iplex_match = get_iplex(str(consensus))
            if len(iplex_match) > 0:  # there is at least 1 match!
                #SeqIO.write(aln, sys.stdout, "clustal")
                for s in iplex_match:

                    iplex_str = filter_iplex(s.group(0))
                    if iplex_str == '':
                        continue

                    passed = True  # sequence is good to use for iplex
                    snp = ''

                    # SNP position of iplex fragment in alignment.
                    # Note the snp and alignment object positions are 0-based,
                    # but the genome and mauve coords are 1-based.
                    snp_pos = iplex_str.index('_')

                    # Value at snp position in sorted order of sample names.
                    snpval = {}
                    for record in new_aln:
                        m = re.match('(.*)/.*', record.id)
                        snpval[m.group(1)] = record.seq[s.start() + snp_pos] 
                    for spl in samples:
                        snp += snpval[reference[spl]]

                    iplex_sequence = iplex_str[:snp_pos] + "[" + '/'.join(set(snp)) + "]" + iplex_str[snp_pos+1:]
                    sys.stderr.write('>>> Alignment {}\n'.format(n))
                    sys.stderr.write('{}\n'.format(str(aln)))
                    sys.stderr.write('iPLEX sequence: [{}-{}] {}\n'.format(s.start(), s.end(), iplex_sequence))
                    sys.stderr.write('SNP: {}\n'.format(snp))
                    if snp.find('-') == -1:  # skip if snp is a gap

                        # Make sure we have the same base across all reads at the snp position.
                        for record in new_aln:

                            snpbase = record.seq[s.start() + snp_pos]

                            # Alignment interval in Mauve coordinates.
                            m = re.match('(.*)/(\d+)-(\d+)', record.id)
                            sample_scaf_file =  m.group(1)
                            aln_start = int(m.group(2)) + 1
                            aln_end = int(m.group(3)) + 1

                            # Find snp position in scaffold/chromosome without gaps.
                            seq_prefix = record.seq[:(s.start()+snp_pos)]
                            num_gaps = seq_prefix.count('-')
                            snp_pos_aln_nogap = s.start() + snp_pos - num_gaps

                            sample = ''
                            for sp, fa in reference.iteritems():
                                if fa == sample_scaf_file:
                                    sample = sp
                                    break
                            sys.stderr.write('{}|'.format(sample))
                            c = []  # snp position in genome
                            direction = get_direction(sample, aln_start + snp_pos_aln_nogap, record.seq[s.start():s.end()+1])
                            if direction == 1:
                                sys.stderr.write('+')
                                c = get_scaffold_coord(sample, aln_start + snp_pos_aln_nogap)
                            elif direction == -1:
                                sys.stderr.write('-')
                                c = get_scaffold_coord(sample, aln_end - 1 - snp_pos_aln_nogap)
                            else:
                                # Cannot find iplex segment in scaffold, possibly due
                                # to iplex sequence spanning two concatenated scaffolds in Mauve.
                                sys.stderr.write('Failed to find iPLEX sequence in reference\n')
                                passed = False
                                break
                            snp_scaf_coord = '{}:{}-{}'.format(c[0], c[1], c[1])
                            sys.stderr.write('{} '.format(snp_scaf_coord))

                            # Check if all reads of at least MIN_REPS replicates have the same
                            # variant, and no other replicates have a different variant.
                            num_rep_with_coverage = 0
                            for bam in bam_files[sample]:
                                samfile = pysam.AlignmentFile(bam, "rb")
                                for pileupcolumn in samfile.pileup(region=snp_scaf_coord):
                                    if pileupcolumn.pos == c[1]-1:
                                        num_rep_with_coverage += 1
                                        for pileupread in pileupcolumn.pileups:
                                            if not pileupread.is_del and not pileupread.is_refskip:
                                                readbase = pileupread.alignment.query_sequence[pileupread.query_position]
                                                sys.stderr.write(readbase)
                                                if direction == 1 and readbase != snpbase:
                                                    passed = False
                                                    break
                                                elif direction == -1 and readbase != dna_complement[snpbase]:
                                                    passed = False
                                        sys.stderr.write(':')
                                samfile.close()
                            sys.stderr.write('\n')
                            if len(bam_files[sample]) >= MIN_REPS and num_rep_with_coverage < MIN_REPS:
                                sys.stderr.write('Less than {} replicates have coverage\n'.format(MIN_REPS))
                                passed = False
                            
                    else:
                        sys.stderr.write('SNP is a deletion\n')
                        passed = False

                    # Print iPLEX compatible sequence if everything looks good
                    if passed:
                        sys.stderr.write('PASS\n')
                        print 'Alignment{}:{}-{}\t{}\t{}'.format(n, s.start(), s.end(), iplex_sequence, snp)
                    else:
                        sys.stderr.write('FAIL\n')
        n += 1

    sys.stderr.write('{} alignment records processed.\n'.format(n))
    return(0)

if __name__ =='__main__':
    main()
