# WisecondorX

import logging

import numpy as np
import pysam

'''
Converts bam file to numpy array by transforming
individual reads to counts per bin.
'''


def convert_bam(args):
    bins_per_chr = dict()
    for chr in range(1, 25):
        bins_per_chr[str(chr)] = None

    logging.info('Importing data ...')

    bam_file = pysam.AlignmentFile(args.infile, 'rb')

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0
    larp = -1
    larp2 = -1

    logging.info('Converting bam ... This might take a while ...')

    for index, chr in enumerate(bam_file.references):

        chr_name = chr
        if chr_name[:3].lower() == 'chr':
            chr_name = chr_name[3:]
        if chr_name not in bins_per_chr and chr_name != 'X' and chr_name != 'Y':
            continue

        logging.info('Working at {}; processing {} bins'
                     .format(chr, int(bam_file.lengths[index] / float(args.binsize) + 1)))
        counts = np.zeros(int(bam_file.lengths[index] / float(args.binsize) + 1), dtype=np.int32)
        bam_chr = bam_file.fetch(chr)

        if chr_name == 'X':
            chr_name = '23'
        if chr_name == 'Y':
            chr_name = '24'

        for read in bam_chr:
            if read.is_paired:
                if not read.is_proper_pair:
                    reads_pairf += 1
                    continue
                if not args.normdup and larp == read.pos and larp2 == read.next_reference_start:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / args.binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                larp2 = read.next_reference_start
                reads_seen += 1
                larp = read.pos
            else:
                if not args.normdup and larp == read.pos:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / args.binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                reads_seen += 1
                larp = read.pos

        bins_per_chr[chr_name] = counts
        reads_kept += sum(counts)

    qual_info = {'mapped': bam_file.mapped,
                 'unmapped': bam_file.unmapped,
                 'no_coordinate': bam_file.nocoordinate,
                 'filter_rmdup': reads_rmdup,
                 'filter_mapq': reads_mapq,
                 'pre_retro': reads_seen,
                 'post_retro': reads_kept,
                 'pair_fail': reads_pairf}
    return bins_per_chr, qual_info
