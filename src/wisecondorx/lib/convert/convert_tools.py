import logging
import re
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pysam


def convert_reads(
    binsize: int, infile: Path, reference: Path, normdup: bool
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Converts aligned reads file to numpy array by transforming
    individual reads to counts per bin.
    """
    bins_per_chr = dict()
    for chrom in range(1, 25):
        bins_per_chr[str(chrom)] = None

    logging.info("Importing data ...")

    reads_file: pysam.AlignmentFile
    if infile.suffix == ".bam":
        reads_file = pysam.AlignmentFile(infile, "rb")
    elif infile.suffix == ".cram":
        if reference:
            if not reference.exists():
                logging.error("Invalid reference file")
                sys.exit(1)
            reads_file = pysam.AlignmentFile(
                infile, "rc", reference_filename=reference
            )
        else:
            logging.error(
                "Cram support requires a reference file, please use the --reference argument"
            )
            sys.exit(1)
    else:
        logging.error(
            "Unsupported input file type. Make sure your input filename has a correct extension ( bam or cram)"
        )
        sys.exit(1)

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0
    larp = -1
    larp2 = -1

    logging.info("Converting aligned reads ... This might take a while ...")

    chr_regex = re.compile(r"^chr", re.IGNORECASE)
    for index, chrom in enumerate(reads_file.references):
        chr_name = chr_regex.sub("", chrom)
        if (
            chr_name not in bins_per_chr
            and chr_name != "X"
            and chr_name != "Y"
        ):
            continue

        logging.info(
            f"Working at {chrom}; processing {int(reads_file.lengths[index] / float(binsize) + 1)} bins"
        )
        counts = np.zeros(
            int(reads_file.lengths[index] / float(binsize) + 1), dtype=np.int32
        )
        bam_chr = reads_file.fetch(chrom)

        if chr_name == "X":
            chr_name = "23"
        if chr_name == "Y":
            chr_name = "24"

        for read in bam_chr:
            if read.is_paired:
                if not read.is_proper_pair:
                    reads_pairf += 1
                    continue
                if (
                    not normdup
                    and larp == read.pos
                    and larp2 == read.next_reference_start
                ):
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                larp2 = read.next_reference_start
                reads_seen += 1
                larp = read.pos
            else:
                if not normdup and larp == read.pos:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                reads_seen += 1
                larp = read.pos

        bins_per_chr[chr_name] = counts
        reads_kept += sum(counts)

    qual_info = {
        "mapped": reads_file.mapped,
        "unmapped": reads_file.unmapped,
        "no_coordinate": reads_file.nocoordinate,
        "filter_rmdup": reads_rmdup,
        "filter_mapq": reads_mapq,
        "pre_retro": reads_seen,
        "post_retro": reads_kept,
        "pair_fail": reads_pairf,
    }
    logging.debug(f"Quality info: {qual_info}")

    return bins_per_chr, qual_info
