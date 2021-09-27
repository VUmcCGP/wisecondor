#!/usr/bin/env python

import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from wisecondorx.lib.utils import (
    exec_R,
    get_cpa,
    get_median_segment_variance,
    get_z_score,
)


def exec_write_plots(
    outid: str,
    ref_gender: str,
    beta: float,
    zscore: float,
    binsize: int,
    n_reads: int,
    cairo: bool,
    ylim: Tuple[int, int],
    results: Dict[str, Any],
):
    """
    Writes plots.
    """
    json_plot_dir = os.path.abspath(outid + "_plot_tmp")
    json_dict = {
        "R_script": Path(os.path.dirname(os.path.realpath(__file__)))
        / "include"
        / "plotter.R",
        "ref_gender": ref_gender,
        "beta": beta,
        "zscore": zscore,
        "binsize": binsize,
        "n_reads": n_reads,
        "cairo": cairo,
        "results_r": results["results_r"],
        "results_w": results["results_w"],
        "results_c": results["results_c"],
        "ylim": ylim,
        "infile": f"{json_plot_dir}.json",
        "out_dir": "{outid}.plots",
    }
    exec_R(json_dict)


def generate_output_tables(
    outid: str,
    binsize: int,
    ref_gender: str,
    beta: float,
    zscore: float,
    bins_per_chr: List[int],
    gender: str,
    n_reads: int,
    results: Dict[str, Any],
) -> None:
    """
    Calculates zz-scores, marks aberrations and writes tables.
    """
    _generate_bins_bed(outid=outid, binsize=binsize, results=results)
    _generate_segments_and_aberrations_bed(
        outid=outid,
        binsize=binsize,
        ref_gender=ref_gender,
        beta=beta,
        zscore=zscore,
        results=results,
    )
    _generate_chr_statistics_file(
        outid=outid,
        binsize=binsize,
        bins_per_chr=bins_per_chr,
        gender=gender,
        n_reads=n_reads,
        results=results,
    )


def _generate_bins_bed(
    outid: str, binsize: int, results: Dict[str, Any]
) -> None:
    bins_file = open(f"{outid}_bins.bed", "w")
    bins_file.write("chr\tstart\tend\tid\tratio\tzscore\n")
    results_r = results["results_r"]
    results_z = results["results_z"]
    binsize = binsize

    for chr in range(len(results_r)):
        chr_name = str(chr + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"
        feat = 1
        for i in range(len(results_r[chr])):
            r = results_r[chr][i]
            z = results_z[chr][i]
            if r == 0:
                r = "nan"
            if z == 0:
                z = "nan"
            feat_str = "{}:{}-{}".format(
                chr_name, str(feat), str(feat + binsize - 1)
            )
            row = "\t".join(
                [chr_name, feat, feat + binsize - 1, feat_str, r, z]
            )
            bins_file.write(f"{row}\n")
            feat += binsize
    bins_file.close()


def _generate_segments_and_aberrations_bed(
    outid: str,
    binsize: int,
    ref_gender: str,
    beta: float,
    zscore: int,
    results: Dict[str, Any],
) -> None:
    segments_file = open(f"{outid}_segments.bed", "w")
    abberations_file = open(f"{outid}_aberrations.bed", "w")
    segments_file.write("chr\tstart\tend\tratio\tzscore\n")
    abberations_file.write("chr\tstart\tend\tratio\tzscore\ttype\n")

    for segment in results["results_c"]:
        chr_name = str(segment[0] + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"
        row = [
            chr_name,
            int(segment[1] * binsize + 1),
            int(segment[2] * binsize),
            segment[4],
            segment[3],
        ]
        segments_file.write("{}\n".format("\t".join([str(x) for x in row])))

        ploidy = 2
        if (chr_name == "X" or chr_name == "Y") and ref_gender == "M":
            ploidy = 1
        if beta:
            if float(segment[4]) > __get_aberration_cutoff(beta, ploidy)[1]:
                abberations_file.write(
                    "{}\tgain\n".format("\t".join([str(x) for x in row]))
                )
            elif float(segment[4]) < __get_aberration_cutoff(beta, ploidy)[0]:
                abberations_file.write(
                    "{}\tloss\n".format("\t".join([str(x) for x in row]))
                )
        elif isinstance(segment[3], str):
            continue
        else:
            if float(segment[3]) > zscore:
                abberations_file.write(
                    "{}\tgain\n".format("\t".join([str(x) for x in row]))
                )
            elif float(segment[3]) < -zscore:
                abberations_file.write(
                    "{}\tloss\n".format("\t".join([str(x) for x in row]))
                )

    segments_file.close()
    abberations_file.close()


def __get_aberration_cutoff(beta: float, ploidy: int) -> Tuple[float, float]:
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def _generate_chr_statistics_file(
    outid: str,
    binsize: int,
    bins_per_chr: List[int],
    gender: str,
    n_reads: int,
    results: Dict[str, Any],
):
    stats_file = open(f"{outid}_statistics.txt", "w")
    stats_file.write("chr\tratio.mean\tratio.median\tzscore\n")
    chr_ratio_means = [
        np.ma.average(
            results["results_r"][chr], weights=results["results_w"][chr]
        )
        for chr in range(len(results["results_r"]))
    ]
    chr_ratio_medians = [
        np.median([x for x in results["results_r"][chr] if x != 0])
        for chr in range(len(results["results_r"]))
    ]

    results_c_chr = [
        [x, 0, bins_per_chr[x] - 1, chr_ratio_means[x]]
        for x in range(len(results["results_r"]))
    ]

    msv = round(
        get_median_segment_variance(
            results["results_c"], results["results_r"]
        ),
        5,
    )
    cpa = round(get_cpa(results["results_c"], binsize), 5)
    chr_z_scores = get_z_score(results_c_chr, results)

    for chrom in range(len(results["results_r"])):

        chr_name = str(chrom + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"

        row = [
            chr_name,
            chr_ratio_means[chrom],
            chr_ratio_medians[chrom],
            chr_z_scores[chrom],
        ]

        stats_file.write("\t".join([str(x) for x in row]) + "\n")

    stats_file.write(
        f"Gender based on --yfrac (or manually overridden by --gender): {gender}\n"
    )

    stats_file.write(f"Number of reads: {n_reads}\n")

    stats_file.write(
        f"Standard deviation of the ratios per chromosome: {round(np.nanstd(chr_ratio_means), 5)}\n"
    )

    stats_file.write(
        f"Median segment variance per bin (doi: 10.1093/nar/gky1263): {msv}\n"
    )

    stats_file.write(
        f"Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4): {cpa}\n"
    )

    stats_file.close()
