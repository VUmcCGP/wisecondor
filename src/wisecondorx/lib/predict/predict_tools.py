#!/usr/bin/env python

import os
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
from wisecondorx.lib.utils import exec_R, get_z_score
from scipy.stats import norm
from sklearn.decomposition import PCA


def predict_gender(sample: np.ndarray.item, trained_cutoff: float) -> str:
    """
    Returns gender based on Gaussian mixture
    model trained during newref phase.
    """
    y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    if y_fraction > trained_cutoff:
        return "M"
    else:
        return "F"


def coverage_normalize_and_mask(sample: np.ndarray.item, ref_file, ap):
    """
    Normalize sample for read depth and apply mask.
    """
    by_chr = []

    chrs = range(1, len(ref_file["bins_per_chr{}".format(ap)]) + 1)

    for chr in chrs:
        this_chr = np.zeros(
            ref_file["bins_per_chr{}".format(ap)][chr - 1], dtype=float
        )
        min_len = min(
            ref_file["bins_per_chr{}".format(ap)][chr - 1],
            len(sample[str(chr)]),
        )
        this_chr[:min_len] = sample[str(chr)][:min_len]
        by_chr.append(this_chr)
    all_data = np.concatenate(by_chr, axis=0)
    all_data = all_data / np.sum(all_data)
    masked_data = all_data[ref_file["mask{}".format(ap)]]

    return masked_data


def project_pc(sample_data, ref_file, ap):
    """
    Project test sample to PCA space.
    """
    pca = PCA(n_components=ref_file["pca_components{}".format(ap)].shape[0])
    pca.components_ = ref_file["pca_components{}".format(ap)]
    pca.mean_ = ref_file["pca_mean{}".format(ap)]

    transform = pca.transform(np.array([sample_data]))

    reconstructed = np.dot(transform, pca.components_) + pca.mean_
    reconstructed = reconstructed[0]
    return sample_data / reconstructed


def get_optimal_cutoff(ref_file, repeats):
    """
    Defines cutoff that will add bins to a blacklist depending on the within reference distances.
    """
    distances = ref_file["distances"]
    cutoff = float("inf")
    for i in range(0, repeats):
        mask = distances < cutoff
        average = np.average(distances[mask])
        stddev = np.std(distances[mask])
        cutoff = average + 3 * stddev
    return cutoff


def normalize_repeat(test_data, ref_file, optimal_cutoff, ct, cp, ap):
    """
    Within sample normalization.
    Cycles through a number of repeats where z-scores define whether a bin is seen as
    normal' in a sample (in most cases this means non-aberrant').
    If it is, it can be used as a reference for the other bins.
    """
    results_z = None
    results_r = None
    ref_sizes = None
    test_copy = np.copy(test_data)
    for i in range(3):
        results_z, results_r, ref_sizes = _normalize_once(
            test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
        )

        test_copy[ct:][np.abs(results_z) >= norm.ppf(0.99)] = -1
    m_lr = np.nanmedian(np.log2(results_r))
    m_z = np.nanmedian(results_z)

    return results_z, results_r, ref_sizes, m_lr, m_z


def _normalize_once(
    test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
):
    masked_bins_per_chr = ref_file["masked_bins_per_chr{}".format(ap)]
    masked_bins_per_chr_cum = ref_file["masked_bins_per_chr_cum{}".format(ap)]
    results_z = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    results_r = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    ref_sizes = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    indexes = ref_file["indexes{}".format(ap)]
    distances = ref_file["distances{}".format(ap)]

    i = ct
    i2 = 0
    for chr in list(range(len(masked_bins_per_chr)))[cp:]:
        start = masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
        end = masked_bins_per_chr_cum[chr]
        chr_data = np.concatenate(
            (
                test_copy[
                    : masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
                ],
                test_copy[masked_bins_per_chr_cum[chr] :],
            )
        )

        for index in indexes[start:end]:
            ref_data = chr_data[index[distances[i] < optimal_cutoff]]
            ref_data = ref_data[ref_data >= 0]
            ref_stdev = np.std(ref_data)
            results_z[i2] = (test_data[i] - np.mean(ref_data)) / ref_stdev
            results_r[i2] = test_data[i] / np.median(ref_data)
            ref_sizes[i2] = ref_data.shape[0]
            i += 1
            i2 += 1

    return results_z, results_r, ref_sizes


def get_weights(ref_file, ap):
    """
    The means of sets of within-sample reference distances can serve as inverse weights for CBS, Z-scoring and plotting.
    """
    inverse_weights = [
        np.mean(np.sqrt(x)) for x in ref_file["distances{}".format(ap)]
    ]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


def inflate_results(results, mask):
    """
    Unmasks results array.
    """
    temp = [0 for x in mask]
    j = 0
    for i, val in enumerate(mask):
        if val:
            temp[i] = results[j]
            j += 1
    return temp


def log_trans(results, log_r_median):
    """
    Log2-transforms results_r. If resulting elements are infinite, all corresponding possible positions
    (at results_r, results_z and results_w are set to 0 (blacklist)).
    """
    for chr in range(len(results["results_r"])):
        results["results_r"][chr] = np.log2(results["results_r"][chr])

    results["results_r"] = [x.tolist() for x in results["results_r"]]

    for c in range(len(results["results_r"])):
        for i, rR in enumerate(results["results_r"][c]):
            if not np.isfinite(rR):
                results["results_r"][c][i] = 0
                results["results_z"][c][i] = 0
                results["results_w"][c][i] = 0
            if results["results_r"][c][i] != 0:
                results["results_r"][c][i] = (
                    results["results_r"][c][i] - log_r_median
                )


def apply_blacklist(
    blacklist_bed: Path, binsize: int, results: Dict[str, Any]
):
    """
    Applies additional blacklist to results_r, results_z and results_w if requested.
    """
    blacklist = _import_bed(blacklist_bed=blacklist_bed, binsize=binsize)

    for chr in blacklist.keys():
        for s_e in blacklist[chr]:
            for pos in range(s_e[0], s_e[1]):
                if len(results["results_r"]) < 24 and chr == 23:
                    continue
                if pos >= len(results["results_r"][chr]) or pos < 0:
                    continue
                results["results_r"][chr][pos] = 0
                results["results_z"][chr][pos] = 0
                results["results_w"][chr][pos] = 0


def _import_bed(blacklist_bed, binsize):
    bed = {}
    for line in open(blacklist_bed):
        chr_name, s, e = line.strip().split("\t")
        if chr_name[:3] == "chr":
            chr_name = chr_name[3:]
        if chr_name == "X":
            chr_name = "23"
        if chr_name == "Y":
            chr_name = "24"
        chr = int(chr_name) - 1
        if chr not in bed.keys():
            bed[chr] = []
        bed[chr].append(
            bed[chr].append([int(int(s) / binsize), int(int(e) / binsize) + 1])
        )
    return bed


def exec_cbs(
    alpha: float,
    ref_gender: str,
    binsize: int,
    outid: str,
    results: Dict[str, Any],
):
    """
    Executed CBS on results_r using results_w as weights. Calculates segmental zz-scores.
    """
    json_cbs_dir = os.path.abspath(outid + "_CBS_tmp")

    json_dict = {
        "R_script": Path(os.path.dirname(os.path.realpath(__file__)))
        / "include"
        / "CBS.R",
        "ref_gender": ref_gender,
        "alpha": alpha,
        "binsize": binsize,
        "results_r": results["results_r"],
        "results_w": results["results_w"],
        "infile": Path(f"{json_cbs_dir}_01.json"),
        "outfile": Path(f"{json_cbs_dir}_02.json"),
    }

    results_c = _get_processed_cbs(exec_R(json_dict))
    segment_z = get_z_score(results_c, results)
    results_c = [
        results_c[i][:3] + [segment_z[i]] + [results_c[i][3]]
        for i in range(len(results_c))
    ]
    return results_c


def _get_processed_cbs(cbs_data) -> Tuple[int, int, int, float]:
    results_c = []
    for i, segment in enumerate(cbs_data):
        chr = int(segment["chr"]) - 1
        s = int(segment["s"])
        e = int(segment["e"])
        r = segment["r"]
        results_c.append((chr, s, e, r))

    return results_c
