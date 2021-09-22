#!/usr/bin/env python

from typing import List

from wisecondorx.lib.predict.predict_tools import (
    coverage_normalize_and_mask,
    get_optimal_cutoff,
    get_weights,
    inflate_results,
    normalize_repeat,
    project_pc,
)


def normalize(sample, ref_file, ref_gender, maskrepeats):
    """
    Control function that executes following normalization strategies:
    - coverage normalization
    - between-sample normalization
    - within-sample normalization
    """
    if ref_gender == "A":
        ap = ""
        cp = 0
        ct = 0
    else:
        ap = ".{}".format(ref_gender)
        cp = 22
        ct = ref_file["masked_bins_per_chr_cum{}".format(ap)][cp - 1]

    sample = coverage_normalize_and_mask(sample, ref_file, ap)
    sample = project_pc(sample, ref_file, ap)
    results_w = get_weights(ref_file, ap)[ct:]
    optimal_cutoff = get_optimal_cutoff(ref_file, maskrepeats)
    results_z, results_r, ref_sizes, m_lr, m_z = normalize_repeat(
        sample, ref_file, optimal_cutoff, ct, cp, ap
    )

    return results_r, results_z, results_w, ref_sizes, m_lr, m_z


def get_post_processed_result(
    minrefbins, result, ref_sizes, mask, bins_per_chr: List[int]
):
    """
    Function processes a result (e.g. results_r) to an easy-to-interpret format. Bins without information are set to 0.
    """
    infinite_mask = ref_sizes < minrefbins
    result[infinite_mask] = 0
    inflated_results = inflate_results(result, mask)

    final_results = []
    for chrom in range(len(bins_per_chr)):
        chr_data = inflated_results[
            sum(bins_per_chr[:chrom]) : sum(bins_per_chr[: chrom + 1])
        ]
        final_results.append(chr_data)

    return final_results
