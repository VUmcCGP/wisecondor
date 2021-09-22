#!/usr/bin/env python

import logging
import sys
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import typer
from convert_tools import convert_reads
from newref_control import (
    tool_newref_main,
    tool_newref_merge,
    tool_newref_prep,
)
from newref_tools import get_mask, train_gender_model
from overall_tools import gender_correct, scale_sample
from predict_control import get_post_processed_result, normalize
from predict_output import exec_write_plots, generate_output_tables
from predict_tools import apply_blacklist, exec_cbs, log_trans, predict_gender


class Gender(str, Enum):
    F = "F"
    M = "M"


class LogLevel(str, Enum):
    info = "info"
    warning = "warning"
    debug = "debug"
    error = "error"
    critical = "critical"


app = typer.Typer(
    help="WisecondorX — An evolved WISECONDOR",
)


@app.command("convert")
def convert(
    infile: Path = typer.Argument(
        ..., help="aligned reads input for conversion"
    ),
    outfile: Path = typer.Argument(..., help="Output .npz file"),
    reference: Optional[Path] = typer.Option(
        None, help="Fasta reference to be used during cram conversion"
    ),
    binsize: int = typer.Option(5000, help="Bin size (bp)"),
    normdup: bool = typer.Option(
        False, "--normmdup", help="Avoid remove duplicates"
    ),
) -> None:
    """
    Convert and filter a aligned reads to .npz
    """
    typer.echo("Starting conversion")

    sample, qual_info = convert_reads(
        infile=infile, binsize=binsize, reference=reference, normdup=normdup
    )
    np.savez_compressed(
        outfile, binsize=binsize, sample=sample, quality=qual_info
    )

    typer.echo("Finished conversion")


@app.command("newref")
def newref(
    infiles: List[Path] = typer.Argument(
        ...,
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    ),
    outfile: Path = typer.Argument(
        ...,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
    ),
    nipt: bool = typer.Option(False, "--nipt", help="Use flag for NIPT"),
    yfrac: float = typer.Option(
        ...,
        help="Use to manually set the y read fraction cutoff, which defines gender",
    ),
    yfracplot: Optional[Path] = typer.Option(
        None,
        help="Path to yfrac .png plot for --yfrac optimization (e.g. path/to/plot.png); software will stop after "
        "plotting after which --yfrac can be set manually",
    ),
    refsize: int = typer.Option(
        300, help="Amount of reference locations per target"
    ),
    binsize: int = typer.Option(
        100000,
        help="Scale samples to this binsize, multiples of existing binsize only",
    ),
    cpus: int = typer.Option(
        1, help="Use multiple cores to find reference bins"
    ),
) -> None:
    """
    Create a new reference using healthy reference samples
    """
    typer.echo("Creating new reference")
    if yfrac and (0 < yfrac <= 1):
        raise typer.BadParameter(
            "Parameter --yfrac should be a positive number lower than or equal to 1"
        )

    base_path = outfile.rstrip(".npz")
    partfile = f"{base_path}_part"
    prepfile = f"{base_path}_prep.npz"

    samples = []
    logging.info("Importing data ...")
    for infile in infiles:
        logging.info("Loading: {}".format(infile))
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        binsize = int(npzdata["binsize"])
        logging.info("Binsize: {}".format(int(binsize)))
        samples.append(scale_sample(sample, binsize, binsize))

    samples = np.array(samples)
    genders, trained_cutoff = train_gender_model(samples, yfrac, yfracplot)

    if genders.count("F") < 5 and nipt:
        logging.warning(
            "A NIPT reference should have at least 5 female feti samples. Removing --nipt flag."
        )
        nipt = False
    if not nipt:
        for i, sample in enumerate(samples):
            samples[i] = gender_correct(sample, genders[i])

    total_mask, bins_per_chr = get_mask(samples)
    if genders.count("F") > 4:
        mask_f, _ = get_mask(samples[np.array(genders) == "F"])
        total_mask = total_mask & mask_f
    if genders.count("M") > 4 and not nipt:
        mask_m, _ = get_mask(samples[np.array(genders) == "M"])
        total_mask = total_mask & mask_m

    outfiles = []
    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        tmpoutfile = f"{base_path}.tmp.A.npz"
        outfiles.append(tmpoutfile)
        tool_newref_prep(
            prepfile=prepfile,
            binsize=binsize,
            samples=samples,
            gender="A",
            mask=total_mask,
            bins_per_chr=bins_per_chr,
        )
        logging.info("This might take a while ...")
        tool_newref_main(
            partfile=partfile,
            prepfile=prepfile,
            tmpoutfile=tmpoutfile,
            refsize=refsize,
            cpus=cpus,
        )
    else:
        logging.critical(
            "Provide at least 10 samples to enable the generation of a reference."
        )
        sys.exit()

    if genders.count("F") > 4:
        logging.info("Starting female gonosomal reference creation ...")
        tmpoutfile = "{}.tmp.F.npz".format(base_path)
        outfiles.append(tmpoutfile)
        tool_newref_prep(
            prepfile=prepfile,
            binsize=binsize,
            samples=samples[np.array(genders) == "F"],
            gender="F",
            mask=total_mask,
            bins_per_chr=bins_per_chr,
        )
        logging.info("This might take a while ...")
        tool_newref_main(
            partfile=partfile,
            prepfile=prepfile,
            tmpoutfile=tmpoutfile,
            refsize=refsize,
            cpus=1,
        )
    else:
        logging.warning(
            "Provide at least 5 female samples to enable normalization of female gonosomes."
        )

    if not nipt:
        if genders.count("M") > 4:
            logging.info("Starting male gonosomal reference creation ...")
            tmpoutfile = "{}.tmp.M.npz".format(base_path)
            outfiles.append(tmpoutfile)
            tool_newref_prep(
                prepfile=prepfile,
                binsize=binsize,
                samples=samples[np.array(genders) == "M"],
                gender="M",
                mask=total_mask,
                bins_per_chr=bins_per_chr,
            )
            tool_newref_main(
                partfile=partfile,
                prepfile=prepfile,
                tmpoutfile=tmpoutfile,
                refsize=refsize,
                cpus=1,
            )
        else:
            logging.warning(
                "Provide at least 5 male samples to enable normalization of male gonosomes."
            )

    tool_newref_merge(
        outfile=outfile,
        outfiles=outfiles,
        trained_cutoff=trained_cutoff,
        nipt=nipt,
    )

    typer.echo("Finished creating reference")


def predict_zscore_callback(zscore: int) -> int:
    if zscore <= 0:
        raise typer.BadParameter(
            "Parameter --zscore should be a strictly positive number"
        )
    return zscore


def predict_beta_callback(beta: float) -> float:
    if 0 <= beta < 1:
        raise typer.BadParameter(
            "Parameter --beta should be a strictly positive number lower than or equal to 1"
        )
    return beta


def predict_alpha_callback(alpha: float) -> float:
    if 0 <= alpha < 1:
        logging.critical(
            "Parameter --alpha should be a strictly positive number lower than or equal to 1"
        )
    return alpha


@app.command("predict")
def predict(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
    outid: str = typer.Argument(
        ...,
        help="Basename (w/o extension) of output files (paths are allowed, e.g. path/to/ID_1)",
    ),
    minrefbins: int = typer.Option(
        150, help="Minimum amount of sensible reference bins per target bin"
    ),
    maskrepeats: int = typer.Option(
        5,
        help="Regions with distances > mean + sd * 3 will be masked. Number of masking cycles",
    ),
    alpha: float = typer.Option(
        1e-4,
        help="p-value cut-off for calling a CBS breakpoint",
        callback=predict_alpha_callback,
    ),
    zscore: float = typer.Option(
        5,
        help="z-score cut-off for aberration calling",
        callback=predict_zscore_callback,
    ),
    beta: float = typer.Option(
        ...,
        help="When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations. Beta is a "
        "number between 0 (liberal) and 1 (conservative) and is optimally close to the purity",
        callback=predict_beta_callback,
    ),
    blacklist: Optional[Path] = typer.Option(
        ...,
        help="Blacklist that masks regions in output, structure of header-less file: chr...(/t)startpos(/t)endpos(/n)",
    ),
    gender: Gender = typer.Option(
        ...,
        help="Force WisecondorX to analyze this case as a male (M) or a female (F)",
    ),
    ylim: Tuple[int, int] = typer.Option(
        ..., help="y-axis limits for plotting. e.g. (-2,2)"
    ),
    bed: bool = typer.Option(
        True,
        help="Outputs tab-delimited .bed files, containing the most important information",
    ),
    plot: bool = typer.Option(False, help="Outputs .png plots"),
    cairo: bool = typer.Option(
        False,
        "--cairo",
        help="Uses cairo bitmap type for plotting. Might be necessary for certain setups",
    ),
):
    """
    Find copy number aberrations
    """
    if not bed and not plot:
        raise typer.BadParameter(
            "No output format selected. "
            "Select at least one of the supported output formats (--bed, --plot)"
        )

    typer.echo("Starting CNA prediction")
    typer.echo("Importing data ...")
    ref_file = np.load(reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(infile, encoding="latin1", allow_pickle=True)

    sample = sample_file["sample"].item()
    n_reads = sum([sum(sample[x]) for x in sample.keys()])

    sample = scale_sample(
        sample, int(sample_file["binsize"].item()), int(ref_file["binsize"])
    )

    predicted_gender = (
        gender
        if gender
        else predict_gender(sample, ref_file["trained_cutoff"])
    )
    if not ref_file["is_nipt"]:
        sample = gender_correct(sample, predicted_gender)
        ref_gender = predicted_gender
    else:
        ref_gender = "F"

    typer.echo("Normalizing autosomes ...")

    results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
        sample=sample,
        ref_file=ref_file,
        ref_gender="A",
        maskrepeats=maskrepeats,
    )

    if not ref_file["is_nipt"]:
        if not ref_file["has_male"] and gender == "M":
            logging.warning(
                "This sample is male, whilst the reference is created with fewer than 5 males. "
                "The female gonosomal reference will be used for X predictions. Note that these might "
                "not be accurate. If the latter is desired, create a new reference and include more "
                "male samples."
            )
            ref_gender = "F"

        elif not ref_file["has_female"] and gender == "F":
            logging.warning(
                "This sample is female, whilst the reference is created with fewer than 5 females. "
                "The male gonosomal reference will be used for XY predictions. Note that these might "
                "not be accurate. If the latter is desired, create a new reference and include more "
                "female samples."
            )
            ref_gender = "M"

    typer.echo("Normalizing gonosomes ...")

    null_ratios_aut_per_bin = ref_file["null_ratios"]
    null_ratios_gon_per_bin = ref_file["null_ratios.{}".format(ref_gender)][
        len(null_ratios_aut_per_bin) :
    ]

    results_r_2, results_z_2, results_w_2, ref_sizes_2, _, _ = normalize(
        sample=sample,
        ref_file=ref_file,
        ref_gender=ref_gender,
        maskrepeats=maskrepeats,
    )

    results_r = np.append(results_r, results_r_2)
    results_z = np.append(results_z, results_z_2) - m_z
    results_w = np.append(
        results_w * np.nanmean(results_w_2),
        results_w_2 * np.nanmean(results_w),
    )
    results_w = results_w / np.nanmean(results_w)

    if np.isnan(results_w).any() or np.isinf(results_w).any():
        logging.warning(
            "Non-numeric values found in weights -- reference too small. Circular binary segmentation and z-scoring "
            "will be unweighted "
        )
        results_w = np.ones(len(results_w))

    ref_sizes = np.append(ref_sizes, ref_sizes_2)

    null_ratios = np.array(
        [x.tolist() for x in null_ratios_aut_per_bin]
        + [x.tolist() for x in null_ratios_gon_per_bin]
    )

    results = {
        "results_r": results_r,
        "results_z": results_z,
        "results_w": results_w,
        "results_nr": null_ratios,
    }

    for result in results.keys():
        results[result] = get_post_processed_result(
            minrefbins=minrefbins,
            result=results[result],
            ref_sizes=ref_sizes,
            mask=ref_file[f"mask.{ref_gender}"],
            bins_per_chr=ref_file[f"bins_per_chr.{ref_gender}"],
        )

    log_trans(results, m_lr)

    if blacklist:
        typer.echo("Applying blacklist ...")
        apply_blacklist(
            blacklist_bed=blacklist,
            binsize=int(ref_file["binsize"]),
            results=results,
        )

    typer.echo("Executing circular binary segmentation ...")

    results["results_c"] = exec_cbs(
        alpha=alpha,
        ref_gender=ref_gender,
        binsize=int(ref_file["binsize"]),
        outid=outid,
        results=results,
    )

    if bed:
        typer.echo("Writing tables ...")
        generate_output_tables(
            outid=outid,
            binsize=int(ref_file["binsize"]),
            ref_gender=ref_gender,
            beta=beta,
            zscore=zscore,
            gender=predicted_gender,
            n_reads=n_reads,
            bins_per_chr=ref_file[f"bins_per_chr.{ref_gender}"],
            results=results,
        )

    if plot:
        typer.echo("Writing plots ...")
        exec_write_plots(
            outid=outid,
            ref_gender=ref_gender,
            binsize=int(ref_file["binsize"]),
            beta=beta,
            zscore=zscore,
            n_reads=n_reads,
            cairo=cairo,
            ylim=ylim,
            results=results,
        )

    logging.info("Finished prediction")


@app.command("gender")
def gender(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
) -> str:
    """
    Returns the gender of a .npz resulting from convert, based on a Gaussian mixture model trained during newref
    """

    ref_file = np.load(reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(infile, encoding="latin1", allow_pickle=True)
    if (
        predict_gender(
            sample_file["sample"].item(), ref_file["trained_cutoff"]
        )
        == "M"
    ):
        print("male")
    else:
        print("female")


@app.callback()
def main(
    loglevel: LogLevel = typer.Option("info", help="Set the logging level")
) -> None:
    """
    WisecondorX — An evolved WISECONDOR
    """
    logging.basicConfig(
        format="[%(levelname)s - %(asctime)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=getattr(logging, loglevel.upper(), None),
    )


if __name__ == "__main__":
    app()
