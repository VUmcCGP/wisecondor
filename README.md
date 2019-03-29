# Background  
After extensively comparing different (shallow) WGS-based CNA tools, including [WISECONDOR](https://github.com/VUmcCGP/wisecondor),
[QDNAseq](https://github.com/ccagc/QDNAseq), [CNVkit](https://github.com/etal/cnvkit), [Control-FREEC](https://github.com/BoevaLab/FREEC),
[BIC-seq2](http://compbio.med.harvard.edu/BIC-seq/) and [cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html),
WISECONDOR appeared to normalize sWGS copy number data in the most consistent way, as shown by [our paper](https://www.ncbi.nlm.nih.gov/pubmed/30566647). Nevertheless,
as is the case with every tool, WISECONDOR has limitations of its own: the Stouffer's z-score approach is error-prone when
dealing with large amounts of aberrations, the algorithm is extremely slow (24h) when using small bin sizes (15 kb) and
sex chromosomes are not included in the analysis. Here, we present WisecondorX, an evolved WISECONDOR that aims at dealing with
previous difficulties. Main adaptations include the additional (and consistent) analysis of the X and Y chromosomes (for
NIPT exclusively X), a weighted CBS-based segmentation technique, a gender prediction algorithm and a custom plotter, resulting
in overall superior results and significantly lower computing times, allowing daily diagnostic use. WisecondorX is meant
to be applicable not only to NIPT, but also to gDNA, PGD, FFPE, LQB, ... etc.

# Manual

## Mapping

We found superior results through WisecondorX when using the following mapping and deduplication strategy.

```bash

bowtie2 --local -p n --fast-local -x index -U input | bamsormadup inputformat=sam threads=n SO=coordinate outputformat=bam indexfilename=sample.bam.bai > sample.bam
```

I would recommend using the latest version of the human reference genome (GRCh38). Note that it is very important that **no** read
quality filtering is executed prior the running WisecondorX: this software requires low-quality reads to distinguish informative
bins from non-informative ones.

## WisecondorX

### Installation

Stable releases can be installed using [Conda](https://conda.io/docs/). This option takes care of all necessary
dependencies.
```bash

conda install -f -c conda-forge -c bioconda wisecondorx
```

Alternatively, WisecondorX can be installed through pip install. This option ascertains the latest version is downloaded,
yet it does not install R dependencies.  
```bash

pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX
```

### Running WisecondorX

There are three main stages (converting, reference creating & predicting) for using WisecondorX:  
- Convert .bam files to .npz files (both reference and test samples)  
- Create a reference (using reference .npz files)  
    - **Important notes**
        - WisecondorX will internally generate a male and female gonosomal reference. It is advised that both male and female
        samples (for NIPT, this means male and female feti) are represented in the reference set. Furthermore, they
        should be more or less equally included.  
        - Gender prediction is based on a Gaussian mixture model of the Y-read fraction. This means that, if only few samples (<30)
        are included during reference creation, this process might not be accurate, and gonosomal predictions could fail.  
        - It is of paramount importance that the reference set consists of exclusively healthy samples that originate from
        the same sequencer, mapper, reference genome, type of material, ... etc, as the test samples. As a rule of thumb,
        think of all laboratory and in silico pre-processing steps: the more sources of bias that can be omitted,
        the better.  
        - Try to include at least 100 samples per reference. The more the better, yet, from 500 on it is unlikely to observe
        additional improvement concerning normalization.  
- Predict CNAs (using the reference and test .npz cases of interest)  

### Stage (1) Convert .bam to .npz

```bash

WisecondorX convert input.bam output.npz [--optional arguments]
```

<br>Optional argument <br><br> | Function  
:--- | :---  
`--binsize x` | Size per bin in bp, the reference bin size should be a multiple of this value (default: x=5e3)  
`--paired` | Enables conversion for paired-end reads  


&rarr; Bash recipe (example for NIPT) at `./pipeline/convert.sh`

### Stage (2) Create reference

```bash

WisecondorX newref reference_input_dir/*.npz reference_output.npz [--optional arguments]
```

<br>Optional argument<br><br> | Function
:--- | :---  
`--nipt` | **Always include this flag for the generation of a NIPT reference**  
`--binsize x` | Size per bin in bp, defines the resolution of the output (default: x=1e5)  
`--refsize x` | Amount of reference locations per target (default: x=300)  
`--cpus x` | Number of threads requested (default: x=1)  

&rarr; Bash recipe (example for NIPT) at `./pipeline/newref.sh`

### Stage (3) Predict CNAs  

```bash

WisecondorX predict test_input.npz reference_input.npz output_id [--optional arguments]
```
  
<br>Optional argument &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Function  
:--- | :---  
`--minrefbins x` | Minimum amount of sensible reference bins per target bin (default: x=150)  
`--maskrepeats x` | Regions with distances > mean + sd * 3 in the reference will be masked, number of masking cycles (default: x=5)  
`--zscore x` | Z-score cut-off to call segments as aberrations (default: x=3)  
`--alpha x` | P-value cut-off for calling a CBS breakpoint (default: x=1e-4)  
`--beta x` | When beta is given, --zscore is ignored. Beta sets a ratio cut-off for aberration calling. It's a number between 0 (liberal) and 1 (conservative) and is optimally close to the purity (e.g. fetal/tumor fraction) (default: x=None)  
`--blacklist x` | Blacklist that masks additional regions in output, requires header-less .bed file. This is particularly useful when the reference set is a too small to recognize some obvious regions (such as centromeres; example at `./example.blacklist/centromere.hg38.txt`) (default: x=None)  
`--gender x` | Force WisecondorX to analyze this file as a male (M) or female (F). Useful when dealing with a whole-chromosome Y-deletion (choices: x=F or x=M)
`--ylim [a,b]` | Force WisecondorX to use y-axis interval [a,b] for plotting, e.g. [-2,2]. (default: a and b selected based on data)  
`--bed` | Outputs tab-delimited .bed files (trisomy 21 NIPT example at `./example.bed`), containing all necessary information  **(\*)**  
`--plot` | Outputs custom .png plots (trisomy 21 NIPT example at `./example.plot`), directly interpretable  **(\*)**  
`--ciaro` | Some operating systems require the cairo bitmap type to write plots  

<sup>**(\*)** At least one of these output formats should be selected</sup>  

&rarr; Bash recipe (example for NIPT) at `./pipeline/predict.sh`

### Additional functionality

```bash

WisecondorX gender test_input.npz reference_input.npz
```

Returns gender based Gaussian bimodel trained during reference creation.

# Parameters

The default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x depth; sWGS) and reference bin sizes 
ranging from 50 to 500 kb. When increasing the reference bin size (`--binsize`), I recommend lowering the reference locations 
per target (`--refsize`) and the minimum amount of sensible reference bins per target bin (`--minrefbins`). Further note that a
reference bin size lower than 15 kb is not advised.  
**Important note**  
Concerning the vast majority of applications, the `--alpha` parameter should not be tweaked. The `--beta` parameter on the contrary
should depend on your type of analysis. For NIPT, its default value should be fine. However, for gDNA, when mosaicisms are of no interest,
it could be increased to its maximum, being 1. When the fetal (NIPT) or tumor (LQB, fresh material, FFPE, ...) fraction is known, this parameter is optimally
close to this fraction. If you have any doubts about this argument, a default `--beta` should still be fine when a good and large reference set was created,
irrespective of the type of analysis.  
If the interest is research, and not diagnostics, feel free to use the z-scores by defining a cut-off of interest.

# Underlying algorithm

To understand the underlying algorithm, I highly recommend reading [Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809); and its
update shortly introduced in [Huijsdens-van Amsterdam et al (2018)](https://www.nature.com/articles/gim201832.epdf).
Numerous adaptations to this algorithm have been made, yet the central principles remain. Changes include e.g. the inclusion of a gender
prediction algorithm, gender handling prior to normalization (ultimately enabling X and Y predictions), extensive
blacklisting options, inclusion of a weighted CBS algorithm and improved codes for outputting tables and plots.

# Interpretation results

## Plots

Every dot represents a bin. The dots range across the X-axis from chromosome 1 to X (or Y, in case of a male). The vertical 
position of a dot represents the ratio between the observed number of reads and the expected number of reads, the latter being 
the 'healthy' state. As these values are log2-transformed, 'healthy dots' should be close-to 0. Importantly, notice that 
the dots are always subject to Gaussian noise. Therefore, segments, indicated by horizontal grey bars, cover bins of 
predicted equal copy number. The size of the dots represent the variability at the reference set. Thus, the size increases 
with the certainty of an observation. The same goes for line width of the segments. 
Vertical grey bars represent the blacklist, which will match hypervariable loci and repeats. Finally, the vertical colored
dotted lines show where the constitutional 1n and 3n states are expected (when constitutional DNA is at 100% purity).
Often, an aberration does not surpass these limits, which has several potential causes: depending on your type of analysis, 
the sample could be subject to tumor fraction, fetal fraction, mosaicisms, etc ... Sometimes, the segments do surpass these 
limits: here it's likely you are dealing with 4n, 5n, 6n, ...

## Tables

### ID_bins.bed

This file contains all bin-wise information. When data is 'NaN', the corresponding bin is included in the blacklist. 
The Z-scores are calculated as default using the within-sample reference bins as a null set.

### ID_segments.bed

This file contains all segment-wise information. A combined Z-score is calculated using a between-sample z-scoring
technique (the test case vs the reference cases).

### ID_aberrations.bed

This file contains segments with a ratio surpassing a certain cutoff value, defined by the `--beta` parameter.

### ID_chr_statistics.bed

This file contains some interesting statistics for each chromosome. The definition of the Z-scores matches the one from 
the 'ID_segments.bed'. Particularly interesting for NIPT.

# Dependencies

- R (v3.4) packages
    - jsonlite (v1.5)
- R Bioconductor (v3.5) packages
    - DNAcopy (v1.50.1)
- Python (v3.6) libraries
	- scipy (v1.1.0)
    - scikit-learn (v0.20.0)
    - pysam (v0.15.1)
    - numpy (v1.15.2)

And of course, other versions are very likely to work as well.  
