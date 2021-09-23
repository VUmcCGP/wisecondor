# Background

After extensively comparing different (shallow) whole-genome sequencing-based copy number detection tools,
including [WISECONDOR](https://github.com/VUmcCGP/wisecondor), [QDNAseq](https://github.com/ccagc/QDNAseq),
[CNVkit](https://github.com/etal/cnvkit), [Control-FREEC](https://github.com/BoevaLab/FREEC),
[BIC-seq2](http://compbio.med.harvard.edu/BIC-seq/) and
[cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html),
WISECONDOR appeared to normalize sequencing data in the most consistent way, as shown by
[our paper](https://www.ncbi.nlm.nih.gov/pubmed/30566647). Nevertheless, WISECONDOR has limitations:
Stouffer's Z-score approach is error-prone when dealing with large amounts of aberrations, the algorithm
is extremely slow (24h) when operating at small bin sizes (15 kb), and sex chromosomes are not part of the analysis.
Here, we present WisecondorX, an evolved WISECONDOR that aims at dealing with previous difficulties, resulting
in overall superior results and significantly lower computing times, allowing daily diagnostic use. WisecondorX is
meant to be applicable not only to NIPT, but also gDNA, PGT, FFPE, LQB, ... etc.

# Manual

## Mapping

We found superior results through WisecondorX when using [bowtie2](https://github.com/BenLangmead/bowtie2) as a mapper.
Note that it is important that **no** read quality filtering is executed prior to running WisecondorX: this software
requires low-quality reads to distinguish informative bins from non-informative ones.

## WisecondorX

### Installation

Stable releases can be installed through pip install. This option ascertains the latest version is
downloaded, however, it does not install R [dependencies](#dependencies).

```bash

pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX
```

Alternatively, [Conda](https://conda.io/docs/) additionally installs all necessary [depedencies](#dependencies),
however, the latest version might not be downloaded.

```bash

conda install -f -c conda-forge -c bioconda wisecondorx
```

### Running WisecondorX

There are three main stages (converting, reference creating and predicting) when using WisecondorX:

- Convert aligned reads to .npz files (for both reference and test samples)
- Create a reference (using reference .npz files)
  - **Important notes**
    - Automated gender prediction, required to consistently analyze sex chromosomes, is based on a Gaussian mixture
      model. If few samples (<20) are included during reference creation, or not both male and female samples (for
      NIPT, this means male and female feti) are represented, this process might not be accurate. Therefore,
      alternatively, one can manually tweak the [`--yfrac`](#stage-2-create-reference) parameter.
    - It is of paramount importance that the reference set consists of exclusively negative control samples that
      originate from the same sequencer, mapper, reference genome, type of material, ... etc, as the test samples.
      As a rule of thumb, think of all laboratory and _in silico_ steps: the more sources of bias that can be omitted,
      the better.
    - Try to include at least 50 samples per reference. The more the better, yet, from 500 on it is unlikely to
      observe additional improvement concerning normalization.
- Predict copy number alterations (using the reference file and test .npz cases of interest)

### Stage (1) Convert aligned reads (bam/cram) to .npz

```bash
 $ WisecondorX convert --help
Usage: WisecondorX convert [OPTIONS] INFILE OUTFILE

  Convert and filter a aligned reads to .npz

Arguments:
  INFILE   aligned reads input for conversion  [required]
  OUTFILE  Output .npz file  [required]

Options:
  --reference PATH   Fasta reference to be used during cram conversion
  --binsize INTEGER  Bin size (bp)  [default: 5000]
  --normmdup         Avoid remove duplicates  [default: False]
  --help             Show this message and exit.

```

&rarr; Bash recipe at `docs/include/pipeline/convert.sh`

### Stage (2) Create reference

```bash
$ WisecondorX newref --help 
Usage: WisecondorX newref [OPTIONS] INFILES... OUTFILE

  Create a new reference using healthy reference samples

Arguments:
  INFILES...  Path to all reference data files (e.g. path/to/reference/*.npz)
              [required]

  OUTFILE     Path and filename for the reference output (e.g.
              path/to/myref.npz)  [required]


Options:
  --nipt             Use flag for NIPT  [default: False]
  --yfrac FLOAT      Use to manually set the y read fraction cutoff, which
                     defines gender  [required]

  --yfracplot PATH   Path to yfrac .png plot for --yfrac optimization (e.g.
                     path/to/plot.png); software will stop after plotting
                     after which --yfrac can be set manually

  --refsize INTEGER  Amount of reference locations per target  [default: 300]
  --binsize INTEGER  Scale samples to this binsize, multiples of existing
                     binsize only  [default: 100000]

  --cpus INTEGER     Use multiple cores to find reference bins  [default: 1]
  --help             Show this message and exit.

```
&rarr; Bash recipe at `docs/include/pipeline/newref.sh`

### Stage (3) Predict copy number alterations

```bash
$ WisecondorX predict --help
Usage: WisecondorX predict [OPTIONS] INFILE REFERENCE OUTID

  Find copy number aberrations

Arguments:
  INFILE     .npz input file  [required]
  REFERENCE  Reference .npz, as previously created with newref  [required]
  OUTID      Basename (w/o extension) of output files (paths are allowed, e.g.
             path/to/ID_1)  [required]


Options:
  --minrefbins INTEGER         Minimum amount of sensible reference bins per
                               target bin  [default: 150]

  --maskrepeats INTEGER        Regions with distances > mean + sd * 3 will be
                               masked. Number of masking cycles  [default: 5]

  --alpha FLOAT                p-value cut-off for calling a CBS breakpoint
                               [default: 0.0001]

  --zscore FLOAT               z-score cut-off for aberration calling
                               [default: 5]

  --beta FLOAT                 When beta is given, --zscore is ignored and a
                               ratio cut-off is used to call aberrations. Beta
                               is a number between 0 (liberal) and 1
                               (conservative) and is optimally close to the
                               purity  [required]

  --blacklist PATH             Blacklist that masks regions in output,
                               structure of header-less file:
                               chr...(/t)startpos(/t)endpos(/n)  [required]

  --gender [F|M]               Force WisecondorX to analyze this case as a
                               male (M) or a female (F)  [required]

  --ylim <INTEGER INTEGER>...  y-axis limits for plotting. e.g. (-2,2)
                               [required]

  --bed / --no-bed             Outputs tab-delimited .bed files, containing
                               the most important information  [default: True]

  --plot / --no-plot           Outputs .png plots  [default: False]
  --cairo                      Uses cairo bitmap type for plotting. Might be
                               necessary for certain setups  [default: False]

  --help                       Show this message and exit.
```
&rarr; Bash recipe at `docs/include/pipeline/predict.sh`

### Additional functionality

```bash
$ WisecondorX gender --help 
Usage: WisecondorX gender [OPTIONS] INFILE REFERENCE

  Returns the gender of a .npz resulting from convert, based on a Gaussian
  mixture model trained during newref

Arguments:
  INFILE     .npz input file  [required]
  REFERENCE  Reference .npz, as previously created with newref  [required]

Options:
  --help  Show this message and exit.


```

Returns gender according to the reference.

# Parameters

The default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x coverage) and reference bin
sizes ranging from 50 to 500 kb.

# Underlying algorithm

To understand the underlying algorithm, I highly recommend reading
[Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809); and its update shortly introduced in
[Huijsdens-van Amsterdam et al (2018)](https://www.nature.com/articles/gim201832.epdf). Numerous adaptations to this
algorithm have been made, yet the central normalization principles remain. Changes include e.g. the inclusion of a gender
prediction algorithm, gender handling prior to normalization (ultimately enabling X and Y predictions), between-sample
Z-scoring, inclusion of a weighted circular binary segmentation algorithm and improved codes for obtaining tables and
plots.

# Interpretation results

## Plots

Every dot represents a bin. The dots range across the X-axis from chromosome 1 to X (or Y, in case of a male). The
vertical position of a dot represents the ratio between the observed number of reads and the expected number of reads,
the latter being the 'normal' state. As these values are log2-transformed, copy neutral dots should be close-to 0.
Importantly, notice that the dots are always subject to Gaussian noise. Therefore, segments, indicated by horizontal
white lines, cover bins of predicted equal copy number. The size of the dots represents the variability at the reference
set. Thus, the size increases with the certainty of an observation. The same goes for the line width of the segments.
Vertical grey bars represent the blacklist, which matches mostly hypervariable loci and repeats. Finally, the horizontal
colored dotted lines show where the constitutional 1n and 3n states are expected (when constitutional DNA is at 100%
purity). Often, an aberration does not reach these limits, which has several potential causes: depending on your type
of analysis, the sample could be subject to tumor fraction, fetal fraction, a mosaicism, ... etc. Sometimes, the
segments do surpass these limits: here it's likely you are dealing with 0n, 4n, 5n, 6n, ...

## Tables

### ID_bins.bed

This file contains all bin-wise information. When data is 'NaN', the corresponding bin is included in the blacklist.
The Z-scores are calculated as default using the within-sample reference bins as a null set.

### ID_segments.bed

This file contains all segment-wise information. A combined Z-score is calculated using a between-sample Z-scoring
technique (the test case vs the reference cases).

### ID_aberrations.bed

This file contains aberrant segments, defined by the [`--beta`](#stage-3-predict-copy-number-alterations) or
[`--zscore`](#stage-3-predict-copy-number-alterations) parameters.

### ID_statistics.bed

This file contains some interesting statistics (per chromosome and overall). The definition of the Z-scores matches the one from
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
  - matplotlib (v2.2.3)
  - typer (v0.4.0)

And of course, other versions are very likely to work as well.
