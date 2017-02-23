# WISECONDOR
### WIthin-SamplE COpy Number aberration DetectOR
#### Detect fetal trisomies and smaller CNV's in a maternal plasma sample using whole-genome data.
###### Author: [Roy Straver](https://github.com/rstraver)

The [legacy version](https://github.com/rstraver/wisecondor/tree/legacy) of the algorithm used is described in the paper:  
[Nucl. Acids Res. (2014) 42 (5): e31. doi: 10.1093/nar/gkt992](http://nar.oxfordjournals.org/content/42/5/e31).  
The current version is not published separately in a paper.

## CHANGELIST
This version of WISECONDOR is significantly different from previously published versions. While the general idea of a within-sample comparison is still applied, implementation is re-done from the ground up. If you prefer to stick with the latest version of WISECONDOR as it was, please use the [legacy branch](https://github.com/rstraver/wisecondor/tree/legacy). The last known [stable release](https://github.com/rstraver/wisecondor/releases) for the legacy version is [v2.0.1](https://github.com/rstraver/wisecondor/releases/tag/v2.0.1).  

A quick overview of changes made in general compared to the legacy version:
- Different license, especially for commercial use. Beware!
- Subparsers, now you can tell `wisecondor.py` to solve a task rather than finding the right script. This approach is alike using bamtools etc.
- Pysam, read bam files directly instead of piping samtools output into consam.
- Numpy, use lots of much faster calculations instead of slow python loops. Also replaces pickle files with npz, a zipped numpy format for smaller files.
- PCA, instead of GC-Correction, faster and solves some other effects. PCA is fit once on training data and applied to all training and test data.
- Z-score based segmentation.  
-- The z-score threshold is now dependent on the amount of tests within a sample.  
-- The single bin and set-size windowed approaches were replaced by calling the most significant sequence of bins per chromosome, until no more significant calls can be made.  
-- This method is slightly oversensitive and may call near chromosome wide read depth changes less than 1%, you are expected to filter such calls afterward. Report and plot functions have optional threshold arguments (`-mineffect`, default 1.5%) to do so.
- Specific chromosomes can now be plotted, perfect for those who want to see only 13, 18 and 21.
- Results saved in npz, obtain desired information from output directly and cleanly instead of parsing stdout. Stdout does not show results anymore.
- Runtime specific information is now saved in output files. No more need for searching through messy log files to find out what happened.
- Colorblind support, replaced the plot color scheme by something that is more easily distinguishable.
- Many awesome technical tricks to improve speed and decrease memory usage.  
-- Parallel processing is possible for reference creation. While crudely implemented for HPC usage, a single computer with multiple cores may also use this by specifying the number of available CPUs.  
-- By default bam files are converted to 50kb bins, but reference creation and analysis is run at 250kb bins. This is possible through scaling the data, allowing to run some tricks on higher resolutions without re-processing the same bam files. We do suggest using more reference samples if you aim to increase resolution. Although probably overdone, we used nearly 600 samples for our 250kb resolution reference set.

What was not changed:
- Chromosomes X and Y are still not part of analysis. I believe these require additional work for consistent results. However, you are free to implement your own solutions.

What is lost:
- DEFRAG, if you do want to use this you can obtain the script from the legacy branch.
- HMTL  Viewer, not that useful with the new segmentation approach. Feel free to re-implement the required JSON export functionality if you desire to use it anyway.

## GETTING HELP

If anything is not clear, try the [wiki](https://github.com/rstraver/wisecondor/wiki) first.
If your problem or question is still unsolved, feel free to post a question or submit a bug report on the [issue tracker](https://github.com/rstraver/wisecondor/issues).
If you do not get any response, the contact information [on my page](https://github.com/rstraver) should help you to get in touch with me.

## QUICKSTART GUIDE
Obtain the required modules:  
`pip install sklearn numpy scipy matplotlib pysam futures`  
You may need to install pip first, and depending on preferences you may add sudo at the beginning of this command.

To get to work without too much reading just use `./run.sh` and follow the directions provided.  
Binsize parameters can be changed by altering the variables near the start of the script.

## REQUIREMENTS
WISECONDOR was developed and tested using Python2.7. Using any other version may cause errors or faulty results. This version uses pysam to read .bam files created by BWA.  

The list of packages required and tested versions as reported by `pip freeze`:  
sklearn==0.0  
numpy==1.8.1  
scipy==0.17.0  
matplotlib==1.3.1  
pysam==0.8.3  
futures==3.0.5


## PREPARING FOR TESTING
To start testing for aberrations, WISECONDOR needs to learn how genomic regions behave when compared to each other. To do so, reference data should be prepared in the same way as the test sample data later on, using the same genomic reference, etc.


### External tools
Settings we used for external tools:  
`bwa aln -n 0 -k 0`  
`bwa samse -n -1`  
`picardtools sortsam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true`  

Make sure your data is sorted as specified here. Results can be odd and meaningless otherwise.

### File conversion
To convert a .bam file to a file for analysis, use the convert tool in wisecondor.py:  
`python wisecondor.py convert ./input/sample.bam ./output/sample.npz -binsize 50000`

The binsize specified here is 50kb. You can take any value here (assuming it fits in memory) and upscale data to multiples of this value later on. The final scale is determined by the scale set in the reference creation step.  
When not specified, this binsize will be set to 1mb, the default for the legacy version of WISECONDOR.


### Reference creation
All reference files should be fed into the reference creation tool. Move the reference samples into a separate folder and start the reference-build script:  
`python wisecondor.py newref ./referenceFiles/*.npz ./dataFiles/reference.npz -binsize 250000`

The binsize specified here is 250kb. Not specifying this will assume all files provided are at the desired resolution already. This value can only be a multiple of the binsize used for any sample in the reference data.  
This is the resolution samples will be scaled to during testing as well. For example, converting a sample at 50kb bins and testing it using a reference at 1mb will provide correct results at 1mb resolution.


## RUNNING TESTS
You need to feed a single converted file (npz) as explained in 'PREPARING FOR TESTING, File conversion', the path/filename where to save the output data, and the reference as explained and created in 'PREPARING FOR TESTING, Reference creation':  
`python wisecondor.py test ./testSamples/sample.npz ./testSamples/sample_out.npz ./dataFiles/reference.npz`

This creates a new npz file, in this example at `./testSamples/sample_out.npz`. This file contains the results but is not directly readable by the user. Instead it contains a lot of information that may or may not be required for reports, as well as information that should not be shown at all in diagnostics.


### Outputting results
This step is optional in the sense that you are expected to write your own report functionality. In the legacy version people needed to parse the stdout text to obtain results, making integration with other systems difficult. This version provides a simplified way to allow extracting output for third party applications.

If you want a textual report, use the report function and provide both the npz created by conversion and the output npz files:  
`python wisecondor.py report ./testSamples/sample.npz ./testSamples/sample_out.npz`


### Plotting results
The file created by the test step can be used as input for the plot tool. This tool turns the prepared data into a visualization, which can be replaced to accommodate for personal preferences without the need to make changes to the original algorithm.  
`python wisecondor.py plot ./testSamples/sample_out.npz ./testSamples/sample_plot`

Be aware that things have changed since previous versions:  
- There is no such thing as a windowed method anymore, now one size fits all
- Calls are marked if significant by z-score, independent of effect size and length
- Calls are orange if likely fetal and green if likely maternal, as determined by their effect size
- Marked opacity depends on strength of aberration
- Calls are marked with a number, this number is the z-score for the marked region it is above or below

#### Adding cytoband information
The plot tool can add a simple cytoband near the bottom of every plot for visual help. The file used for Hg19 can be found at [UCSC](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz). Unpack the archive and feed the text file to the plot tool using `-cytofile cytoBand.txt`. As long as you keep the same format you can use another cytoband file that matches your reference genome.


## TIPS

### Pitfalls
Do not use reference data from one lab to test samples from another. Every reference file, laboratory and sequencing machine has its own effect on how read depth per bin behaves. Any results obtained by combining files from different origins are unreliable.

### In general
To improve your results you probably want to change a few parameters. Most of the values used in WISECONDOR can be altered using arguments. To find out what arguments can be passed into any script, try running it with -h as argument, for example:  
`python wisecondor.py -h`  
or:  
`python wisecondor.py convert -h`

### Multi-core newref
Add `-cpus 4` to `newref` to enable multi core reference creation using 4 cores. Replace `4` with any integer matching the amount of cores you have available for a drastic speed up.

### Cluster newref
Submit reference creation parts as jobs on a compute cluster using three steps instead of `newref`. These steps consist of a single-core preparation step, a multi-core step, and finally another single-core step to combine data again. While overkill for 1mb bins, and likely overdone for 250kb bins, this approach becomes more useful when increasing resolution even further (i.e. 50kb).  
1. Create a combined file using all reference samples and applies PCA using a single CPU.  
`newrefprep ./refSamples/*.npz ./dataFiles/refprep.npz`  
2. Take the prepared file and obtains reference bins for the first 1/4th of all bins. It solves splitting chromosomes etc by itself. Replace 1 with 2, 3 and 4 to run all parts in this example. Of course, using a cluster would allow you to use any value instead of 4: start as many jobs as you want and provide the right index of the job to the script.   
`newrefpart ./dataFiles/refprep.npz ./dataFiles/refpart 1 4`  
3. Combine results from all jobs into a single reference file by providing the amount of parts used (`4`).   
`newrefpost ./dataFiles/refprep.npz ./dataFiles/refpart 4 ./dataFiles/reference.npz`  

### Multi-core vs single-core
The single core approach is currently implemented to call all steps in the multicore process in sequential order. The `-cpus` option in newref simply applies the cluster approach by itself. Therefore, all approaches use the same code and differences between single an multi core reference creation should be non-existent.

### Plotting specific chromosomes
If you do not want to see any results in the plots than aberrations on chromosomes 13 18 and 21, as some medical centers appear to prefer instead of viewing the whole genome, you can add the following options to the `plot` tool:  
`-chromosomes 13 18 21` Force only showing these three chromosomes. Replace at will with any selection of autosomal chromosomes you want.  
`-columns 1` Use the whole page width for chromosomes rather than two columns


### Customizing reports
If you open the script wisecondor.py and find `def toolReport(args):` you can see how all information can be obtained from the npz files directly. However, part of the arguments saved in an npz file are functions referenced directly and will be saved as a weird pointer in the npz. If you load one of these files without importing the functions first you may run into some errors. Either fill the missing functions variables with empty stubs or import the actual functions: `from wisetools import *`.
