WISECONDOR
(WIthin-SamplE COpy Number aberration DetectOR):
Detect fetal trisomies and smaller CNV's in a maternal plasma sample using whole-genome data.



1. REQUIREMENTS

WISECONDOR was developed and tested using Python2.7. Using any other version may cause errors or faulty results. The working version is tested using SAMTOOLS on .bam files created by BWA. Original work was based on SOAP2 files, and any SOAP2 mapped input stream should work as well although it requires additional arguments when using consam.py

Using Ubuntu, this should obtain the required packages:
	sudo apt-get install python-numpy python-matplotlib python-biopython

For any other unix based system, if easy_install is installed, you should be able to install all of these using this command in a terminal:
	sudo easy_install numpy matplotlib biopython



2. PREPARING FOR TESTING

To start testing for aberrations, WISECONDOR needs to learn how genomic regions behave when compared to eachother. To do so, reference data should be prepared in the same way as the test sample data later on, using the same genomic reference, the same GC-correction, etc.

Settings we used for external tools:
	bwa aln -n 0 -k 0
	bwa samse -n -1
	picardtools sortsam	SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
	samtools rmdup -s $INFILE - | samtools view - -q 1

To create a GC-count table for the genomic reference you use, use countgc.py:
	python countgc.py path/to/reference.fa ./ref/gccount

Applying the GC-correction should be done for any sample file you have.
	python ./gcc.py ./in/sample.pickle ./ref/gccount ./in/sample.gcc

Then, all reference files should be fed into the newref.py script. The current implementation used all samples in a single folder, so after applying all your corrections to your data, move the reference samples into a separate folder and start the reference-build script:
	python newref.py ./in/refs/ ./ref/reference



3. RUNNING TESTS

An example pipeline implementation can be found in pipeline.sh. This script will look for .pickle files in ./in/ and write output to ./out/. Additionally, details showing how to prepare your data can be found in this script but are not carried out as these require BWA and SAMTOOLS to be set up correctly.
To improve your results you probably want to change a few parameters. Most of the values used in WISECONDOR can be altered using arguments. To find out what arguments can be passed into any script, try running it with -h as argument, for example:
	python test.py -h



4. STATISTICS REPORT

If a statistics report should be created, the output files of above steps are needed. The script will look for the required files (.pickle, .bam and .bai files) in the directory indicated:
	cd to ./qc_scripts
	bash ./wc_qc.sh $INPUT_DIRECTORY
	python ./wc_makeReport.py $INPUT_DIRECTORY $OUTPUT_DIRECTORY
The first command creates a text file with statistics calculated on the bam-file, the second command creates the report using this text file, in combination with the .pdf and .txt file created by test.py.



4. PITFALLS

Do not use any of the example data for your own tests, if available. Every reference file, laboratory and sequencing machine has its own effect on how read depth per bin behaves. Any results obtained by combining files from different origins are unreliable.
