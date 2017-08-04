#!/usr/bin/env python
# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (github.com/rstraver)
#
# This file is part of WISECONDOR
# WISECONDOR is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



import argparse
import sys
from scipy.stats import norm
import math
import os

from wisetools import *

def toolConvert(args):
	converted, qual_info = convertBam(args.infile, binsize=args.binsize, minShift=args.retdist, threshold=args.retthres)
	np.savez_compressed(args.outfile,
						arguments=vars(args),
						runtime=getRuntime(),
						sample=converted,
						quality=qual_info)
	print 'Conversion finished'


def toolNewref(args):
	splitPath = list(os.path.split(args.outfile))
	if splitPath[-1][-4:] == '.npz':
		splitPath[-1] = splitPath[-1][:-4]
	basePath = os.path.join(splitPath[0], splitPath[1])

	# Add single thread information used for parallel processing functions
	args.prepfile = basePath + "_prep.npz"
	args.partfile = basePath + "_part"
	args.parts = max(args.parts, args.cpus)
	#args.part = [1, args.parts]

	# Run parallel processing functions
	if not os.path.isfile(args.prepfile):
		toolNewrefPrep(args)

	# Use multiple cores if requested
	if args.cpus != 1:
		import concurrent.futures
		import copy
		with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpus) as executor:
			for part in xrange(1, args.parts + 1):
				if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
					thisArgs = copy.copy(args)
					thisArgs.part = [part, args.parts]
					executor.submit(toolNewrefPart, thisArgs)
			executor.shutdown(wait=True)
	else:
		for part in xrange(1, args.parts + 1):
			if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
				args.part = [part, args.parts]
				toolNewrefPart(args)

	# Put it all together
	toolNewrefPost(args)

	# Remove parallel processing temp data
	os.remove(args.prepfile)
	for part in xrange(1, args.parts + 1):
		os.remove(args.partfile + '_' + str(part) + '.npz')


def toolNewrefPrep(args):
	samples = []
	binsizes = set()
	for infile in args.infiles:  # glob.glob(args.infiles):
		print 'Loading:', infile,
		npzdata = np.load(infile)
		print ' \tbinsize:', int(npzdata['arguments'].item()['binsize'])
		samples.append(scaleSample(npzdata['sample'].item(), npzdata['arguments'].item()['binsize'], args.binsize))
		binsizes.add(npzdata['arguments'].item()['binsize'])

	if args.binsize is None and len(binsizes) != 1:
		print 'ERROR: There appears to be a mismatch in binsizes in your dataset:', binsizes
		print 'Either remove the offending sample or use -binsize to scale all samples'
		exit(1)

	binsize = args.binsize
	if args.binsize is None:
		binsize = binsizes.pop()

	maskedData, chromosomeBins, mask = toNumpyArray(samples)
	del samples
	maskedChromBins = [sum(mask[sum(chromosomeBins[:i]):sum(chromosomeBins[:i]) + x]) for i, x in
						enumerate(chromosomeBins)]
	maskedChromBinSums = [sum(maskedChromBins[:x + 1]) for x in range(len(maskedChromBins))]
	correctedData, pca = trainPCA(maskedData)
	np.savez_compressed(args.prepfile,
						arguments=vars(args),
						runtime=getRuntime(),
						binsize=binsize,
						chromosomeBins=chromosomeBins,
						maskedData=maskedData,
						mask=mask,
						maskedChromBins=maskedChromBins,
						maskedChromBinSums=maskedChromBinSums,
						correctedData=correctedData,
						pca_components=pca.components_,
						pca_mean=pca.mean_)


def toolNewrefPart(args):
	if args.part[0] > args.part[1]:
		print 'ERROR: Part should be smaller or equal to total parts:', args.part[0], '>', args.part[
			1], 'is wrong'
		exit(1)
	if args.part[0] < 0:
		print 'ERROR: Part should be at least zero:', args.part[0], '<', 0, 'is wrong'
		exit(1)

	npzdata = np.load(args.prepfile)
	correctedData = npzdata['correctedData']
	maskedChromBins = npzdata['maskedChromBins']
	maskedChromBinSums = npzdata['maskedChromBinSums']

	indexes, distances = getReference(correctedData, maskedChromBins, maskedChromBinSums,
										selectRefAmount=args.refsize, part=args.part[0], splitParts=args.part[1])

	np.savez_compressed(args.partfile + '_' + str(args.part[0]) + '.npz',
						arguments=vars(args),
						runtime=getRuntime(),
						indexes=indexes,
						distances=distances)


def toolNewrefPost(args):
	# Load prep file data
	npzdata = np.load(args.prepfile)
	maskedChromBins = npzdata['maskedChromBins']
	chromosomeBins = npzdata['chromosomeBins']
	mask = npzdata['mask']
	pca_components = npzdata['pca_components']
	pca_mean = npzdata['pca_mean']
	binsize = npzdata['binsize'].item()

	# Load and combine part file data
	bigIndexes = []
	bigDistances = []
	for part in xrange(1, args.parts + 1):  # glob.glob(args.infiles):
		infile = args.partfile + '_' + str(part) + '.npz'
		print 'Loading:', infile
		npzdata = np.load(infile)
		bigIndexes.extend(npzdata['indexes'])
		bigDistances.extend(npzdata['distances'])

		print part, npzdata['indexes'].shape

	indexes = np.array(bigIndexes)
	distances = np.array(bigDistances)

	np.savez_compressed(args.outfile,
						arguments=vars(args),
						runtime=getRuntime(),
						binsize=binsize,
						indexes=indexes,
						distances=distances,
						chromosome_sizes=chromosomeBins,
						mask=mask,
						masked_sizes=maskedChromBins,
						pca_components=pca_components,
						pca_mean=pca_mean)


# Most of this "tool" is actually basic functionality and should be in wisetools.py instead
def toolTest(args):
	# Reference data handling
	mask_list = []
	referenceFile = np.load(args.reference)
	binsize = referenceFile['binsize'].item()
	indexes = referenceFile['indexes']
	distances = referenceFile['distances']
	chromosome_sizes = referenceFile['chromosome_sizes']
	mask = referenceFile['mask']
	mask_list.append(mask)
	masked_sizes = referenceFile['masked_sizes']
	maskedChromBinSums = [sum(masked_sizes[:x + 1]) for x in range(len(masked_sizes))]

	pca_mean = referenceFile['pca_mean']
	pca_components = referenceFile['pca_components']

	del referenceFile

	# Test sample data handling
	sampleFile = np.load(args.infile)
	sample = sampleFile['sample'].item()
	sampleBinSize = sampleFile['arguments'].item()['binsize']
	sample = scaleSample(sample, sampleBinSize, binsize)
	# sample = np.load(args.infile)['sample'].item()

	testData = toNumpyRefFormat(sample, chromosome_sizes, mask)
	testData = applyPCA(testData, pca_mean, pca_components)
	optimalCutoff, cutOffMask = getOptimalCutoff(distances, 3)  # TODO: Stuff this mask into repeatTest instead of value

	num_tests = sum(masked_sizes)
	z_threshold = norm.ppf(1 - 1. / (num_tests * 0.5 * args.multitest))

	if args.minzscore is not None:
		z_threshold = args.minzscore
	print 'Per bin z-score threshold for first testing cycles:', z_threshold

	testCopy = np.copy(testData)
	resultsZ, resultsR, refSizes, stdDevAvg = repeatTest(testCopy, indexes, distances, masked_sizes, maskedChromBinSums, optimalCutoff, z_threshold, args.repeats)
	print 'ASDES:', stdDevAvg, '\nAASDEF:', stdDevAvg * z_threshold

	# Get rid of infinite values caused by having no reference bins or only zeros in the reference
	infinite_mask = (refSizes >= args.minrefbins) #np.isfinite(propZ)
	mask_list.append(infinite_mask)
	cleanedR = resultsR[infinite_mask]
	cleanedZ = resultsZ[infinite_mask]

	cleanedBinSums = [np_sum(infinite_mask[:val]) for val in maskedChromBinSums]
	cleanedBins = [cleanedBinSums[i] - cleanedBinSums[i - 1] for i in range(1, len(cleanedBinSums))]
	cleanedBins.insert(0, cleanedBinSums[0])

	# And now for Stouffers approach
	timeStartTest = time.time()
	shifter_chromsizes = [0]
	shifter_chromsizes.extend(chromosome_sizes)
	shifter = np.ones(cleanedZ.shape, dtype=bool)
	shifter_inflated = inflateArrayMulti(shifter, mask_list)
	stouffCalls = []
	chromWide = []

	for i in [x-1 for x in args.chromosomes]:
		start = sum(cleanedBins[:i])
		end = sum(cleanedBins[:i + 1])
		zTriangle = fillTriMin(cleanedZ[start:end],cleanedR[start:end],args.mineffectsize)
		chromWide.append(zTriangle.getValue(0,end-start-1))
		stoseg = zTriangle.segmentTri(z_threshold, 3)
		for seg in stoseg:

			perc = (seg[1][1] - seg[1][0]) / float(end - start) * 100
			shifter_start = sum(chromosome_sizes[:i])
			filledBins = 0
			while filledBins <= seg[1][0]:
				filledBins += shifter_inflated[shifter_start] != 0
				shifter_start += 1
			shifter_start -= 1
			shifter_end = shifter_start
			while filledBins <= seg[1][1]:
				filledBins += shifter_inflated[shifter_end] != 0
				shifter_end += 1
			shifter_start -= sum(chromosome_sizes[:i])
			shifter_end -= sum(chromosome_sizes[:i])
			#print i + 1, str(perc)[:5], '%\t', end - start, seg[1], seg[0], np_median(cleanedR[start + seg[1][0]:start + seg[1][
			#	1]+1]) - 1, shifter_start, shifter_end, shifter_start * binsize, shifter_end * binsize
			stouffCalls.append([i + 1, shifter_start, shifter_end,
								seg[0], np_median(cleanedR[start + seg[1][0]:start + seg[1][1]+1]) - 1])
	print 'Time spent on obtaining stouffers z-scores:', int(time.time() - timeStartTest), 'seconds'

	resultsZ = []
	resultsR = []
	inflatedZ = inflateArrayMulti(cleanedZ, mask_list)
	inflatedR = inflateArrayMulti(cleanedR - 1, mask_list)
	for chrom in xrange(len(chromosome_sizes)):
		chromData = inflatedZ[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
		resultsZ.append(chromData)
		chromData = inflatedR[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
		resultsR.append(chromData)

	np.savez_compressed(args.outfile,
						arguments=vars(args),
						runtime=getRuntime(),
						binsize=binsize,
						results_r=resultsR,
						results_z=resultsZ,
						results_cwz=chromWide,
						results_calls=stouffCalls,
						threshold_z=z_threshold,
						asdef=stdDevAvg,
						aasdef=stdDevAvg*z_threshold)
	exit(0)


def toolPlot(args):
	resultFile = np.load(args.infile)
	name = args.infile.split('/')[-1].split('.')[0]

	plotLines(resultFile['results_z'], resultFile['results_calls'], resultFile['threshold_z'],
		sampleName=name,
		minEffect=args.mineffect,
		binsize=resultFile['binsize'],
		cytoFile=args.cytofile,
		chromosomes=args.chromosomes,
		columns=args.columns,
		size=args.size).savefig(args.outfile+'_z.'+args.filetype)
	# plotLines(resultFile['results_r'], [], 0.1,
	# 	sampleName=name,
	# 	cytoFile=args.cytofile,
	# 	chromosomes=args.chromosomes,
	# 	columns=args.columns,
	# 	size=args.size).savefig(args.outfile+'_r.'+args.filetype)


def toolReport(args):
	test_file = np.load(args.testfile)
	result_file = np.load(args.resultfile)
	#print result_file.keys()
	binSize = result_file['binsize'].item()

	print '\n# Arguments used in convert: #'
	test_args = test_file['arguments'].item()
	for key in test_args:
		print key,'=',test_args[key]
	print '\n# Arguments used in test: #'
	result_args = result_file['arguments'].item()
	for key in result_args:
		print key,'=',result_args[key]

	test_quality = test_file['quality'].item()
	print '\n# BAM information: #'
	print 'Reads mapped:  \t',test_quality['mapped']
	print 'Reads unmapped:\t',test_quality['unmapped']
	print 'Reads nocoord: \t',test_quality['no_coordinate']
	print 'Reads rmdup:   \t',test_quality['filter_rmdup']
	print 'Reads lowqual: \t',test_quality['filter_mapq']

	#print 'Reads filtered:\t',test_quality['pre_retro']-test_quality['post_retro']
	print '\n# RETRO filtering: #'
	print 'Reads in:     \t',test_quality['pre_retro']-test_quality['no_coordinate']+test_quality['filter_rmdup']+test_quality['filter_mapq']
	print 'Reads removed:\t',test_quality['pre_retro']-test_quality['no_coordinate']+test_quality['filter_rmdup']+test_quality['filter_mapq']-test_quality['post_retro']
	print 'Reads out:    \t',test_quality['post_retro']

	print '\n# Z-Score checks: #'
	print 'Z-Score used:\t',"{:.2f}".format(result_file['threshold_z'].item())
	print 'AvgStdDev:   \t',"{:.2f}".format(result_file['asdef']*100)+'%'
	print 'AvgAllStdDev:\t',"{:.2f}".format(result_file['aasdef']*100)+'%'

	print '\n# Test results: #'
	print 'z-score\teffect\tmbsize\tlocation'
	for result in result_file['results_calls']:
		if args.mineffect < abs(result[4]*100):
			print "{:.2f}\t{:.2f}\t{:.2f}\t{:.0f}:{:.0f}-{:.0f}".format(result[3],result[4]*100,(result[2]-result[1]+1) * binSize/1e6,result[0],result[1]*binSize,(result[2]+1)*binSize)


def main():
	parser = argparse.ArgumentParser(
		description="WISECONDOR (WIthin-SamplE COpy Number aberration DetectOR)")
	subparsers = parser.add_subparsers()

	# File conversion
	parser_convert = subparsers.add_parser('convert',
		description='Convert and filter a bam file to an npz')
	parser_convert.add_argument('infile',
		type=str,
		help='Bam input file for conversion')
	parser_convert.add_argument('outfile',
		type=str,
		help='File to write binned information to')
	parser_convert.add_argument('-binsize',
		type=float, default=1e6,
		help='Size per bin in bp')
	parser_convert.add_argument('-retdist',
		type=int, default=4,
		help='Maximum amount of base pairs difference between sequential reads to consider them part of the same tower')
	parser_convert.add_argument('-retthres',
		type=int, default=4,
		help='Threshold for when a group of reads is considered a tower and will be removed')
	parser_convert.set_defaults(func=toolConvert)

	# New reference creation
	parser_newref = subparsers.add_parser('newref',
		description='Create a new reference using healthy reference samples')
	parser_newref.add_argument('infiles',
		type=str, nargs='*',
		help='Path and all to reference data files (i.e. ./reference/*.npz)')
	parser_newref.add_argument('outfile',
		type=str,
		help='Path and filename for the reference output (i.e. ./reference/myref.npz)')
	parser_newref.add_argument('-refsize',
		type=int, default=100,
		help='Amount of reference locations per target')
	parser_newref.add_argument('-binsize',
		type=int, default=None,
		help='Try to scale samples to this binsize, multiples of existing binsize only')
	parser_newref.add_argument('-cpus',
		type=int, default=1,
		help='EXPERIMENTAL: Use multiple cores to find reference bins')
	parser_newref.add_argument('-parts',
		type=int, default=1,
		help='Split reference finding in this many jobs, only used if > cpus')
	parser_newref.set_defaults(func=toolNewref)

	# New reference creation prep
	parser_newrefprep = subparsers.add_parser('newrefprep',
		description='Prepare creation of new reference split over several processes')
	parser_newrefprep.add_argument('infiles',
		type=str, nargs='*',
		help='Path and all to reference data files (i.e. ./reference/*.npz)')
	parser_newrefprep.add_argument('prepfile',
		type=str,
		help='Path and filename for the prep output (i.e. ./reference/myref_prep.npz)')
	parser_newrefprep.add_argument('-binsize',
		type=int, default=None,
		help='Try to scale samples to this binsize, multiples of existing binsize only')
	parser_newrefprep.set_defaults(func=toolNewrefPrep)

	# New reference creation by part
	parser_newrefpart = subparsers.add_parser('newrefpart',
		description='Creation of new reference split over several processes')
	parser_newrefpart.add_argument('prepfile',
		type=str,
		help='Path and filename for the prep file  (i.e. ./reference/myref_prep.npz)')
	parser_newrefpart.add_argument('partfile',
		type=str,
		help='Path and basename for the reference part output (i.e. ./reference/myref_part), receives _m.npz extension from part argument')
	parser_newrefpart.add_argument('part',
		type=int, default=[0, 1], nargs=2,
		help='Used to split data for multi core reference making, represents part m out of n parts, appends _m.npz to file name')
	parser_newrefpart.add_argument('-refsize',
		type=int, default=100,
		help='Amount of reference locations per target')
	parser_newrefpart.set_defaults(func=toolNewrefPart)

	# New reference creation post
	parser_newrefpost = subparsers.add_parser('newrefpost',
		description='Combine creation of new reference split over several processes')
	parser_newrefpost.add_argument('prepfile',
		type=str,
		help='Path and filename for the prep file  (i.e. ./reference/myref_prep.npz)')
	parser_newrefpost.add_argument('partfile',
		type=str,
		help='Path and basename to reference part data files (i.e. ./reference/myref_part), appends _m.npz depending on parts')
	parser_newrefpost.add_argument('parts',
		type=int, default=1,
		help='Used combine all data parts, represents n parts previously specified')
	parser_newrefpost.add_argument('outfile',
		type=str,
		help='Path and filename for the reference output (i.e. ./reference/myref.npz)')
	parser_newrefpost.set_defaults(func=toolNewrefPost)

	# Sample testing
	parser_test = subparsers.add_parser('test',
		description='Test sample for Copy Number Aberrations')
	parser_test.add_argument('infile',
		type=str,
		help='Sample to test')
	parser_test.add_argument('outfile',
		type=str,
		help='Basename of files to write')
	parser_test.add_argument('reference',
		type=str,
		help='Reference as previously created')
	parser_test.add_argument('-minzscore',
		type=float, default=None,
		help='Minimum absolute z-score')
	parser_test.add_argument('-chromosomes', 
		help="Integer of every chromosome to plot, comma delimited",
		type=(lambda x: map(int,x.split(','))),
		default=range(1,23) ),
	parser_test.add_argument('-mineffectsize',
		type=float, default=0,
		help='Minimum absolute relative change in read depth')
	parser_test.add_argument('-multitest',
		type=float, default=1000,
		help='Increase z-score to compensate for multiple sample testing')
	parser_test.add_argument('-minrefbins',
		type=int, default=25,
		help='Minimum amount of sensible ref bins per target bin')
	parser_test.add_argument('-repeats',
		type=int, default=5,
		help='Repeats when calling')
	parser_test.set_defaults(func=toolTest)

	# Result plotting
	parser_plot = subparsers.add_parser('plot',
		description='Plot results produced by sample testing')
	parser_plot.add_argument('infile',
		type=str,
		help='Tested sample to plot')
	parser_plot.add_argument('outfile',
		type=str,
		help='Basename of files to write')
	parser_plot.add_argument('-cytofile',
		type=str, default=None,
		help='File containing cytoband information')
		# Obtained here: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
	parser_plot.add_argument('-chromosomes', 
		help="Integer of every chromosome to plot, comma delimited",
		type=(lambda x: map(int,x.split(','))),
		default=range(1,23))
	parser_plot.add_argument('-columns',
		type=int, default = 2,
		help='Amount of columns to distribute plots over')
	parser_plot.add_argument('-filetype',
		type=str, default='pdf',
		help='Export image as filetype')
	parser_plot.add_argument('-size',
		type=float, nargs=2, default = [11.7, 8.3],
		help='Size of plot in inches, default is A4 landscape')
	parser_plot.add_argument('-mineffect',
		type=float, default = 1.5,
		help='Minimal percentual change in read depth of a call to mark it')
	parser_plot.set_defaults(func=toolPlot)

	# Report results
	parser_report = subparsers.add_parser('report',
		description='Report results produced by sample testing')
	parser_report.add_argument('testfile',
		type=str,
		help='Input of test')
	parser_report.add_argument('resultfile',
		type=str,
		help='Output of test')
	parser_report.add_argument('-mineffect',
		type=float, default = 1.5,
		help='Minimal percentual change in read depth of a call to report')
	parser_report.set_defaults(func=toolReport)

	args = parser.parse_args(sys.argv[1:])
	printArgs(args)
	args.func(args)


if __name__ == '__main__':
	main()
