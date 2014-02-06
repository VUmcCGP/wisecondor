##############################################################################
#																			#
#	Test a maternal DNA sample for fetal Copy Number Aberrations.		   #
#	Copyright(C) 2013  TU Delft & VU University Medical Center Amsterdam	#
#	Author: Roy Straver, r.straver@vumc.nl								  #
#																			#
#	This file is part of WISECONDOR.										#
#																			#
#	WISECONDOR is free software: you can redistribute it and/or modify	  #
#	it under the terms of the GNU General Public License as published by	#
#	the Free Software Foundation, either version 3 of the License, or	   #
#	(at your option) any later version.									 #
#																			#
#	WISECONDOR is distributed in the hope that it will be useful,		   #
#	but WITHOUT ANY WARRANTY; without even the implied warranty of		  #
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		   #
#	GNU General Public License for more details.							#
#																			#
#	You should have received a copy of the GNU General Public License	   #
#	along with WISECONDOR.  If not, see <http://www.gnu.org/licenses/>.	 #
#																			#
##############################################################################


import cutoff
import gcc

import argparse
import math
import logging
import pickle
import pprint
import sys

import matplotlib
import numpy

numpy.seterr('ignore')

class WiseCondorTest(object):

	def __init__(self, samplefile, referencefile, binsize=1000*1000, verbose=True, *args, **kwargs):

		self.samplefile = samplefile
		self.sample = self.readpickle( samplefile )
		self.reference = self.readpickle( referencefile )

		self.refmaxbin = kwargs.get('refmaxbin')
		self.refminbin = kwargs.get('refminbin')
		self.maxrounds = kwargs.get('maxrounds')
		self.trithres = kwargs.get('trithres')
		self.window = kwargs.get('window')

		self.STOUFFER_CUTOFF = 3
		self.BINSIZE = binsize
		self.verbose = verbose
		
		# variables to fill during the run:
		
		self.markedbins = []
		self.zscores = {}
		self.markedsmoothbins = []
		self.zsmooth = {}
		self.blinds = {}
		

	def readpickle(self, picklefile):
		fd = self.readfile(picklefile)
		contents = pickle.load( fd )
		fd.close()
		return contents

	def readfile(self, filename, mode='rb'):
		self.file = open(filename, mode)
		return self.file

	def getReference(self, chrom, bin, markedBins, maxDist):
		"""
			Get reference values for a sample, only return values that are not marked
		"""
		reference = []

		if len(self.reference[chrom]) <= bin:
			logger.debug('ERROR: Unexpected end of bin list: ' + chrom + ':' + str(bin))
			return reference

		for value in self.reference[chrom][bin]:
			if [int(value[0]),value[1]] in [marked[:2] for marked in markedBins]:
				"""
					Skip processing is a certain bin combination is already in the markedBins
				"""
				logger.debug('Exclusion found!')
				continue

			if len(self.sample[value[0]]) > value[1]:
				# Only add bin if the distance is small enough
				if value[2] <= maxDist:
					reference.append(self.sample[value[0]][value[1]])
				else:
					break # Stop trying, only worse bins to come

		# Ignore bin because of too few reference bins
		if len(reference) < self.refminbin:
			return []
		# Ignore bins after maxBins is reached
		return reference[:self.refmaxbin]

	def checkAverageDeviation(self, maxDist):
		sample = self.sample
		
		deviations = []
		for chrom in range(1,23):
			for bin in range(0, len(self.sample[str(chrom)])-1 ):
				reference = self.getReference(str(chrom), bin, [], maxDist)
				dev = numpy.std(reference)/numpy.median(reference)
				if not math.isnan(dev):
					deviations.append(dev)

		avgDev = numpy.average(deviations)
		if self.verbose:
			devString = 'Average allowed deviation: ' + str(avgDev*100*3)[:5] + '%'
			if avgDev*3 > 0.05: # Unlikely to call anything sensible when over 5% deviation is considered normal
				devString += '\tWARNING: High value (>5%) calls are unreliable'
			logger.warn(devString)
		return avgDev*3


	def getZScore(self, value, reference):
		average	= numpy.average(reference)
		stddev	= numpy.std(reference)

		if stddev == 0:
			return 0
		Z = (value - average) / stddev
		return Z


	def markBins(self, maxDist, smoothRange=None):
		totalBins = sum([len(self.sample[str(chrom)]) for chrom in range(1,23)])
		prevMarks = [[0,0,0]]
		markedBins = []
		rounds = 1
		zScoresDict = dict()
		zSmoothDict = dict()
		blindsDict = dict()
		
		if not smoothRange:
			smoothRange = self.window
		
	
		while ([marked[:2] for marked in prevMarks] != [marked[:2] for marked in markedBins]) and rounds <= self.maxrounds:
			logger.info('\tRound: ' + str(rounds) + '\tMarks: ' + str(len(markedBins)))
			
			proBins = 0
			rounds += 1
			prevMarks = markedBins
			markedBins = []
			storeLen = 0

			for chrom in range(1,23):
				startBin = 0
				endBin = 0
				zScores = []
				blinds = []
				chromFound = 0
				avgBins = 0
				average = 0

				for bin in range(0,len(self.sample[str(chrom)])-1):
					value = self.sample[str(chrom)][bin]
					reference = self.getReference(str(chrom),bin,prevMarks,maxDist)
					zValue = self.getZScore(value, reference)

					if abs(zValue) > 0:
						avgBins += 1
						average += zValue

					if reference == []:
						blinds.append(bin)

					if (abs(zValue) >= 3):# and not ([chrom,bin] in [marked[:2] for marked in markedBins]):
						markedBins.append([chrom,bin,zValue])

					zScores.append(zValue)

				zScoresDict[str(chrom)] = zScores
				blindsDict[str(chrom)] = blinds

		if self.verbose:
			print 'Stopped\tMarks: ' + str(len(markedBins))


		# This part is for the windowed and chromosome wide scoring
		markedSmoothBins = []
		for chrom in zScoresDict:
			zSmooth = [1] * len(zScoresDict[str(chrom)])
			for bin in range(len(zScoresDict[str(chrom)])):
				"""
					Slice a window out of the zScoresDict to apply smooting on
				"""
				slice_start = max(0,bin-smoothRange)
				slice_end = min(bin+smoothRange+1,len(zScoresDict[str(chrom)]))
				temp = zScoresDict[str(chrom)][slice_start:slice_end]
				temp = numpy.array(temp)
				# Remove NaN values from the set
				temp = temp[~numpy.isnan(temp)]
				temp.sort()
				temp = temp[1:-1] # why strip the first and last item?
				zSmooth[bin] = numpy.sum(temp) / numpy.sqrt(len(temp))
				zSmoothDict[chrom] = zSmooth

			for bin in range(0,len(zSmooth)):

				if (abs(zSmooth[bin]) >= 3):
					markedSmoothBins.append([chrom,bin,zSmooth[bin]])

		markedBins.sort()
		
		self.markedbins = markedBins
		self.zscores = zScoresDict
		self.markedsmoothbins = markedSmoothBins
		self.zsmooth = zScoresDict
		self.blinds = blindsDict
		
		return markedBins, zScoresDict, markedSmoothBins, zSmoothDict, blindsDict

	def getMulti(self,chrom,start,end):
		totalVals = []
		for mark in range(start,end+1):
			binVals = []
			if len(self.reference[str(chrom)]) > mark:
				for lookUp in self.reference[str(chrom)][mark]:
					if len(self.sample[lookUp[0]]) > lookUp[1] and not self.sample[lookUp[0]][lookUp[1]] == 0:
						binVals.append(self.sample[lookUp[0]][lookUp[1]])
				if len(self.reference[str(chrom)][mark]) > 0:
					binAvg = self.sample[str(chrom)][mark]/numpy.average(binVals)
					totalVals.append(binAvg)
		totalAvg = numpy.average(totalVals)
		return totalAvg


	def testBins(self, markedBins,maxBinSkip,minBinLength):
		found = 0
		results = []
		count = 0
		prevChrom = 1
		curStrand = [[-1,-1,-1]]

		def finalize(curStrand):
			results.append(curStrand)

		for bin in markedBins:
			emptyBins = sum([val == [] for val in self.reference[str(bin[0])][bin[1]:curStrand[-1][1]]])
			distance = (bin[1] - curStrand[-1][1]) - emptyBins
			if (bin[0] == curStrand[0][0]) and (distance <= maxBinSkip) and (bin[2] * curStrand[0][2] > 0):
				curStrand.append(bin)
			else:
				finalize(curStrand)
				curStrand = []
				curStrand.append(bin)
		finalize(curStrand)

		results.pop(0)
		filtered = 0
		kept = []
	
		if self.verbose:
			print "\tChange\tMulti\tZ-Score\tSize\tPosition"
		for result in results:
			if len(result) >= minBinLength:
				found += 1
				chrom = result[0][0]
				startBin = result[0][1]
				endBin = result[-1][1]
				kept.append([chrom,startBin,endBin])

				zAverage = numpy.average([value[2] for value in result])
				variation = '-'
				if zAverage > 0:
					variation = '+'

				if self.verbose:
					print '\t' + variation + '\t' + str(self.getMulti(chrom,startBin,endBin))[:5] + \
						'\t' + str(zAverage)[:5] + '\t' + str((endBin-startBin + 1)*binSize/1000000) + \
						'Mb\tchr' + str(chrom) + ':' + str(startBin*binSize) + '-' + str((endBin+1)*binSize)
			else:
				filtered += 1
		if self.verbose:
			print 'Found:\t' + str(found) + '\tFiltered:\t' + str(filtered)

		return kept

	def testTrisomyAlt(self,kept):
	
		resultList = []
		for chrom in range(1,23):
			tempKept = []
			for val in kept:
				if val[0] == str(chrom):
					tempKept.append(val)

			marked = 0
			for val in tempKept:
				marked += val[2]-val[1]+1

				for blind in self.blinds[str(chrom)]:
					if blind >= val[1] and blind <= val[2]:
						marked -= 1

			# Skip testing if there are no testable bins on the chromosome
			if len(self.zsmooth[str(chrom)]) - len(self.blinds[str(chrom)]) == 0:
				continue

			result = marked / float(len(self.zsmooth[str(chrom)]) - len(self.blinds[str(chrom)]))
			if result > self.trithres:
				if self.verbose:
					print '\tchr' + str(chrom) + ': ' + str(round(result*100)) + '% marked'
				resultList.append([str(chrom), round(result*100)])

		return resultList



	def TrisomyStoufferDirect(self, zScoresDict):
		"""
			Test for Trisomy based on Stouffers Direct score
		"""
		self.dStouffDirect = {}

		for chromosome in zScoresDict.keys():
			zScoresDict[chromosome] = numpy.array(zScoresDict[chromosome])
			# First remove the NaN scores from the zScoresDict
			zscores = zScoresDict[chromosome][~numpy.isnan(zScoresDict[chromosome])]
			zscores.sort()

			behead = int( 0.05 * len(zscores) )
			behead = behead and 1 or behead
			selected_zscores = zscores[behead:-behead]
			self.dStouffDirect[ chromosome ] = (numpy.sum(selected_zscores)/numpy.sqrt(len(selected_zscores)))

			if self.verbose and abs( self.dStouffDirect[chromosome] ) > self.STOUFFER_CUTOFF:
				print "\tChr{}\t{:0.12f}".format( chromosome, self.dStouffDirect[chromosome] )
		return self.dStouffDirect




# --- MAIN ---
if __name__ == "__main__":
	logger = logging.getLogger('WiseCondorTest')
	logger.setLevel(logging.DEBUG)

	parser = argparse.ArgumentParser(description='WISECONDOR \
		(WIthin-SamplE COpy Number aberration DetectOR): \
		Detect fetal trisomies and smaller CNV\'s in a maternal plasma sample',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('sample', type=str,
					   help='sample to be tested (pickled)')

	parser.add_argument('reference', type=str,
					   help='reference table used for within sample comparison')
	parser.add_argument('outfile', type=str,
					   help='output pickled results to this file for plotting later on')

	parser.add_argument('-maxrounds', default=5, type=int,
					   help='maximum amount of rounds used to test for aberrations')
	parser.add_argument('-binsize', default=1000000, type=int,
					   help='binsize used for samples (equals arg used in reference)')
	parser.add_argument('-refminbin', default=10, type=int,
					   help='minimum number of reference bins, ignore target bin if not there are less reference bins available')
	parser.add_argument('-refmaxbin', default=100, type=int,
					   help='maximum number of reference bins, ignore any reference bin after')

	parser.add_argument('-refmaxval', default=1000000, type=int,
					   help='start cutoff value for determining good quality reference bins')
	parser.add_argument('-refmaxrep', default=3, type=int,
					   help='amount of improval rounds for determining good quality reference bins, use 0 to ignore this and just take -refmaxval as cutoff')

	parser.add_argument('-smaxskip', default=2, type=int,
					   help='single-bin-testing: ignore gaps of specified amount of bins')
	parser.add_argument('-sminbins', default=10, type=int,
					   help='single-bin-testing: ignore aberrations smaller than specified amount of bins')

	parser.add_argument('-window', default=5, type=int,
					   help='window size for sliding window approach, number of bins is considered in each direction (i.e. using 3 results in using 3+1+3=7 bins per call)')
	parser.add_argument('-wmaxskip', default=1, type=int,
					   help='window-testing: ignore gaps of specified amount of bins')
	parser.add_argument('-wminbins', default=10, type=int,
					   help='window-testing: ignore aberrations smaller than specified amount of bins')

	parser.add_argument('-trithres', default=0.5, type=float,
						help='threshold value for determining aneuploidy')

	parser.add_argument('-v', '--verbose', action="store_true",
						help='Enables verbose output to the stdout')

	args = parser.parse_args()

	if args.verbose:
		logger.info('# Script information:')
		logger.info('\n# Settings used:')
		argsDict = args.__dict__
		argsKeys = argsDict.keys()
		argsKeys.sort()
		for arg in argsKeys:
			logger.info('\t'.join([arg,str(argsDict[arg])]))

	# instantiate object
	options = {
		'samplefile': args.sample,
		'referencefile': args.reference,
		'binsize': args.binsize,
		'refminbin': args.refminbin,
		'refmaxbin': args.refmaxbin,
		'maxrounds': args.maxrounds,
		
		# Single bin testing
		'smaxskip': args.smaxskip,
		'sminbins': args.sminbins,
		
		# Window bin testing
		'window': args.window,
		'wmaxskip': args.wmaxskip,
		'wminbins': args.wminbins,

		# Tri threshold
		'trithres': args.trithres,

		'verbose': args.verbose,
	}
	
	wc = WiseCondorTest( **options )


	if args.verbose:
		print '\n# Processing:'
		print 'Loading:\tSample:\t' + args.sample
		print 'Loading:\tReference Table\t' + args.reference
	sample = wc.sample
	lookUpTable = wc.reference
	binSize 	= int(args.binsize)


	if args.verbose:
		print '\nDetermining reference cutoffs'
	maxDist = cutoff.getOptimalCutoff(lookUpTable, args.refmaxrep, args.refmaxval)

	if args.verbose:
		print '\tCutoff determined:\t' + str(maxDist)
		print
	avgDev = wc.checkAverageDeviation( maxDist )

	markedBins,zScoresDict,markedSmoothBins,zSmoothDict,blindsDict = wc.markBins(maxDist)

	if args.verbose:
		print '\nUncallable bins:\t' + str(sum([len(blindsDict[key]) for key in blindsDict.keys()]) \
			/float(sum([len(sample[key]) for key in sample.keys()]))*100)[:5] + '%'

	if args.verbose:
		print '\n\n# Results:'
		print '\nSingle bin, bin test:'
	kept = wc.testBins(markedBins,args.smaxskip,args.sminbins)

	if args.verbose:
		print '\nSingle bin, aneuploidy test:'
	singlebin_test = wc.testTrisomyAlt(kept)

	if args.verbose:
		if len(singlebin_test) == 0:
			print 'Nothing found'


	if args.verbose:
		print '\nWindowed, bin test:'
	kept2 = wc.testBins(markedSmoothBins,args.wmaxskip,args.wminbins)

	if args.verbose:
		print '\nWindowed, aneuploidy test:'
	windowed_test = wc.testTrisomyAlt(kept2)

	if args.verbose:
		if len(windowed_test) == 0:
			print 'Nothing found'
	
	if args.verbose:
		print '\nChromosome wide, aneuploidy test:'
	wc.TrisomyStoufferDirect(zScoresDict)
	

	if args.verbose:
		print '\n\n# Script information:\n'
		print '\nComputing additional data for plots'
	wastedBins = dict()
	for chrom in range(1,23):
		wastedBins[str(chrom)] = []
		for bin in range(len(sample[str(chrom)])-1):
			wastedBins[str(chrom)].append(len(wc.getReference(str(chrom),bin,[],1)) <= 3)

	if args.verbose:
		print '\nStoring data for creating plots'
	outputData=dict()
	outputData['sample']=sample
	outputData['runtime_arguments'] = args
	outputData['avgDev'] = avgDev


	outputData['markedBins']=markedBins

	outputData['kept']=kept
	outputData['singlebin_test'] = singlebin_test

	outputData['kept2']=kept2
	outputData['windowed_test'] = windowed_test

	#outputData['outputFile']=outputFile
	outputData['zScoresDict']=zScoresDict
	outputData['zSmoothDict']=zSmoothDict
	outputData['blindsDict']=blindsDict
	outputData['wastedBins']=wastedBins
#	pickle.dump(outputData,open(args.outfile,'wb'))

	if args.verbose:
		print '\n# Finished'
