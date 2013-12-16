##############################################################################
#                                                                            #
#    Test a maternal DNA sample for fetal Copy Number Aberrations.           #
#    Copyright(C) 2013  TU Delft & VU University Medical Center Amsterdam    #
#    Author: Roy Straver, r.straver@vumc.nl                                  #
#                                                                            #
#    This file is part of WISECONDOR.                                        #
#                                                                            #
#    WISECONDOR is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    WISECONDOR is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with WISECONDOR.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                            #
##############################################################################



import sys
import pickle
import gcc

import numpy
import math
import cutoff

import matplotlib
import argparse

numpy.seterr('ignore')

binSize = 0
lookUpTableDir = None
lookUpTable = None

def getZScore(value, reference):
	average	= numpy.average(reference)
	stddev	= numpy.std(reference)

	if stddev == 0:
		return 0

	Z = (value - average) / stddev

	return Z

def getReference(sample,chrom,bin,markedBins,minBins,maxBins,maxDist):
	reference = []
	if len(lookUpTable[chrom]) <= bin:
		#print 'ERROR: Unexpected end of bin list: ' + chrom + ':' + str(bin)
		return reference
	for value in lookUpTable[chrom][bin]:
		if [int(value[0]),value[1]] in [marked[:2] for marked in markedBins]:
			#print 'Exclusion found!'
			continue
		if len(sample[value[0]]) > value[1]:
			# Only add bin if the distance is small enough
			if value[2] <= maxDist:
				reference.append(sample[value[0]][value[1]])
			else:
				break # Stop trying, only worse bins to come

	# Ignore bin because of too few reference bins
	if len(reference) < minBins:
		return []
	# Ignore bins after maxBins is reached
	return reference[:maxBins]

def checkAverageDeviation(sample,minBins,maxBins,maxDist):
	deviations = []
	for chrom in range(1,23):
		for bin in range(0,len(sample[str(chrom)])-1):
			reference = getReference(sample,str(chrom),bin,[],minBins,maxBins,maxDist)				
			dev = numpy.std(reference)/sample[str(chrom)][bin]
			if not math.isnan(dev):
				deviations.append(dev)

	avgDev = numpy.average(deviations)
	devString = 'Average allowed deviation: ' + str(avgDev*100*3)[:5] + '%'
	if avgDev*3 > 0.05: # Unlikely to call anything sensible when over 5% deviation is considered normal
		devString += '\tWARNING: High value (>5%) calls are unreliable'
	print devString
	return avgDev*3
	#quit()

def markBins(sample,maxRounds,minBins,maxBins,maxDist,smoothRange):
	#return 0

	totalBins = sum([len(sample[str(chrom)]) for chrom in range(1,23)])

	prevMarks = [[0,0,0]]
	markedBins = []
	rounds = 1
	zScoresDict = dict()
	zSmoothDict = dict()
	blindsDict = dict()
	
	while ([marked[:2] for marked in prevMarks] != [marked[:2] for marked in markedBins]) and rounds <= maxRounds:

		print '\tRound: ' + str(rounds) + '\tMarks: ' + str(len(markedBins))

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


			for bin in range(0,len(sample[str(chrom)])-1):
				value = sample[str(chrom)][bin]
				reference = getReference(sample,str(chrom),bin,prevMarks,minBins,maxBins,maxDist)
				zValue = getZScore(value, reference)

				if abs(zValue) > 0:
					avgBins += 1
					average += zValue

				if reference == []:
					blinds.append(bin)

				if (abs(zValue) >= 3):# and not ([chrom,bin] in [marked[:2] for marked in markedBins]):
					markedBins.append([chrom,bin,zValue])

				zScores.append(zValue)

				#mult = value/numpy.average(reference)
				#if abs(mult - 1) > 0.05 and mult != 0.0:
				#	print str(chrom) + '\t' + str(zValue)[:5] + '\t' + str((mult-1)*100)[:5]
					#markedBins.append([chrom,bin,zValue])

			zScoresDict[str(chrom)] = zScores
			blindsDict[str(chrom)] = blinds

			#print str(chrom) + '\tAverage Deviation: ' + str(numpy.average(binDevs)*100)[:5] + '%'


		#print markedBins
	print 'Stopped\tMarks: ' + str(len(markedBins))

	markedSmoothBins = []
	for chrom in zScoresDict:
		zSmooth = [1] * len(zScoresDict[str(chrom)])
		#smoothRange = 5
		for bin in range(len(zScoresDict[str(chrom)])):
			temp = zScoresDict[str(chrom)][max(0,bin-smoothRange):min(bin+smoothRange+1,len(zScoresDict[str(chrom)]))]
			# Remove NaN values from the set
			temp = [val for val in temp if not math.isnan(val)]
			# Get lost, frigging peak.
			temp.sort()
			temp = temp[1:-1]
			zSmooth[bin] = numpy.sum(temp)/numpy.sqrt(len(temp))
			zSmoothDict[chrom] = zSmooth

		for bin in range(0,len(zSmooth)):

			if (abs(zSmooth[bin]) >= 3):# and not ([chrom,bin] in [marked[:2] for marked in markedBins]):
				markedSmoothBins.append([chrom,bin,zSmooth[bin]])

	markedBins.sort()
	return markedBins, zScoresDict, markedSmoothBins, zSmoothDict, blindsDict

# TODO: Only give multi for enough reference bins, otherwise nan?
def getMulti(sample,chrom,start,end):
	totalVals = []
	for mark in range(start,end+1):
		binVals = []
		if len(lookUpTable[str(chrom)]) > mark:
			for lookUp in lookUpTable[str(chrom)][mark]:
				if len(sample[lookUp[0]]) > lookUp[1] and not sample[lookUp[0]][lookUp[1]] == 0:
					binVals.append(sample[lookUp[0]][lookUp[1]])
			if len(lookUpTable[str(chrom)][mark]) > 0:
				binAvg = sample[str(chrom)][mark]/numpy.average(binVals)
				totalVals.append(binAvg)
	totalAvg = numpy.average(totalVals)
	return totalAvg


def testBins(sample,markedBins,maxBinSkip,minBinLength):
	found = 0
	results = []

	count = 0
	prevChrom = 1
	curStrand = [[-1,-1,-1]]

	def finalize(curStrand):
		results.append(curStrand)

	for bin in markedBins:
		emptyBins = sum([val == [] for val in lookUpTable[str(bin[0])][bin[1]:curStrand[-1][1]]])
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
			
			print '\t' + variation + '\t' + str(getMulti(sample,chrom,startBin,endBin))[:5] + \
					'\t' + str(zAverage)[:5] + '\t' + str((endBin-startBin + 1)*binSize/1000000) + \
					'Mb\tchr' + str(chrom) + ':' + str(startBin*binSize) + '-' + str((endBin+1)*binSize)
		else:
			filtered += 1
	print 'Found:\t' + str(found) + '\tFiltered:\t' + str(filtered)

	return kept

def plotResults(sample, markedBins, kept, kept2, outputFile, zScoresDict,zSmoothDict,blindsDict):
	import matplotlib.pyplot as plt
	from matplotlib.collections import BrokenBarHCollection

	#(1,0,0.6)#=geenstijlroze
	#(0.4,0.76,0.13)#=dumpertgroen

	colorSample			= 'blue'
	colorReference 		= 'red'
	colorMarkedBin		= (0.5,1,0.35)
	colorMarkedSure		= (0.3,0.6,0.07)#'green'
	colorMarkedSmooth	= (1,0.6,1)
	colorMarkedSureSm	= (0.4,0.2,0.4)
	colorHorzHelper		= (0.7,0.7,0.7)
	colorHorzMarker		= 'orange'
	colorBlinds			= (0.85,0.85,0.85)
	colorWaste			= (0.7,0.7,0.7)

	binScalar = 320
	edgeSize = 0.15

	wastedBins = dict()
	for chrom in range(1,23):
		wastedBins[str(chrom)] = []
		for bin in range(len(sample[str(chrom)])-1):
			wastedBins[str(chrom)].append(len(getReference(sample,str(chrom),bin,[],0,4,1)) <= 3)

	def ideogram(chrom):

		color_lookup = {
			'gneg': (1., 1., 1.),
			'gpos25': (.6, .6, .6),
			'gpos50': (.4, .4, .4),
			'gpos75': (.2, .2, .2),
			'gpos100': (0., 0., 0.),
			'acen': (.8, .4, .4),
			'gvar': (.8, .8, .8),
			'stalk': (.9, .9, .9),
			'select': (0., 0., 1.)
		}

		fin = open('./ideogram.txt')
		prevChrom = ''
		for line in fin:
			lchrom, start, stop, label, stain = line.strip().split('\t')
			if ('chr'+ str(chrom)) == lchrom:
				start = int(start)/float(binSize)
				stop = int(stop)/float(binSize)
				
				width = stop - start
				center = start + width/2.
				plt.axvline(x=center,  ymin=1-edgeSize, ymax=1, linewidth=width, color=color_lookup[stain])
			elif prevChrom != lchrom and prevChrom == chrom:
				break
			prevChrom = lchrom
		fin.close()

		frame1 = plt.gca()
		frame1.axes.get_xaxis().set_visible(False)
		frame1.axes.get_yaxis().set_visible(False)

	def printMarkers(chrom):
		binWidth = 1/float(len(sample[str(chrom)]))*binScalar
		
		zSmooth = zSmoothDict[str(chrom)]
		for bin in range(len(zSmooth)):
			if abs(zSmooth[bin]) > 3:
				plt.axvline(x=bin, ymin=1-edgeSize*2, ymax=1, linewidth=binWidth, color=colorMarkedSmooth)

		for region in kept2:
			if str(region[0]) == chrom:
				for bin in range(region[1],region[2]+1):
					plt.axvline(x=bin, ymin=1-edgeSize, ymax=1, linewidth=binWidth, color=colorMarkedSureSm)

		marks = [1] * len(sample[str(chrom)])
		for bin in markedBins:
			if str(bin[0]) == chrom:
				plt.axvline(x=bin[1], ymin=0, ymax=edgeSize*2, linewidth=binWidth, color=colorMarkedBin)

		for region in kept:
			if str(region[0]) == chrom:
				for bin in range(region[1],region[2]+1):
					plt.axvline(x=bin, ymin=0, ymax=edgeSize, linewidth=binWidth, color=colorMarkedSure)

	def drawBlinds(chrom):
		binWidth = 1/float(len(sample[str(chrom)]))*binScalar
		for bin in blindsDict[str(chrom)]:
			plt.axvline(x=bin, linewidth=binWidth, color=colorBlinds)
		for bin in range(len(wastedBins[str(chrom)])):
			if wastedBins[str(chrom)][bin]:
				plt.axvline(x=bin, linewidth=binWidth, color=colorWaste)

	def preparePlot(chrom):
		plt.subplot(11,2,chrom)
		drawBlinds(chrom)
		printMarkers(str(chrom))
		
	print 'Plotting Z-Scores'
	ax = plt.figure(2)
	ax.text(0.5, 0.06, 'Chromosomal position in bins', ha='center', va='bottom')
	ax.text(0.05, 0.5, 'Z-score', ha='left', va='center', rotation='vertical')
	ax.text(0.5, 0.93, 'Z-score versus chromosomal position - Sample ' + args.sample.split("/")[-1].split(".")[0], ha='center', va='bottom')

	for chrom in range(1,23):
		preparePlot(chrom)

		zSmooth = zSmoothDict[str(chrom)]
		
		zKeep = []
		for val in zScoresDict[str(chrom)]:
			if not math.isnan(val):
				zKeep.append(val)
		
		zTotal = numpy.sum(zKeep) / numpy.sqrt(len(zKeep))

		plt.axhline(y=0, linewidth=0.75, color=colorHorzHelper)
		plt.axhline(y=3, linewidth=0.50, color=colorHorzHelper)
		plt.axhline(y=-3, linewidth=0.50, color=colorHorzHelper)

		lRed, = plt.plot(zSmooth,colorReference)
		lBlue, = plt.plot(zScoresDict[str(chrom)],colorSample)

		plt.ylabel(chrom)
		plt.xlim(0,len(sample[str(chrom)])-2)
		frame1 = plt.gca()
		#frame1.axes.get_xaxis().set_visible(False)
		frame1.axes.get_xaxis().get_major_ticks()
		frame1.axes.xaxis.set_ticklabels([])

		for tick in frame1.axes.get_yaxis().get_major_ticks():
			tick.label.set_fontsize(4) 

	bLGreen = plt.Rectangle((0, 0), 1, 1, fc=colorMarkedBin)
	bDGreen = plt.Rectangle((0, 0), 1, 1, fc=colorMarkedSure)
	bPink = plt.Rectangle((0, 0), 1, 1, fc=colorMarkedSmooth)
	bPurple = plt.Rectangle((0, 0), 1, 1, fc=colorMarkedSureSm)
	bLGrey = plt.Rectangle((0, 0), 1, 1, fc=colorBlinds)
	bGrey = plt.Rectangle((0, 0), 1, 1, fc=colorWaste)
	
	ax.legend((lRed, lBlue, bLGreen, bDGreen, bPink, bPurple, bLGrey, bGrey), ('Z-score windowed method', 'Z-score individual bin method', 'Detected deviation by individual bin method', 'Called by individual bin method', 'Detected deviation by windowed method', 'Called by windowed method', 'Not enough reference bins available', 'Unmappable region'), 'lower left',prop={'size':6}, mode='expand', ncol=4)
	
	plt.savefig(outputFile+'.zscore.pdf',figsize=(11.7, 8.3),dpi=160)


def testTrisomyAlt(sample,kept,zValues,blindsDict,threshold):
	resultList = []
	for chrom in range(1,23):
		tempKept = []
		for val in kept:
			if val[0] == str(chrom):
				tempKept.append(val)

		marked = 0
		#print tempKept
		for val in tempKept:
			marked += val[2]-val[1]+1

			for blind in blindsDict[str(chrom)]:
				if blind >= val[1] and blind <= val[2]:
					marked -= 1

		# Skip testing if there are no testable bins on the chromosome
		if len(zValues[str(chrom)]) - len(blindsDict[str(chrom)]) == 0:
			continue

		result = marked / float(len(zValues[str(chrom)]) - len(blindsDict[str(chrom)]))
		if result > threshold:
			print '\tchr' + str(chrom) + ': ' + str(round(result*100)) + '% marked'
			resultList.append([str(chrom), round(result*100)])

	return resultList

def testTrisomyStoufferDirect(zScoresDict):
	stouff = []
	for chromInt in range(22):
		chrom = str(chromInt+1)
		temp = [val for val in zScoresDict[chrom] if not math.isnan(val)]
		temp.sort()
		behead = int(0.05 * len(temp))
		if behead == 0:
			behead = 1
		temp = temp[behead:-behead]
			#print str(chrom) + "\t" + str(zScore)
		stouff.append((numpy.sum(temp)/numpy.sqrt(len(temp))))
	for chromInt in range(22):
		if abs(stouff[chromInt]) > 3:
			print "\tChr" + str(chromInt+1) + "\t" + str(stouff[chromInt])
	#print "chromWideStouff:\t" + "\t".join(stouff)

def loadSettings(settingsFile):
	settings = dict()
	f = open(settingsFile, 'r')
	for nextLine in f:
		if len(nextLine.split()) > 0 and nextLine[0] != '#':
			#print nextLine
			value = 0
			splitLine = nextLine.rsplit()
			try:	value = int(splitLine[1])
			except ValueError:
				try:	value = float(splitLine[1])
				except ValueError:
					value = splitLine[1]
			settings[splitLine[0]] = value
	f.close()
	#print settings
	return settings



import argparse
parser = argparse.ArgumentParser(description='WISECONDOR \
	(WIthin-SamplE COpy Number aberration DetectOR): \
	Detect fetal trisomies and smaller CNV\'s in a maternal plasma sample.',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('sample', type=str,
                   help='sample to be tested (pickle)')

parser.add_argument('reference', type=str,
                   help='reference table used for within sample comparison')
parser.add_argument('output', type=str,
                   help='output basename for plots, specify both path and file prefix (pdf)')

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
parser.add_argument('-mpluse', default='agg', type=str, 
					help='make matplotlib use another backend for plotting')

parser.add_argument('-notriper', action='store_true',
                   help='do not test for aneuploidy using the single bin approach')
parser.add_argument('-notriwin', action='store_true',
                   help='do not test for aneuploidy using the window bin approach')
parser.add_argument('-notrichr', action='store_true',
                   help='do not test for aneuploidy using the chromome wide z-score approach')
parser.add_argument('-noplot', action='store_true',
                   help='do not plot results in a pdf, dependent on -notriper and -notriwin, if either is toggled off, no plots will be shown')


args = parser.parse_args()


print '# Script information:'

print '\n# Settings used:'

matplotlib.use(args.mpluse)
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
	print '\t'.join([arg,str(argsDict[arg])])

print '\n# Processing:'

print 'Loading:\tSample:\t' + args.sample
sample 		= pickle.load(open(args.sample,'rb'))
print 'Loading:\tReference Table\t' + args.reference
lookUpTable = pickle.load(open(args.reference,'rb'))
binSize 	= int(args.binsize)
outputBase	= args.output

print '\nDetermining reference cutoffs'
maxDist = cutoff.getOptimalCutoff(lookUpTable, args.refmaxrep, args.refmaxval)
print '\tCutoff determined:\t' + str(maxDist)


print ''
avgDev = checkAverageDeviation(sample,args.refminbin,args.refmaxbin,maxDist)
markedBins,zScoresDict,markedSmoothBins,zSmoothDict,blindsDict = \
	markBins(sample,args.maxrounds,args.refminbin,args.refmaxbin,maxDist,args.window)

print '\nUncallable bins:\t' + str(sum([len(blindsDict[key]) for key in blindsDict.keys()]) \
		/float(sum([len(sample[key]) for key in sample.keys()]))*100)[:5] + '%'
print '\n\n# Results:'

if not args.notriper:
	print '\nSingle bin, bin test:'
	kept = testBins(sample,markedBins,args.smaxskip,args.sminbins)

	print '\nSingle bin, aneuploidy test:'
	if len(testTrisomyAlt(sample,kept,zSmoothDict,blindsDict,args.trithres)) == 0:
		print 'Nothing found'

if not args.notriwin:
	print '\nWindowed, bin test:'
	kept2 = testBins(sample,markedSmoothBins,args.wmaxskip,args.wminbins)

	print '\nWindowed, aneuploidy test:'
	if len(testTrisomyAlt(sample,kept2,zSmoothDict,blindsDict,args.trithres)) == 0:
		print 'Nothing found'

if not args.notrichr:
	print '\nChromosome wide, aneuploidy test:'
	testTrisomyStoufferDirect(zScoresDict)

if not args.notriper and not args.notriwin and not args.noplot:
	print '\n\n# Script information:\n'
	plotResults(sample,markedBins,kept,kept2,outputBase,zScoresDict,zSmoothDict,blindsDict)

