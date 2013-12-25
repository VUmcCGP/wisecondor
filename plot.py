##############################################################################
#                                                                            #
#    Plot test results of WISECONDOR.                                        #
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

def plotResults(sample, markedBins, kept, kept2, outputFile, zScoresDict,zSmoothDict,blindsDict,wastedBins,sampleName):
	import matplotlib.pyplot as plt
	from matplotlib.collections import BrokenBarHCollection

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
	ax.text(0.5, 0.93, 'Z-score versus chromosomal position - Sample ' + sampleName, ha='center', va='bottom')

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

	bLGreen	= plt.Rectangle((0, 0), 1, 1, fc=colorMarkedBin)
	bDGreen	= plt.Rectangle((0, 0), 1, 1, fc=colorMarkedSure)
	bPink	= plt.Rectangle((0, 0), 1, 1, fc=colorMarkedSmooth)
	bPurple	= plt.Rectangle((0, 0), 1, 1, fc=colorMarkedSureSm)
	bLGrey	= plt.Rectangle((0, 0), 1, 1, fc=colorBlinds)
	bGrey	= plt.Rectangle((0, 0), 1, 1, fc=colorWaste)
	
	ax.legend((lRed, lBlue, bLGreen, bDGreen, bPink, bPurple, bLGrey, bGrey), ('Z-score windowed method', 'Z-score individual bin method', 'Detected deviation by individual bin method', 'Called by individual bin method', 'Detected deviation by windowed method', 'Called by windowed method', 'Not enough reference bins available', 'Unmappable region'), 'lower left',prop={'size':6}, mode='expand', ncol=4)
	
	plt.savefig(outputFile+'.pdf',figsize=(11.7, 8.3),dpi=160)


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Plot results generated by WISECONDOR',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('plotdata', type=str,
					   help='datafile as output by test.py')
	parser.add_argument('outfile', type=str,
					   help='output file to store plot in, .pdf is added as extension')
					   
	parser.add_argument('-mpluse', default='agg', type=str, 
					help='make matplotlib use another backend for plotting')

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
	print 'Loading:\tSample:\t' + args.plotdata
	sampleData = pickle.load(open(args.plotdata,'rb'))
	
	sampleName=args.plotdata.split("/")[-1].split(".")[0]
			
	plotResults( \
		sampleData['sample'], \
		sampleData['markedBins'], \
		sampleData['kept'], \
		sampleData['kept2'], \
		args.outfile, \
		sampleData['zScoresDict'], \
		sampleData['zSmoothDict'], \
		sampleData['blindsDict'], \
		sampleData['wastedBins'], \
		sampleName
		)
	print '\n# Finished'
