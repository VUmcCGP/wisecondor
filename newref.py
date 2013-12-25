##############################################################################
#                                                                            #
#    Find optimal cutoff for 'good' bins for a certain reference set.        #
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



import glob
import sys
import pickle
import gcc
import argparse

parser = argparse.ArgumentParser(description='Create a new reference table from a set of reference samples, outputs table as pickle to a specified output file',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('refdir', type=str,
					help='directory containing samples to be used as reference (pickle)')
parser.add_argument('refout', type=str,
					help='reference table output, used for sample testing (pickle)')

parser.add_argument('-ignore', type=int, default=0,
					help='ignore x highest scoring control samples per bin distance calculation, 0 to keep all, use if you are unsure about what samples have aberrations and assume at most x have the same aberration, beware: results become less accurate')

args = parser.parse_args()

print '\n# Settings used:'
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
	print '\t'.join([arg,str(argsDict[arg])])

print '\n# Processing:'

# Calculate distances from each bin from chromosome i to every bin of chromosome j
def getDistanceTable(controls,iChrom,jChrom):
	chromosomeDistances = []
	# For each bin...
	iLen = max([len(controls[key][iChrom]) for key in controls.keys()])
	for iBin in range(0,iLen):
		binDistances = []
		# For every other bin...
		jLen = max([len(controls[key][jChrom]) for key in controls.keys()])
		for jBin in range(0,jLen):
			# Take the difference over all samples...
			distance = 0
			distances = []
			for control in controls:
				# Get iBin
				iVal = 0
				if iBin < len(controls[control][iChrom]):
					iVal = controls[control][iChrom][iBin]
				# Get jBin
				jVal = 0
				if jBin < len(controls[control][jChrom]):
					jVal = controls[control][jChrom][jBin]
				# Calculate the distance
				distance += pow((iVal - jVal),2)
				# Don't try to match with 0 bins later on, -1 if iBin or jBin is zerobin
				if (iVal == 0) or (jVal == 0):
					distance = -1
					distances = [-1]
					break
				distances.append(distance)
			# Append new found distance
			if args.ignore > 0 and len(distances) > 0 and distance != -1:
				distances.sort() # yeah
				binDistances.append(sum(distances[:-args.ignore]))
			else:
				binDistances.append(distance)
		# Append all distances of this bin to the complete set of distances per chromosome
		chromosomeDistances.append(binDistances)
	return chromosomeDistances

# Load the reference samples
print 'Loading reference samples'
controls	= dict()
referenceFiles = glob.glob(args.refdir + '/*.gcc')
for refFile in referenceFiles:
	print '\tLoading:\t' + refFile
	curFile = pickle.load(open(refFile,'rb'))
	controls[refFile] = curFile

print 'Building reference table'
chromList = [str(chrom) for chrom in range(1,23)]
refTable = dict()
for chrom in chromList:
	refTable[chrom] = []

for iChrom in chromList:
	print '\tTargeting chromosome:\t' + str(iChrom)
	jChroms = dict()
	for jChrom in chromList:
		print '\t\tCalculating distances to chromosome:\t' + str(jChrom)
		jChroms[jChrom] = getDistanceTable(controls,iChrom,jChrom)

	# Remove own chromosome from the list of referable chromosomes
	iChromBins = jChroms.pop(iChrom)
	print '\t\tPicking reference bins'
	for iBin in range(len(iChromBins)):
		# Tuple: [Value,Bin,Chromosome]
		topRanks = [[sys.maxint,-1,'']] * 250 # plenty of spots to go around...

		def updateBestBins(bins,chrom):
			# Return position of the worst bin in our list
			def getWorstPos():
				worstPos = 0
				for i in range(len(topRanks)):
					if topRanks[i][0] > topRanks[worstPos][0]:
						worstPos = i
				return worstPos

			worstPos = getWorstPos()
			for i in range(len(bins)):
				# Zero bins are marked by -1, ignore them
				if float(bins[i]) >= 0:
					if float(bins[i]) < topRanks[worstPos][0]:
						topRanks[worstPos] = [float(bins[i]),i,chrom]
						worstPos = getWorstPos()

		# Read and process line for each file we have of autosomal chromosomes
		for jChrom in jChroms.keys():
			updateBestBins(jChroms[jChrom][iBin],jChrom)

		# Don't take bins close to eachother, take the best one instead
		minDistance = 2
		remove = []
		lastChrom = 0
		lastBin = 0
		tempList = []
		for bin in topRanks:
			if (bin[2] == lastChrom) and (abs(bin[1]-lastBin) < minDistance):
				tempList.append(bin)
			else:
				bestPos = 0
				for temp in range(len(tempList)):
					if tempList[bestPos][0] > tempList[temp][0]:
						bestPos = temp
				# remove best bin from list to remove
				if tempList != []:
					tempList.pop(bestPos) 
				# add other bins to the remove list
				remove += tempList
				tempList = [bin]
			lastBin = bin[1]
			lastChrom = bin[2]

		for bin in remove:
			topRanks.remove(bin)

		#print '\t\tSorting reference bins'
		topRanks.sort()
		tempTop = []
		for jBin in topRanks:
			if jBin[1] >= 0:
				jBin.reverse()
				tempTop.append(jBin)
		refTable[iChrom].append(tempTop[:100])

print 'Writing reference to file'
pickle.dump(refTable, open(args.refout, 'wb'))
