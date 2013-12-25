##############################################################################
#                                                                            #
#    Convert reference .fa to a pickled list with GC-content per bin.        #
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
import argparse

parser = argparse.ArgumentParser(description='Create GC-count file for GC-corrections. Outputs table as pickle to a specified output file',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('reffile', type=str,
					help='reference file used to map reads to (hg19 etc)')
parser.add_argument('gccout', type=str,
					help='gc-count output file')

parser.add_argument('-binsize', type=int, default=1000000,
					help='binsize used for samples')

args = parser.parse_args()

filename	= args.reffile
binSize = args.binsize

gcCounts = dict()
for chrom in range(1,23):
	gcCounts[str(chrom)] = []
gcCounts['X'] = []
gcCounts['Y'] = []

nCounts = dict()
for chrom in range(1,23):
	nCounts[str(chrom)] = []
nCounts['X'] = []
nCounts['Y'] = []


key = ''
start = 0
gcCount = 0.
nCount = 0.
totalCount = 0


def finishBin():
	global gcCount
	global totalCount
	global nCount
	global start
	global key
	if key in gcCounts:
		binNumber = (start+totalCount) / binSize - 1
		while len(gcCounts[key]) <= binNumber:
			gcCounts[key].append(0)
			nCounts[key].append(0)
		gcCounts[key][binNumber] = gcCount
		nCounts[key][binNumber] = nCount

		print "\tChr: "+str(key)+"\tBin: "+str(binNumber)+"\t-\tGCC: "+str(gcCount)+"\tNC: "+str(nCount)+"\tTotal: "+str(totalCount)#+"\tStart: "+str(start)
	gcCount = 0.
	nCount = 0.

f = open(filename, 'r')
for nextLine in f:
	if nextLine[0] == '>':
		finishBin()
		totalCount = 0
		splitLine = nextLine.split()
		key   = splitLine[0][1:]
		if key[:3] == 'chr':
			key = key[3:]
		start = 0
		print '\nWorking on:\t' + key
	elif key in gcCounts.keys():
		rstripLine = nextLine.rstrip()
		for char in rstripLine:
			if char in 'GgCc': # Beware the softmasked fasta files...
				gcCount += 1
			elif char in 'Nn':
				nCount += 1
			totalCount += 1
			if (start + totalCount) % binSize == 0:
				finishBin()
f.close()
finishBin()

# Now put it all together in a single dict for ease of use
for chrom in nCounts.keys():
	gcCounts['N'+chrom]=nCounts[chrom]
	
# And dump it to a file
pickle.dump(gcCounts, open(args.gccout, 'wb'))
