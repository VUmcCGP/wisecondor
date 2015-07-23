##############################################################################
#                                                                            #
#    Correct read depths for GC-content using a LOWESS function.             #
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



from math import *
import sys
import argparse
import Bio.Statistics.lowess as biostat
import numpy as np

def correct(sample,gcCount,binSize,maxN=0.1,minRD=0.0001,fVal=0.1,iVal=3):
	allX = []
	allY = []

	chroms = sample.keys()

	for chrom in chroms:
		for bin in range(min(len(gcCount[chrom]),len(sample[chrom]))):
			if gcCount['N'+chrom][bin] < binSize * maxN and sample[chrom][bin] > binSize * minRD:
				allX.append(gcCount[chrom][bin])
				allY.append(sample[chrom][bin])

	allX = np.array(allX,np.float)
	allY = np.array(allY,np.float)
	lowessCurve = biostat.lowess(allX,allY,f=fVal, iter=iVal).tolist()
	
	correctedSample = dict()
	for chrom in chroms:
		correctedSample[chrom] = []
		for bin in range(min(len(gcCount[chrom]),len(sample[chrom]))):
			if gcCount['N'+chrom][bin] < binSize * maxN and sample[chrom][bin] > binSize * minRD:
				correctedValue = sample[chrom][bin]/lowessCurve.pop(0)
				correctedSample[chrom].append(correctedValue)
			else:
				correctedSample[chrom].append(0)

	return correctedSample

if __name__ == "__main__":
	import pickle
	import argparse
	parser = argparse.ArgumentParser(description='Correct a sample for GC-Content using a LOESS function',
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('infile', type=str,
						help='directory containing samples to be used as reference (pickle)')
	parser.add_argument('gccount', type=str,
						help='gc-counts file used for gc-correction (pickle)')
	parser.add_argument('outfile', type=str,
						help='gc-corrected sample data')

	parser.add_argument('-binsize', default=1000000, type=int,
						help='binsize used for samples')
	parser.add_argument('-maxn', default=0.1, type=float,
						help='maximum relative amount of unknown (n) bases in bin used for gc-correction (equals arg used in test)')
	parser.add_argument('-minrd', default=0.0001, type=float,
						help='minimum relative amount of reads in bin used for gc-correction (equals arg used in test)')
	parser.add_argument('-fval', default=0.1, type=float,
						help='width of data used in loess function used for gc-correction (equals arg used in test)')
	parser.add_argument('-ival', default=3, type=int,
						help='amount of fitting iterations in loess function used for gc-correction (equals arg used in test)')

	args = parser.parse_args()

	print '\n# Settings used:'
	argsDict = args.__dict__
	argsKeys = argsDict.keys()
	argsKeys.sort()
	for arg in argsKeys:
		print '\t'.join([arg,str(argsDict[arg])])

	print '\n# Processing:'
	print 'Loading:\tSample:\t' + args.infile
	sample = pickle.load(open(args.infile,'rb'))
	print 'Loading:\tGC-Count:\t' + args.gccount
	gcCount = pickle.load(open(args.gccount,'rb'))
	print 'Correcting'
	corrected = correct(sample,gcCount,args.binsize,args.maxn,args.minrd,args.fval,args.ival)
	print 'Writing to file'
	pickle.dump(corrected,open(args.outfile,'wb'))
	print '\n# Finished'
