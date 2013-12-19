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



import sys
import numpy
import argparse

def getReference(lookUp, cutOff):
	reference = []
	removed = 0
	for chrom in lookUp:
		for bin in lookUp[chrom]:
			if len(bin) > 0:
				if float(bin[0][2]) < cutOff:
					reference.append(float(bin[0][2]))
				else:
					removed += 1
	return reference,removed

def getOptimalCutoff(lookUp, repeats, optimalCutoff):
	for i in range(0,repeats):
		reference,removed = getReference(lookUp, optimalCutoff)
		average	= numpy.average(reference)
		stddev	= numpy.std(reference)
		optimalCutoff = average + 3 * stddev
	return optimalCutoff

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Determine optimal cutoff value for the reference table provided',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('reference', type=str,
		               help='reference table to work on (pickle)')
	parser.add_argument('-refmaxval', default=1000000, type=int,
		               help='start cutoff value for determining good quality reference bins')
	parser.add_argument('-refmaxrep', default=3, type=int,
		               help='amount of improval rounds for determining good quality reference bins')
	args = parser.parse_args()

	getOptimalCutoff(args.reference,args.refmaxval,args.refmaxrep)
