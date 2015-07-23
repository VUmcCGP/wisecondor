import sys
import json
import pickle
import numpy

sampleData = pickle.load(open(sys.argv[1],'rb'))
#print sampleData.keys()

for chrom in sampleData["zScoresDict"]:
	for i,val in enumerate(sampleData["zScoresDict"][chrom]):
		if numpy.isnan(val):
			sampleData["zScoresDict"][chrom][i] = 0
for chrom in sampleData["zSmoothDict"]:
	for i,val in enumerate(sampleData["zSmoothDict"][chrom]):
		if numpy.isnan(val):
			sampleData["zSmoothDict"][chrom][i] = 0		



print json.dumps(sampleData)
