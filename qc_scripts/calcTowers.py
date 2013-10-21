import sys
import pickle

val = 0
curFile = pickle.load(open(sys.argv[1],'rb'))
for chrom in curFile.keys():
	val += sum(curFile[chrom])
	
print val
