###############################################################################
#                                                                             #
#    DEFRAG	(DEtection of fetal FRaction And Gender)        				  #
#    Copyright(C) 2014  VU University Medical Center Amsterdam    			  #
#    Authors: 																  #
#	 Daphne van Beek, d.vanbeek@vumc.nl  									  #
#	 Roy Straver, r.straver@vumc.nl                                 		  #
#                                                                             #
#    This script is supplementary to WISECONDOR.                       		  #
#                                                                             #
#    WISECONDOR is free software: you can redistribute it and/or 	  		  #
#	 modify it under the terms of the GNU General Public License as 		  #
#	 published by the Free Software Foundation, either version 3 of the 	  #
# 	 License, or (at your option) any later version.                          #
#                                                                             #
#    WISECONDOR is distributed in the hope that it will be useful,     		  #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with WISECONDOR. If not, see <http://www.gnu.org/licenses/>.       #
#                                                                             #
###############################################################################

import glob
import pickle
import matplotlib
matplotlib.use('Agg')
from pylab import *
import argparse
import numpy as np
import os
import time
import re
    
def getCoverage(sample):
	return sum([sum(sample[chrom]) for chrom in sample])


def getYPerc(sample):
	return sum(sample[testChrom])/sum([sum(sample[chrom]) for chrom in sample])


def getYPercMean(sampleList):
	values=[]
	for sample in sampleList:
		values.append(getYPerc(sampleList[sample]))
	#print "Upper limit: " + str(max(values))
	#print "Lower limit: " + str(min(values))
	#return median(values)
	return mean(values)
	
	
def getYPercGrand(sampleList):
	values=[]
	for sample in sampleList:
		values.append(getYPerc(sampleList[sample]))

	return values
	
	
def getGender(prediction):
	#print prediction
	if prediction == [1]:
		return "Male"
	elif prediction == [0]:
		return "Female"
	else:
		return None
	
	
def getColor(gender):
	colors = ["Blue", "HotPink","Cyan", "Violet", "Grey"] 
	pos = ["Male", "Female", "Probably male", "Probably female", "Unknown"].index(gender)
	return colors[pos]
	
	
def solveFetalFraction(percYMales, percYFemales, percYSample):
	#Based on: %chrY sample = meanY% males * FF + meanY% women with female fetusses * (1 - FF)
	#Taken from: Chiu et al, Non-invasive prenatal assessment of trisomy 21 by multiplexed maternal plasma DNA sequencing: large scale validity study, 2011
	return (percYSample - percYFemales) / (percYMales - percYFemales)
	#return (percYSample - percYFemales) / (percYMales)
			
default_fig = "./DEFRAG_out"
print >> sys.stderr, default_fig

parser = argparse.ArgumentParser(description='DEFRAG \
	(DEtection of fetal FRaction And Gender): \
	Determine fetal gender and fraction in a maternal plasma sample. \
	Can be used with or without a pool of male reference samples. \
	It is recommended to use a male reference set that is processed in the normal labflow to get the best results. \
	\nThis tool can be used in addition to WISECONDOR, as it uses two types of WISECONDOR output as input.\n\n \
	Please set up your reference sets as follows:\n \
	\tCreate two/three directories and place the corresponding .gcc and .pickle files in these directories. You should provide:\n \
	\t- Directory with normal pregnancy samples with male fetus\n \
	\t- Directory with normal pregnancy samples with female fetusses\n \
	\t- Optional directory with male reference samples',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('boydir', type=str,
					help='Directory containing fetal boy samples to be used as reference (.gcc and .pickle)')
parser.add_argument('girldir', type=str,
					help='Directory containing fetal girl samples to be used as reference (.gcc and .pickle)')
parser.add_argument('--maledir', type=str,
					help='Directory containing full male samples to be used as reference (.ggc and .pickle)')	
parser.add_argument('--scalingFactor', type=str, help='Factor that is used for correcting the calculated fetal fraction')	
parser.add_argument('--percYonMales', type=str, help='Percentage of reads that is mapped on Y in males')
parser.add_argument('testdir', type=str,
					help='Directory containing test samples (gcc and pickle)')
parser.add_argument('outputfig', type=str, default=default_fig, help='prefix of output figure (extension is added by script)')


args = parser.parse_args()


print '# Script information:'
print '\n# Settings used:'
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
	print '\t'.join([arg,str(argsDict[arg])])


## Load the reference data

print >> sys.stderr, '\n# Processing:'
print >> sys.stderr, 'Loading reference samples'
boySamples	= dict()
boySamplesPickle = dict()
boyFiles = glob.glob(args.boydir + '/*.gcc')
for boyFile in boyFiles:
	print >> sys.stderr, '\tLoading boy gcc:\t' + boyFile
	curFile = pickle.load(open(boyFile,'rb'))
	boySamples[boyFile] = curFile
	pic = os.path.splitext(boyFile)[0] + ".pickle"
	print >> sys.stderr, '\tLoading boy pickle:\t' + pic
	curFile = pickle.load(open(pic,'rb'))
	boySamplesPickle[pic] = curFile

girlSamples	= dict()
girlSamplesPickle = dict()
girlFiles = glob.glob(args.girldir + '/*.gcc')
for girlFile in girlFiles:
	print >> sys.stderr, '\tLoading girl gcc:\t' + girlFile
	curFile = pickle.load(open(girlFile,'rb'))
	girlSamples[girlFile] = curFile
	pic = os.path.splitext(girlFile)[0] + ".pickle"
	print >> sys.stderr, '\tLoading girl pickle:\t' + pic
	curFile = pickle.load(open(pic,'rb'))
	girlSamplesPickle[pic] = curFile

if args.maledir:	
	print >> sys.stderr, 'Found directory with male reference samples.'	
	maleSamples	= dict()
	maleSamplesPickle = dict()
	maleFiles = glob.glob(args.maledir + '/*.gcc')
	for maleFile in maleFiles:
		print >> sys.stderr, '\tLoading man gcc:\t' + maleFile
		curFile = pickle.load(open(maleFile,'rb'))
		maleSamples[maleFile] = curFile
		pic = os.path.splitext(maleFile)[0] + ".pickle"
		print >> sys.stderr, '\tLoading man pickle:\t' + pic
		curFile = pickle.load(open(pic,'rb'))
		maleSamplesPickle[pic] = curFile


## Determine the subset of Y that is used in one of the fetal fraction determinations

testChrom = 'Y'

minLen = min([len(girlSamples[girlSample][testChrom]) for girlSample in girlSamples])
girlData = []
for i in range(minLen):
	girlData.append([girlSamples[girlSample][testChrom][i] for girlSample in girlSamples])

boyData = []
for i in range(minLen):
	boyData.append([boySamples[boySample][testChrom][i] for boySample in boySamples])

removables=[]
keepers=[]
for pos,values in enumerate(girlData):
	if median(values) != 0 or sum(boyData[pos]) == 0: #sum or median
		removables.append(pos)
	else:
		keepers.append(pos)

for i in reversed(removables):
	boyData.pop(i)
	girlData.pop(i)

print >> sys.stderr, 'Bins that are kept for subset Y analysis:'
print >> sys.stderr, keepers


# Load the test data

print >> sys.stderr, 'Loading test samples'
testSamples	= dict()
testSamplesPickle = dict()
testFiles = glob.glob(args.testdir + '/*.gcc')
for testFile in testFiles:
	print >> sys.stderr, '\tLoading test gcc:\t' + testFile
	samplename = os.path.splitext(testFile)[0]
	curFile = pickle.load(open(testFile,'rb'))
	testSamples[samplename] = curFile
	pic = samplename + ".pickle"
	print >> sys.stderr, '\tLoading test pickle:\t' + pic
	curFile = pickle.load(open(pic,'rb'))
	testSamplesPickle[samplename] = curFile

	
## Determine the backgroud values used for correction of the whole chrY fetal fraction determination

percYBoys = getYPercMean(boySamplesPickle)
percYGirls = getYPercMean(girlSamplesPickle)

if args.maledir:
	percYMales = getYPercMean(maleSamplesPickle)
	maleCorMedian = []
	for male in maleSamples:
		corrMales=[maleSamples[male][testChrom][pos] for pos in keepers]
		maleCorMedian.append(median(corrMales))
	corrMalesMedian = mean(maleCorMedian)
else:
	corrMalesMedian = 0.412516803449
	percYMales= 0.00278246251169
	if args.scalingFactor:
		corrMalesMedian = float(args.scalingFactor)
	if args.percYonMales:
		percYMales = float(args.percYonMales)
print >> sys.stderr, 'percYMales:\t' + str(percYMales)
print >> sys.stderr, 'corrMalesMedian:\t' + str(corrMalesMedian)

## Build trainingset for gender determination

training = getYPercGrand(girlSamplesPickle)[:]
training.extend(getYPercGrand(boySamplesPickle)[:])
training = np.array([[x] for x in training])
targets = [0] * len(getYPercGrand(girlSamplesPickle))
targets.extend([1] * len(getYPercGrand(boySamplesPickle)))
targets = np.array(targets)
#from sklearn.naive_bayes import GaussianNB
#gnb = GaussianNB()
#from sklearn.qda import QDA
#gnb = QDA()
from sklearn.neighbors import KNeighborsClassifier
gnb = KNeighborsClassifier(5)
#from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
#gnb= AdaBoostClassifier()
#gnb= RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
#from sklearn.tree import DecisionTreeClassifier
#gnb = DecisionTreeClassifier()
gnb.fit(training, targets)
y_pred = gnb.predict(training)
print >> sys.stderr, "Testing classifier on trainingset.\nNumber of mislabeled points : %d" % (targets != y_pred).sum()

## Start analyzing the test samples and create plots

f, ((ax1,ax2),(ax3,ax4)) = subplots(nrows=2, ncols=2)
axislabels = 8
labelSize = 9
plt.rc('font', **{'size':'8'})

ax1.set_title("Training Set Gender Determination", fontsize=labelSize)
#ax1.set_ylim(0.95 * (min(training))*100, 1.05 * (max(training))*100)
ax1.set_ylim(0.95 * (min(training))*100, 0.05)
ax1.set_ylabel("% of reads on Y chromosome", fontsize=labelSize)
ax1.tick_params(labelsize=axislabels)
for index,val in enumerate(training):
	col = ['HotPink', 'Blue']
	ax1.scatter(0, val*100, c=col[targets[index]])
	#print str(val) + str(targets[index])

sortedList=testSamples.keys()
sortedList.sort()
pdfData = []
headerLine = ["Sample", "DEFRAG subset ChrY", "DEFRAG whole ChrY", "Determined Gender", "Total number of reads", "Cluster", "% on Y"]
print "\t".join(headerLine)
pdfData.append(headerLine)

for index,testSample in enumerate(sortedList):
	result=[testSamples[testSample][testChrom][pos] for pos in keepers]
	votesBoy=len([x for x in result if x != 0])
	votesGirl=len([x for x in result if x == 0])
	
	#Use classifier to predict gender
	prediction = gnb.predict(getYPerc(testSamplesPickle[testSample]))

	#Based on: %chrY sample = meanY% males * FF + meanY% women with female fetusses * (1 - FF)
	daphGender = solveFetalFraction(percYMales, percYGirls, getYPerc(testSamplesPickle[testSample]))

	if median(result)/corrMalesMedian == 0.0 and getGender(prediction) == 'Male':
		cluster = "BAD"
	elif median(result)/corrMalesMedian == 0.0:
		cluster = "Girls"
	else:
		cluster = "Boys"

	lines = [testSample.split('/')[-1], str((median(result)/corrMalesMedian)*100), str(daphGender*100), str(getGender(prediction)), str(getCoverage(testSamplesPickle[testSample])), cluster, str(getYPerc(testSamplesPickle[testSample])*100)]
	print "\t".join(lines)
	pdfData.append(lines)
	
	color = getColor(getGender(prediction))
			
	ax2.set_title("Test Samples Gender Determination", fontsize=labelSize)
	ax2.set_ylabel("% of reads on Y chromosome", fontsize=labelSize)
	ax2.tick_params(labelsize=axislabels)
	ax2.set_ylim(0,0.05)
	ax2.scatter(0, getYPerc(testSamplesPickle[testSample])*100, c=color, marker='o')
	
	ax3.set_title("DEFRAG Script", fontsize=labelSize)
	ax3.set_ylabel("Fetal Fraction (%) on subset of Y", fontsize=labelSize)	#Underestimation
	ax3.set_xlabel("Fetal Fraction (%) on whole Y", fontsize=labelSize)
	ax3.tick_params(labelsize=axislabels)
	ax3.set_xlim(-5,50)
	ax3.set_ylim(-5,50)
	ax3.scatter(daphGender*100, median(result)/corrMalesMedian*100, c=color, marker='o')
	
	ax4.set_title("Fetal Fraction (%) on whole chr Y vs. Coverage", fontsize=labelSize)
	ax4.set_ylabel("Fetal Fraction (%) on whole chr Y", fontsize=labelSize)	#Underestimation
	ax4.set_xlabel("Reads left after filtering", fontsize=labelSize)
	ax4.tick_params(labelsize=axislabels)
	#ax4.scatter(getCoverage(testSamplesPickle[testSample]), median(result)/corrMalesMedian*100, c=color, marker='o')
	ax4.scatter(getCoverage(testSamplesPickle[testSample]), daphGender*100, c=color, marker='o')

ax4.axvline(x=8000000, color='r')
ax4.axvline(x=12000000, color='g')

savefig(args.outputfig + ".png")

exit()

