# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (github.com/rstraver)
#
# This file is part of WISECONDOR
# WISECONDOR is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import get_cmap
import numpy as np
import time
import bisect
from sklearn.decomposition import PCA
from sklearn.utils.extmath import fast_dot
import pysam
from triarray import *
import subprocess
import getpass
import datetime
import socket

# Get rid of some useless warnings, I know there are emtpy slices
import warnings
warnings.filterwarnings('ignore', 'Mean of empty slice')
warnings.filterwarnings('ignore', 'Degrees of freedom <= 0 for slice')

curTime = datetime.datetime.now()

np.seterr('ignore')
np_sum = np.sum
np_pow = np.power
np_max = np.argmax
np_mean = np.mean
np_median = np.median
np_std = np.std
np_abs = np.abs
np_sqrt = np.sqrt

find_pos = bisect.bisect


def getRuntime():
	runtime = dict()
	runtime['version']=getVersion()
	runtime['datetime']=curTime
	runtime['hostname']=socket.gethostname()
	runtime['username']=getpass.getuser()
	return runtime


def getVersion():
	version = 'unknown'
	try:
		version = subprocess.check_output(["git", "describe", "--always"]).split()[0]
	except:
		pass
	return version


def printArgs(args):
	argdict=vars(args)
	print 'tool =', str(argdict['func']).split()[1][4:]
	for arg in sorted(argdict.keys()):
		if arg != 'func':
			print arg,'=',argdict[arg]


def loadCytoBands(cytoFile):
	cytoDict = dict()
	curChrom = None
	curCyto = None
	with open(cytoFile, 'r') as cytoData:
		for line in cytoData:
			splitLine = line.split()
			if splitLine[0][3:] != curChrom:
				if curChrom != None:
					cytoDict[curChrom] = curCyto
				curCyto = []
				curChrom = splitLine[0][3:]
			curCyto.append(splitLine[1:])
	return cytoDict


def trainPCA(refData,pcacomp=3):
	tData = refData.T
	pca = PCA(n_components=pcacomp)
	pca.fit(tData)
	PCA(copy=True, whiten=False)
	transformed = pca.transform(tData)
	inversed = pca.inverse_transform(transformed)
	corrected = tData / inversed
	#print pca.n_components_
	#print pca.explained_variance_
	#print pca.explained_variance_ratio_

	return corrected.T, pca


def applyPCA(sampleData, mean, components):
	pca = PCA(n_components=components.shape[0])
	pca.components_ = components
	pca.mean_ = mean

	transform = pca.transform(np.array([sampleData]))

	reconstructed = fast_dot(transform, pca.components_) + pca.mean_
	reconstructed = reconstructed[0]
	return sampleData / reconstructed


def convertBam(bamfile, binsize=1000000, minShift=4, threshold=4, mapq=1, demandPair=False):
	# Prepare the list of chromosomes
	chromosomes = dict()
	for chromosome in range(1, 23):
		chromosomes[str(chromosome)] = None
	chromosomes['X'] = None
	chromosomes['Y'] = None

	# Flush the current stack of reads
	def flush(readBuff, counts):
		stairSize = len(readBuff)
		if stairSize <= threshold or threshold < 0:
			for read in readBuff:
				location = read.pos / binsize
				counts[location] += 1
				#if location >= len(counts):
				#	print read

	sam_file = pysam.AlignmentFile(bamfile, "rb")
	reads_seen = 0
	reads_kept = 0
	reads_mapq = 0
	reads_rmdup = 0
	reads_pairf = 0
	larp = -1 # LAst Read Position...
	larp2 = -1

	for index,chrom in enumerate(sam_file.references):

		chromName = chrom
		if chromName[:3].lower() == 'chr':
			chromName = chromName[3:]
		if chromName not in chromosomes:
			continue

		print chrom,'length:', sam_file.lengths[index], 'bins:', int(sam_file.lengths[index] / float(binsize) + 1)
		counts = np.zeros(int(sam_file.lengths[index] / float(binsize) + 1), dtype=np.int32)

		readBuff = []
		sam_iter = sam_file.fetch(chrom)

		prevRead = sam_iter.next()

		# Split paths here, for-loop was heavily slowed down by if-statements otherwise
		if demandPair:
			for read in sam_iter:
				if ((int(read.pos) - int(prevRead.pos)) > minShift):
					flush(readBuff, counts)
					readBuff = []
				# Normal ndups will be appended here

				if not (read.is_proper_pair and read.is_read1):
					reads_pairf += 1
					continue

				if larp == read.pos and larp2 == read.next_reference_start:
						reads_rmdup += 1
				else:
					if read.mapping_quality >= mapq:
						readBuff.append(read)
						prevRead = read
					else:
						reads_mapq += 1

				larp2 = read.next_reference_start

				reads_seen += 1
				larp = read.pos
		else:
			for read in sam_iter:
				if ((int(read.pos) - int(prevRead.pos)) > minShift):
					flush(readBuff, counts)
					readBuff = []
				# Normal ndups will be appended here

				if larp == read.pos:
						reads_rmdup += 1
				else:
					if read.mapping_quality >= mapq:
						readBuff.append(read)
						prevRead = read
					else:
						reads_mapq += 1

				reads_seen += 1
				larp = read.pos

		# Flush after we're done
		flush(readBuff, counts)
		chromosomes[chromName] = counts
		reads_kept += sum(counts)

	#print reads_seen,reads_kept
	qual_info = {'mapped':sam_file.mapped,
				 'unmapped':sam_file.unmapped,
				 'no_coordinate':sam_file.nocoordinate,
				 'filter_rmdup':reads_rmdup,
				 'filter_mapq':reads_mapq,
				 'pre_retro':reads_seen,
				 'post_retro':reads_kept,
				 'pair_fail':reads_pairf}
	return chromosomes, qual_info


def scaleSample(sample, fromSize, toSize):
	if fromSize == toSize or toSize == None:
		return sample

	if toSize == 0 or fromSize == 0 or toSize < fromSize or toSize % fromSize > 0:
		print 'ERROR: Impossible binsize scaling requested:', fromSize, 'to', toSize
		exit(1)

	returnSample = dict()
	scale = toSize/fromSize
	for chrom in sample:
		chromData = sample[chrom]
		newLen = int(np.ceil(len(chromData)/float(scale)))
		scaledChrom = np.zeros(newLen, dtype=np.int32)
		for i in range(newLen):
			scaledChrom[i] = np_sum(chromData[i*scale:i*scale+scale])
			returnSample[chrom] = scaledChrom
	return returnSample


def toNumpyArray(samples):
	byChrom = []
	chromBins = []
	sampleCount = len(samples)
	for chromosome in range(1, 23):
		maxLen = max([sample[str(chromosome)].shape[0] for sample in samples])
		thisChrom = np.zeros((maxLen, sampleCount), dtype=float)
		chromBins.append(maxLen)
		i = 0
		for sample in samples:
			thisChrom[:, i] = sample[str(chromosome)]
			i += 1
		byChrom.append(thisChrom)
	allData = np.concatenate(byChrom, axis=0)

	sumPerSample = np_sum(allData, 0)
	allData = allData / sumPerSample

	print 'Applying nonzero mask on the data:', allData.shape,
	sumPerBin = np_sum(allData, 1)
	mask = sumPerBin > 0
	maskedData = allData[mask, :]
	print 'becomes',maskedData.shape

	return maskedData, chromBins, mask


def toNumpyRefFormat(sample, chromBins, mask):
	byChrom = []
	for chromosome in range(1, 23):
		thisChrom = np.zeros(chromBins[chromosome - 1], dtype=float)
		minLen = min(chromBins[chromosome - 1], len(sample[str(chromosome)]))
		thisChrom[:minLen] = sample[str(chromosome)][:minLen]
		byChrom.append(thisChrom)
	allData = np.concatenate(byChrom, axis=0)
	allData = allData / np_sum(allData)
	maskedData = allData[mask]

	return maskedData


def inflateArray(array, mask):
	temp = np.zeros(mask.shape[0])
	j = 0
	for i, val in enumerate(mask):
		if val:
			temp[i] = array[j]
			j += 1
	return temp


def inflateArrayMulti(array, mask_list):
	temp = array
	for mask in reversed(mask_list):
		temp = inflateArray(temp, mask)
	return temp


def getRefForBins(amount, start, end, sampleData, otherData):
	refIndexes = np.zeros((end - start, amount), dtype=np.int32)
	refDistances = np.ones((end - start, amount))
	for thisBin in xrange(start, end):
		thisMask = np_sum(np_pow(otherData - sampleData[thisBin, :], 2), 1)

		# 209 seconds on 0.5mb size:
		thisIndexes = [-1 for i in xrange(amount)]
		thisDistances = [1e10 for i in xrange(amount)]
		removeIndex = thisIndexes.pop
		removeDist = thisDistances.pop
		insertIndex = thisIndexes.insert
		insertDist = thisDistances.insert
		i = 0
		curMax = 1e10
		for binVal in thisMask:
			if binVal < curMax:
				pos = find_pos(thisDistances, binVal)
				removeIndex(-1)
				removeDist(-1)
				insertIndex(pos, i)
				insertDist(pos, binVal)
				curMax = thisDistances[-1]
				i += 1

		refIndexes[thisBin - start, :] = thisIndexes
		refDistances[thisBin - start, :] = thisDistances
	return refIndexes, refDistances


def getOptimalCutoff(reference, repeats):
	optimalCutoff = float("inf")
	mask = np.zeros(reference.shape)
	for i in range(0, repeats):
		mask = reference < optimalCutoff
		average = np.average(reference[mask])
		stddev = np.std(reference[mask])
		optimalCutoff = average + 3 * stddev
	return optimalCutoff, mask


# Returns: Chromosome index, startBinNumber, endBinNumber
def splitByChrom(start, end, chromosomeBinSums):
	areas = []
	tmp = [0, start, 0]
	for i, val in enumerate(chromosomeBinSums):
		tmp[0] = i
		if val >= end:
			break
		if start < val < end:
			tmp[2] = val
			areas.append(tmp)
			tmp = [i, val, 0]
		tmp[1] = val
	tmp[2] = end
	areas.append(tmp)
	return areas


# Returns: Start and end bin numbers this instance should work on
def getPart(partnum, outof, bincount):
	startBin = int(bincount / float(outof) * partnum)
	endBin = int(bincount / float(outof) * (partnum + 1))
	return startBin, endBin


def getReference(correctedData, chromosomeBins, chromosomeBinSums, selectRefAmount=100, part=1, splitParts=1):
	timeStartSelection = time.time()
	bigIndexes = []
	bigDistances = []

	bincount = chromosomeBinSums[-1]

	startNum, endNum = getPart(part - 1, splitParts, bincount)
	print 'Working on part', part, 'of', splitParts, 'meaning bins', startNum, 'up to', endNum
	regions = splitByChrom(startNum, endNum, chromosomeBinSums)

	for region in regions:
		chrom = region[0]
		start = region[1]
		end = region[2]

		if startNum > start:
			start = startNum
		if endNum < end:
			end = endNum

		print part, 'Actual Chromosome area', chromosomeBinSums[chrom] - chromosomeBins[chrom], chromosomeBinSums[chrom]
		chromData = np.concatenate((correctedData[:chromosomeBinSums[chrom] - chromosomeBins[chrom], :],
									correctedData[chromosomeBinSums[chrom]:, :]))

		partIndexes, partDistances = getRefForBins(selectRefAmount, start, end, correctedData, chromData)
		bigIndexes.extend(partIndexes)
		bigDistances.extend(partDistances)

		print part, 'Time spent:', int(time.time() - timeStartSelection), 'seconds'

	indexArray = np.array(bigIndexes)
	distanceArray = np.array(bigDistances)

	return indexArray, distanceArray


def prepSample(sample, chromosome_sizes, mask, pca_mean, pca_components):
	testData = toNumpyRefFormat(sample, chromosome_sizes, mask)
	testData = applyPCA(testData, pca_mean, pca_components)
	return testData


def trySample(testData, testCopy, indexes, distances, chromosomeBins, chromosomeBinSums, cutoff):
	bincount = chromosomeBinSums[-1]

	resultsZ = np.zeros(bincount)
	resultsR = np.zeros(bincount)
	refSizes = np.zeros(bincount)
	stdDevSum = 0.
	stdDevNum = 0
	i = 0

	for chrom in xrange(len(chromosomeBins)):
		start = chromosomeBinSums[chrom] - chromosomeBins[chrom]
		end = chromosomeBinSums[chrom]
		chromData = np.concatenate(
			(testCopy[:chromosomeBinSums[chrom] - chromosomeBins[chrom]], testCopy[chromosomeBinSums[chrom]:]))

		for index in indexes[start:end]:
			refData = chromData[index[distances[i] < cutoff]]
			refData = refData[refData >= 0]  # Previously found aberrations may be marked by negative values
			refMean = np_mean(refData)
			refStdDev = np_std(refData)
			if not np.isnan(refStdDev):
				stdDevSum += refStdDev
				stdDevNum += 1
			resultsZ[i] = (testData[i] - refMean) / refStdDev
			resultsR[i] = testData[i] / refMean
			refSizes[i] = refData.shape[0]
			i += 1
	return resultsZ, resultsR, refSizes, stdDevSum/stdDevNum


def repeatTest(testData, indexes, distances, chromosomeBins, chromosomeBinSums, cutoff, threshold, repeats):
	timeStartTest = time.time()
	resultsZ = None
	resultsR = None
	testCopy = np.copy(testData)
	for i in xrange(repeats):
		resultsZ, resultsR, refSizes, stdDevAvg = trySample(testData, testCopy, indexes, distances, chromosomeBins,
												  chromosomeBinSums, cutoff)
		testCopy[np_abs(resultsZ) >= threshold] = -1
	print 'Time spent on obtaining z-scores:', int(time.time() - timeStartTest), 'seconds'
	return resultsZ, resultsR, refSizes, stdDevAvg


# TODO: Take care of regions that flip dup/del
def positionsToStretches(positions, maxDist):
	if positions == []:
		return []

	stretches = []
	start = positions[0]
	for i, val in enumerate(positions[:-1]):
		if positions[i + 1] - val > maxDist:
			stretches.append((start, val))
			start = positions[i + 1]
	stretches.append((start, positions[-1]))
	return stretches


def fillTri(region):
	tri_arr = TriArr(region.shape[0])
	addVal = tri_arr.addValue
	for x in xrange(region.shape[0]):
		for y in xrange(x, region.shape[0]):
			addVal(np_sum(region[x:y+1]) / np_sqrt(y - x + 1))
	return tri_arr


def fillTriMin(regionZ,regionR,threshold):
	if threshold == 0:
		return fillTri(regionZ)

	tri_arr = TriArr(regionZ.shape[0])
	addVal = tri_arr.addValue
	for x in xrange(regionZ.shape[0]):
		for y in xrange(x, regionZ.shape[0]):
			if abs(np_median(regionR[x:y+1])-1) >= threshold:
				addVal(np_sum(regionZ[x:y+1]) / np_sqrt(y - x + 1))
			else:
				addVal(0)
	return tri_arr


def segmentThis(region):
	regionLen = region.shape[0]
	champReg = (0, regionLen)
	champVal = abs(np_mean(region)) * np_sqrt(regionLen)
	up = np.median(region) > 1
	for x in xrange(regionLen + 1):
		for y in xrange(x + 1, regionLen + 1):
			thisMean = np_mean(region[x:y])
			optVal = abs(thisMean) * np_sqrt(y - x)
			if optVal > champVal and thisMean > 1 == up:
				champVal = optVal
				champReg = (x, y)
	return champReg


def stouffSeg(region, threshold):
	myResult = []
	champVal = np_sum(region) / np_sqrt(region.shape[0])
	selection = (0, region.shape[0])
	for x in xrange(region.shape[0]):
		for y in xrange(x + 1, region.shape[0]):
			optVal = np_sum(region[x:y]) / np_sqrt(y - x)
			if abs(optVal) > abs(champVal):
				champVal = optVal
				selection = (x, y)
	if abs(champVal) < threshold:
		return myResult
	if selection[0] > 3:
		myResult.extend(stouffSeg(region[:selection[0]], threshold))
	myResult.append((champVal, selection))
	if selection[1] < region.shape[0] - 3:
		rightEnd = stouffSeg(region[selection[1]:], threshold)
		rightEnd = [(x[0], (x[1][0] + selection[1], x[1][1] + selection[1])) for x in rightEnd]
		myResult.extend(rightEnd)
	return myResult


def plotLines(zscores,marks,threshold,sampleName='',binsize=250000,cytoFile=None, chromosomes=range(1,23), columns=2, size=[11.7, 8.3], minEffect=0):

	rows = int(np.ceil(len(chromosomes)/float(columns)))

	colorHorzHelper = (0.7, 0.7, 0.7)
	colorPalette = [
		(0, 0, 0),  # 0 black
		(0.90, 0.60, 0),  # 1 Orange
		(0.35, 0.70, 0.90),  # 2 Sky blue
		(0, 0.60, 0.50),  # 3 Bluish green
		(0.95, 0.90, 0.25),  # 4 Yellow
		(0, 0.45, 0.70),  # 5 Blue
		(0.80, 0.40, 0),  # 6 Vermillion
		(0.80, 0.60, 0.70),  # 7 Reddish purple
	]

	colorMarked = colorPalette[1]
	colorMaternal = colorPalette[3]

	def preparePlot(index, chromnum):
		plt.subplot(rows,columns,index+1)

		frame1 = plt.gca()
		#frame1.axes.get_xaxis().get_major_ticks()
		frame1.axes.xaxis.set_ticklabels([])
		#frame1.axes.get_yaxis().get_major_ticks()
		#frame1.axes.yaxis.set_ticklabels([])
		for tick in frame1.axes.get_yaxis().get_major_ticks():
			tick.label.set_fontsize(4)

		plt.xlim(0,len(zscores[chromnum]))

		plt.axhline(y=0, linewidth=1, color=colorHorzHelper)
		plt.axhline(y=threshold, linewidth=0.75, color=colorHorzHelper)
		plt.axhline(y=-threshold, linewidth=0.75, color=colorHorzHelper)

		move = 0.5
		if cytoFile is not None:
			cytoDict = loadCytoBands(cytoFile)
			for band in cytoDict[str(chromnum+1)]:
				start = float(band[0])/binsize
				end = float(band[1])/binsize
				height=threshold/2
				bottom=-threshold-threshold/2
				plt.axhline(y=-threshold-threshold/2, linewidth=0.75, color=colorHorzHelper)
				alphascale=0.5
				if band[3][1:4] == 'pos':
					alpha = float(band[3][4:])/100
					frame1.add_patch(
						matplotlib.patches.Rectangle((start, bottom),
													 end-start, height,
													 color='black',
													 alpha=alpha* alphascale,
													 linewidth=0.5))
				elif band[3] == 'acen':
					alpha=1
					frame1.add_patch(
						matplotlib.patches.Rectangle((start, bottom),
													 end - start, height,
													 color='black',
													 alpha=alpha * alphascale,
													 hatch='//',
													 linewidth=0.5))
				elif band[3] == 'gvar':
					alpha=0.75
					frame1.add_patch(
						matplotlib.patches.Rectangle((start, bottom),
													 end - start, height,
													 color='black',
													 alpha=alpha * alphascale,
													 hatch='\\',
													 linewidth=0.5))


		zeros = []
		y = 0
		for x in zscores[chromnum]:
			if x == 0:
				zeros.append(y)
			y += 1
		#print zeros
		zeroPatches = positionsToStretches(zeros,1)

		for zeroPatch in zeroPatches:
			frame1.add_patch(
				matplotlib.patches.Rectangle((zeroPatch[0]-0.5, -threshold),
											 zeroPatch[1]-zeroPatch[0]+1, threshold*2,
											 color=colorPalette[7],
											 alpha=0.5,
											 linewidth=0.5))

		for mark in marks:
			if mark[0] == chromnum+1 and abs(mark[4])*100 > minEffect:
				colorTmp = colorMarked
				if abs(mark[4]) >= 0.3:
					colorTmp = colorMaternal
				plt.axvline(x=mark[1]-move, linewidth=0.5, color=colorTmp)
				plt.axvline(x=mark[2]+move, linewidth=0.5, color=colorTmp)
				#plt.axhline(y=mark[3], linewidth=0.5, color=colorTmp)

				vertical_dir = threshold
				if mark[3] < 0:
					vertical_dir = -threshold
				#frame1.add_patch(matplotlib.patches.Rectangle((mark[1]-move,0), mark[2]-mark[1]+2*move, mark[3], facecolor=colorTmp, alpha=min(1,abs(mark[4])*10)))
				frame1.add_patch(
					matplotlib.patches.Rectangle((mark[1] - move, 0), mark[2] - mark[1] + 2 * move, vertical_dir,
												 facecolor=colorTmp, alpha=min(1, abs(mark[4]) * 10)))

				vertical_place = 'bottom'#threshold * 0.75
				if mark[3] > 0:
					vertical_place = 'top' #*= -1
				plt.text(mark[1]+(mark[2]-mark[1])/2, vertical_dir, "{:.1f}".format(mark[3]), fontsize=8,
						 verticalalignment=vertical_place, horizontalalignment='center')

	print 'Plotting Z-Scores'
	#ax = plt.figure(figsize=(11.69, 8.27))
	ax = plt.figure(figsize=(size[0], size[1]))
	#ax.text(0.5, 0.06, 'Chromosomal position in bins', ha='center', va='bottom')
	#ax.text(0.05, 0.5, 'Z-score', ha='left', va='center', rotation='vertical')
	ax.text(0.5, 0.93, 'Z-score versus chromosomal position - Sample ' + sampleName, ha='center', va='bottom')

	for index,chrom in enumerate(chromosomes):
		preparePlot(index,chrom-1)

		plt.plot(zscores[chrom-1],color=colorPalette[5],linewidth=0.5,alpha=1)
		plt.ylabel(chrom)

	rectCall	= plt.Rectangle((0, 0), 1, 1, fc=colorMarked)
	rectNoCall = plt.Rectangle((0, 0), 1, 1, fc=colorPalette[7])
	rectZScores = plt.Rectangle((0, 0), 1, 1, fc=colorPalette[5])
	rectThreshold	= plt.Rectangle((0, 0), 1, 1, fc=colorHorzHelper)

	ax.legend((rectCall,rectNoCall,rectZScores,rectThreshold),
			  ('Called region','Uncallable region','Z-score per bin','Z-score threshold'),
			  'lower center',prop={'size':8}, ncol=2)
	return plt
