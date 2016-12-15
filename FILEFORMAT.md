# FILE FORMAT

## Preface

Numpy compressed formats are somewhat unusual for Python users. This document is meant to help you understand what data is saved in the files and how to obtain them should you wish to. Either for custom report creation or for additions to the tool.

Import the required module:

	import numpy

Rid us of some function pointing variable errors:

	toolConvert = '<toolConvert>'
	toolTest    = '<toolTest>'
	toolNewref  = '<toolNewref>'

Load a file:  

	theFile = numpy.load('./afile.npz')

Ask what data is in the file:  

	theFile.keys()

Numpy has a limited specification of what is in the file, a few workarounds are needed.
Sometimes it's fine to just access the data directly, sometimes it's not. This depends on whether or not the data accessed is a default numpy array structure. If it is not, use `.item()`:  

	theFile['theData'].item()

This should return the data in it's original structure.



## Converted file
Data contained in the first npz created from a bam file using `convert`:  
- `runtime`, a `dict` specifying information about the system when the conversion was ran.
- `arguments`, a `dict` with all arguments used to call the script.
- `quality`, a `dict` with integers showing the amount of reads lost per filter.
- `sample`, a `dict` with an `int32` array per chromosome, telling the amount of reads per bin.


	runtime {
		'username': 'rstraver',
		'version':  '5cce5e3\n',
		'hostname': 'marvin',
		'datetime': datetime.datetime(2016, 12, 12, 17, 11, 52, 131386)
	}

	arguments {
		'retthres': 4,
		'retdist':  4,
		'binsize':  50000.0,
		'outfile':  './testSamples/sampleFile.npz',
		'func':     '<toolConvert>',
		'infile':   './testSamples/sampleFile.bam'
	}

	quality {
		'filter_rmdup':  1189893,
		'filter_mapq':   3401000,
		'no_coordinate': 12580916L,
		'mapped':        43366013L,
		'unmapped':      0L,
		'post_retro':    38754546,
		'pair_fail':     0,
		'pre_retro':     43365720
	}

	sample {
		'1':  array([0, ..., 0], dtype=int32),
		'2':  array([0, ..., 0], dtype=int32),
		...,
		'22': array([0, ..., 0], dtype=int32),
		'X':  array([0, ..., 0], dtype=int32),
		'Y':  array([0, ..., 0], dtype=int32)
	}



## Results file
Data contained in the npz generated using `test`:

- `runtime`, a `dict` specifying information about the system when the conversion was ran.
- `arguments`, a `dict` with all arguments used to call the script.
- `binsize`, an `int` that represents the binsize this sample was tested at.
- `asdef`, a `float` that represents the average standard deviation in this sample.
- `aasdef`, a `float`, equals `asdef * threshold_z`, relative effect size required on average to make a call for a single bin using the z-score method.
- `threshold_z`, a `float` representing the treshold used to determine significant variations.
- `results_calls`, an `array` of `arrays` of `float` values, every array defines a call.  
format: `[chromosome, startBin, endBin, zScore, effectSize]`. Cast first 3 values to `int` if desired.
- `results_cwz`, `array` of `float` values representing chromosome wide z-scores, using stouffers z-score to combine all z-scores on a chromosome.
- `results_r`, `array` of `arrays` of `float` values, one array per chromosome, representing relative read depth change per bin.
- `results_z`, `array` of `arrays` of `float` values, one array per chromosome, representing z-score per bin.


	runtime {
		'username': 'rstraver',
		'version':  '5cce5e3\n',
		'hostname': 'marvin',
		'datetime':  datetime.datetime(2016, 12, 12, 17, 20, 21, 492569)
	}

	arguments {
		'minzscore':     None,
		'minrefbins':    25,
		'reference':     './dataFiles/reference.npz',
		'outfile':       './testSamples/sampleFile_out.npz',
		'func':          '<toolTest>',
		'mineffectsize': 0,
		'repeats':       5,
		'infile':        './testSamples/sampleFile.npz',
		'multitest':     1000
	}

	binsize
		250000

	asdef
		0.030506466625100436

	aasdef
		0.15509784520217804

	threshold_z
		5.084097319699588

	results_calls [
		[1.00000000e+00, 5.40000000e+01, 5.40000000e+01,
		 1.06117139e+01,   3.82842947e-01],
		...,
		[2.20000000e+01, 8.60000000e+01, 8.60000000e+01,
		-6.10927573e+00,  -1.46670420e-01]
	]

	results_cwz [
		-2.00089991,
		...,
		1.28552298
	]

	results_r [
		array([0.00000000e+00, ..., 0.00000000e+00]),
		...,
		array([0.00000000e+00, ..., 0.00000000e+00])
	]

	results_z [
		array([0.00000000e+00, ..., 0.00000000e+00]),
		...,
		array([0.00000000e+00, ..., 0.00000000e+00])
	]


## Reference file
- `runtime`, a `dict` specifying information about the system when the conversion was ran.
- `arguments`, a `dict` with all arguments used to call the script.
- `binsize`, an `int` that represents the binsize this reference was created at.
- `chromosome_sizes`, an `array` of `int` values to tell the amount of bins per chromosome.
- `masked_sizes`, an `array` of `int` values to tell the amount of bins per chromosome after masking.
- `mask`, an `array` of `bool` values to tell what bins are masked out for further analyses.
- `pca_components`, an `array` of `arrays` of `float` values describing PCA components.
- `pca_mean`, an `array` of `float` values describing PCA means.
- `indexes`, an `array` of `arrays` of `int` values telling the indexes of bins selected for every bin. Every array contains the refset for one bin, all chromosomes are concatenated into one large array.
- `distances`, an `array` of `arrays` of `float` values telling the distance squared between target and reference bin for every bin selected in `indexes`.


	runtime {
		'username': 'rstraver',
		'version':  'b12fadf\n',
		'hostname': 'marvin',
		'datetime':  datetime.datetime(2016, 12, 8, 15, 40, 19, 693540)
	}

	arguments {
		'infiles':  ['./refSamples/firstRefSample.npz',
		             ...,
		             './refSamples/lastRefSample.npz'],
		'partfile': './dataFiles/reference_part',
		'cpus':     1,
		'binsize':  250000,
		'outfile':  './dataFiles/reference.npz',
		'parts':    1,
		'func':     '<toolNewref>',
		'part':     [1, 1],
		'refsize':  100,
		'prepfile': './dataFiles/reference_prep.npz'
	}

	binsize
		250000

	chromosome_sizes [
		998, ..., 206
	]

	masked_sizes [
		913, ..., 141
	]

	mask [
		True, ..., False
	]

	pca_components [
		[0.00395148, ...,  0.01222231],
		...,
		[0.00389488, ...,  0.01998623]
	]

	pca_mean [  
		1.22552300e-05, ...,   8.10593756e-05
	]

	indexes [
		[395, 446, 426 ..., 466, 467, 468],
		...,
		[315, 377, 422 ..., 331, 338, 340]
	]

	distances [
		[0.00000000e+00, ..., 1.23259516e-32],
		...,
		[5.83621524e-24, ..., 5.86302828e-24]
	]
