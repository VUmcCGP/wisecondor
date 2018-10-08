# WisecondorX

import os
import sys
import json
import subprocess
import logging
import numpy as np


'''
Scales the bin size of a sample.npz to the one  
requested for the reference
'''

def scale_sample(sample, from_size, to_size):
	if not to_size or from_size == to_size:
		return sample

	if to_size == 0 or from_size == 0 or to_size < from_size or to_size % from_size > 0:
		logging.critical('Impossible binsize scaling requested: {} to {}'.format(int(from_size),
																				 int(to_size)))
		sys.exit()

	return_sample = dict()
	scale = to_size / from_size
	for chr_name in sample:
		chr_data = sample[chr_name]
		new_len = int(np.ceil(len(chr_data) / float(scale)))
		scaled_chr = np.zeros(new_len, dtype=np.int32)
		for i in range(new_len):
			scaled_chr[i] = np.sum(chr_data[int(i * scale):int(i * scale + scale)])
			return_sample[chr_name] = scaled_chr
	return return_sample


'''
Levels gonosomal reads with the one at the autosomes.
'''

def gender_correct(sample, gender):
	if gender == 'M':
		sample['23'] = sample['23'] * 2
		sample['24'] = sample['24'] * 2

	return sample


'''
Communicates with R. Outputs new json dictionary,
resulting from R, if 'outfile' is a key in the
input json. 'infile' and 'R_script' are mandatory keys
and correspond to the input file required to execute the
R_script, respectively.
'''

def exec_R(json_dict):
	json.dump(json_dict, open(json_dict['infile'], 'w'))

	r_cmd = ['Rscript', json_dict['R_script'],
			 '--infile', json_dict['infile']]
	logging.debug('CBS cmd: {}'.format(r_cmd))

	try:
		subprocess.check_call(r_cmd)
	except subprocess.CalledProcessError as e:
		logging.critical('CBS failed with error {}'.format(e))
		sys.exit()

	os.remove(json_dict['infile'])
	if 'outfile' in json_dict.keys():
		json_out = json.load(open(json_dict['outfile']))
		os.remove(json_dict['outfile'])
		return json_out


'''
Calculates z-score.
'''

def get_z_score(results_c, results_rlm):
	zs = []
	for i in range(len(results_c)):
		segment_rlm = results_rlm[results_c[i][0]][results_c[i][1]:results_c[i][2]]
		null_sd = np.std([x for x in segment_rlm if x != 0])
		zs.append(results_c[i][3] / null_sd)
	return zs


'''
Returns MSV, measure for sample-wise noise.
'''

def get_median_segment_variance(results_c, results_r):
	vars = []
	for segment in results_c:
		segment_r = results_r[segment[0]][int(segment[1]):int(segment[2])]
		segment_r = [x for x in segment_r if x != 0]
		if segment_r:
			var = np.var(segment_r)
			vars.append(var)

	return np.nanmedian(vars)


'''
Converts to string, across different python versions
-- tested on py2.7 vs py3.6
'''

def my_str(chars):
	try:
		chars = str(chars, 'utf-8').rstrip('\x00')
	except:
		chars = str(chars).format('utf-8')
	return chars