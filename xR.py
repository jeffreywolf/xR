#!/usr/bin/env python

"""
xR - a program to measure spatial autocorrelation and spatial 
cross-correlation that uses bootstraping to estimate 95 percent 
confidence intervals.



Author: Jeffrey Wolf
Email: wolfjeff@ucla.edu
Copyright (c) 2014 Jeffrey Wolf, MIT License (MIT)
"""

import argparse, itertools, sys, os, csv, time
import multiprocessing as mp
import numpy as np
from scipy.spatial import distance
from scipy.stats import pearsonr


def getArgs():
	parser = argparse.ArgumentParser(
		description = """Bootstrap correlograms or cross-correlograms"""
	)
	parser.add_argument(
		"-d",
		"--data",
		type = str,
		required = True,
		help = "Data file"
	)
	parser.add_argument(
		"-x",
		type = str,
		required = True,
		help = "The name of the x-coordinate field"
	)
	parser.add_argument(
		"-y",
		type = str,
		required = True,
		help = "The name of the y-coordinate field"
	)
	parser.add_argument(
		"-a",
		"--var1",
		type = str,
		required = True,
		help = "The first variable"
	)
	parser.add_argument(
		"-b",
		"--var2",
		type = str,
		required = True,
		help = "The second variable"
	)
	parser.add_argument(
		"-z",
		"--zsep",
		type = float,
		help = """Separation distance for grid"""
	)
	parser.add_argument(
		"-l",
		"--useLogs",
		action = "store_true",
		default = False,
		help = "Including -l will log-transform non-spatial fields"
	)
	parser.add_argument(
		"-n",
		"-samples",
		type = int,
		help = """Number of bootstrap replicates"""
	)
	parser.add_argument(
		"-o",
		"--output",
		type = str,
		required = True,
		help = """Directory to hold output files."""
	)

	parser.add_argument(
		'-m',
		'--maxsep',
		type = float,
		required = False,
		help = """Used to specify the maximum lag distance analyzed"""
	)

	parser.add_argument(
		'-s',
		'--seed',
		type = int,
		required = False,
		help = """Set seed for the pseudo-random number generator. \
By default system time is used."""
	)

	return parser.parse_args()

def getIndex(header, item):
	for i, elem in enumerate(header):
		if elem.lower() == item.lower():
			return i
	return None


def getData(path, fields):
	""" reads in a CSV file and returns data table for analysis """
	data = []
	with open(path, 'rUb') as f:
		indata = csv.reader(f)
		var_indices = []
		for i, line in enumerate(indata):
			if i == 0:
				header = line
				continue
			if 'NA' in line:
				continue # remove lines with NA's
			data.append(line)
	npdata = np.array(data, dtype=np.float64)
	#print header
	indices = np.array([getIndex(header, item) for item in fields])
	header = [header[i] for i in indices]
	npdata = npdata[:, indices]
	return header, npdata


def set_grid(data, maxsep, zsep):
	"""
	Input: shape (N, 2) numpy array with rows observations
	and columns gx and gy
	
	Return: grid for lag distances of correlogram
	"""
	D = distance.cdist(data, data, 'euclidean')
	#print D
	D_masked = np.ma.masked_equal(D, 0.0, copy=False)
	minimums = np.amin(D_masked, axis=0)
	#maximums = np.amax(D_masked, axis=0)
	#mean_max = np.mean(maximums)
	maximum = np.amax(D_masked)
	if zsep is not None:
		sep_dist = zsep
	else:
		sep_dist = np.mean(minimums)
	if maxsep is not None:
		grid = np.arange(0, maxsep, sep_dist)
	else:
		grid = np.arange(0, maximum, sep_dist)
	return grid


def correlogram(data, grid):
	"""
	Input: shape (N,4) numpy array with rows observations
	and columns gx, gy, feature1, and feature2
	
	Return: Correlogram with separation distance, pairs
	of points, and pearson correlation coefficient (r)
	"""	
	assert data.shape[1] == 4, 'cross correlogram requires N x 4 np.array.'
	# Test for whether to calculate correlogram or cross-correlogram
	if np.array_equal(data[:,2], data[:,3]):
		usage = 'r'
	else:
		usage = 'xr'
	x_coord = 0
	y_coord = 1
	coords = np.array((x_coord,y_coord))
	feature_1_index = 2 # used below to select feature 1 from np.array
	feature_2_index = 3 # used below to select feature 2 from np.array
	result = np.empty(np.array([0,0])) # python list used to return [cross-]correlogram
	# Euclidean Distance Matrix: N x N Matrix
	D = distance.cdist(data[:,coords], data[:,coords], 'euclidean')
	# Set tolerance to half times first separation distance
	# tol is tolerance inclusion when +/- from sep distance h
	tol = 0.5 * grid[1] 
	# traverse the grid and calculate r at each separation interval
	for i, h in enumerate(grid):
		if i == 0 and usage == 'r':
			continue
		# Symmetric boolean mask for distance matrix at specific grid separation distance
		h_mask = ((h-tol) <= D) & (D < (h+tol))
		# Empty feature vectors 
		feature_1 = np.empty([0], dtype = np.float64)
		feature_2 = np.empty([0], dtype = np.float64)
		# j will be the row number of h_mask. Iteration is over rows.
		for j, row in enumerate(h_mask):
			reps = j * row[row == 1]
			feature_1 = np.append(feature_1, data[reps, feature_1_index])
			feature_2 = np.append(feature_2, data[row, feature_2_index])
		pairs = feature_1.shape[0] # could have used feature_2 instead.
		# There is double counting of points
		r = pearsonr(feature_1, feature_2)[0] #0th item is r, 1st item is p-value
		vec = np.array([h, pairs, r], dtype = np.float64)
		if i == 0 and usage =='xr':
			result = vec
			continue
		if i == 1 and usage == 'r':
			result = np.array(vec, dtype =np.float64)
			continue
		result = np.vstack((result, vec))
	return result


def logTransform(data):
	# Implement zero values check here! Cannot log-transform zeroes
	# may be better way to handle zeroes instead of aborting program!
	var1 = 2 # index
	var2 = 3 # index
	indices = np.array([var1,var2])
	try:
		data[:,indices] = np.log(data[:,indices])
	except Exception as e:
		print "Cannot log transform both variables because of zeroes in these data."
		raise Exception
	return data


def getFilenames(header, useLogs, n):
	var1 = 2 # index
	var2 = 3 # index
	print header[var1], header[var2]
	if useLogs == True:
		rf = 'xR-log{0}-log{1}.csv'.format(header[var1],header[var2])
		brf  = 'bxR-log{0}-log{1}-{2}-iters.csv'.format(header[var1],header[var2], n)
		sbrf = 'summary-bxR-log{0}-log{1}-{2}-iters.csv'.format(header[var1],header[var2], n)
	else:
		rf = 'xR-{0}-{1}.csv'.format(header[var1], header[var2])
		brf = 'bxR-{0}-{1}-{2}-iters.csv'.format(header[var1], header[var2], n)
		sbrf = 'summary-bxR-{0}-{1}-{2}-iters.csv'.format(header[var1], header[var2], n)
	return rf, brf, sbrf


def map_func(bootstrap_iter, seed, data, grid):
	print "Bootstrap iteration: {}".format(bootstrap_iter)
	newseed = seed+bootstrap_iter 
	np.random.seed(seed+bootstrap_iter) 
	n = data.shape[0] 
	indices = np.array(np.arange(n))
	resampled_indices = np.random.choice(indices, n, replace = True)
	resampled_data = data[resampled_indices,:]
	result = correlogram(resampled_data, grid)
	nrows = result.shape[0]
	index_vec = np.repeat(bootstrap_iter, nrows).reshape((nrows,1))
	result = np.append(index_vec, result, axis=1)
	return result


def map_func_star(a_b):
	"""mapstar"""
	return map_func(*a_b)


def reduce_func(a, b):
	a = np.vstack((a,b))
	return a


def summarize(data, grid, n, usage):
	data = np.dstack(data)
	bootr= data[:,3,:]
	lower_CI = np.percentile(bootr, 2.5, axis=1)
	upper_CI = np.percentile(bootr, 97.5, axis=1)
	mean = np.mean(bootr, axis=1)
	print data, bootr
	print lower_CI
	print upper_CI
	print mean
	# Correlograms do not have a zero separation distance but cross-correlograms do.
	if usage == 'xR':
		iters = np.repeat(n, grid.shape[0]) # grid.shape[0] give length
		results = np.column_stack((iters, grid, mean, lower_CI, upper_CI))
	else:
		iters = np.repeat(n, grid.shape[0]-1) # grid.shape[0] give length
		results = np.column_stack((iters, grid[1:], mean, lower_CI, upper_CI))
	return results


def writeOut(data, header, filename):
	print "\nInitializing {0} file output".format(filename)
	with open(filename, "w") as f:
		header_str = ",".join(header)+"\n"
		f.write(header_str)
		for line in data: 
			row = ",".join([str(elem) for elem in line])+"\n"
			f.write(row)
	print "Wrote {0} to disk\n".format(filename)


def main():
	#np.random.seed(int(time.time())) # uncomment after testing
	#np.random.seed(0) # comment out after testing

	args = getArgs()

	fields = [
		args.x, 
		args.y, 
		args.var1, 
		args.var2
	]

	header, data = getData(
		args.data, 
		fields
	)

	#print header, data
	initDir = os.getcwd()
	# Change to output directory
	os.chdir(args.output)
	
	# checks to see whether calculating correlogram or cross-correlogram
	if args.var1 == args.var2:
		usage = "R"
	else:
		usage = "xR"

	# Log-transform features
	if args.useLogs == True:
		print "Using logs"
		data = logTransform(data)

	# * Log data seem fine
	# Set names of output files
	if args.n <=100:
		print "Setting number of resamples to 100"
		n = 100
	else:
		n = args.n
	rf, brf, sbrf = getFilenames(header, args.useLogs, n)

	# Set lag distance grid
	# update this function for specifying hmax
	coords = data[:, np.array([0,1])]
	grid = set_grid(coords, args.maxsep, args.zsep)

	# Calculate correlogram
	result = correlogram(data, grid)
	correlogram_header = ["h","pairs","r"]
	writeOut(result, correlogram_header, rf)

	bootstrap_iter = range(n)
	if args.seed is not None:
		seed = [args.seed]*n
	else:
		seed = [int(time.time())]*n
	cores = mp.cpu_count()
	pool = mp.Pool(processes=cores)
	outdata = pool.map(
		map_func_star,
		itertools.izip(
			bootstrap_iter,
			seed,
			itertools.repeat(data),
			itertools.repeat(grid)
		)
	)

	results = reduce(
		reduce_func,
		outdata
	)

	#print results
	bootstraping_header = ["iter","h","pair","r"]
	writeOut(results, bootstraping_header, brf)

	summary_header = ["iters","h","mean","lower_95CI","upper_95CI"]
	summary = summarize(outdata, grid, args.n, usage)
	writeOut(summary, summary_header, sbrf)

	# Change to directory the executable
	os.chdir(initDir)	


if __name__ == '__main__':
	ti = time.time()
	main()
	tf = time.time()
	dt = tf-ti
	print "bootstrap-correlogram.py program is finished."
	print "Total time elapsed is {0} seconds".format(dt)
	print "Cheers!\n"
