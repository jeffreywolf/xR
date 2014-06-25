README 
======

General
-------
This program computes empirical spatial correlograms or cross-correlograms and uses bootstrapping to estimate confidence intervals.  

Spatial correlograms can be used to assess spatial autocorrelation of one variable and spatial cross-correlograms can be used to assess spatial codependence of two variables.

A correlation coefficient is calculated at regularly spaced separation distances.  The separation distance grid is by default the mean minimum separation distance between pairs of points, but there is an option to specify a value. Tolerance for inclusion of point pairs at a separation distance is plus or minus half of the first separation distance. Pearson's correlation coefficient is currently supported. The correlogram is the correlation as a function of euclidean separation distance.

*Bootstrapping*
System time and a unique increment are used to seed each bootstrap resample. There is an option (-s) to set the seed. 

Dependencies
------------
SciPy stack: numpy, scipy.spatial, scipy.stats

Data File
---------
This program takes as input a single comma separated value (CSV) format data file. The data file must contain an x coordinate, a y coordinate, and at least one spatially indexed variable.

See *data.csv* for an example of a data file.

Missing data should be specified by NA. All lines that contain NA in the data file are removed.

The spatially indexed variable can be log-transformed using the -l flag as long as it contains no zeroes or negative values.


Options
-----------
For a description of the command line options type:

`./xR.py -h`

However, make sure that the program is executable before doing any calls.

`chmod u+x xR.py`


Output
-------
This program will output three files in CSV format. 

1. An empirical correlogram or cross-correlogram. This file begins with the prefix "xR" then is followed by the first and second variable names.

2. A summary file containing 95% Confidence Intervals at each separation distance. This file begins with "summary-bxR", is followed by the first and second variable names, then ends with the number of bootstrap resamples.

3. A file containing all bootstrap resampled empirical correlograms or cross-correlograms.  The file begins with "bxR", is followed by the first and second variable names, then ends with the number of bootstrap resamples.

The command line flag for output should specify a directory for where to write the output files.


Example call
------------
Calculate a spatial correlogram using default separation distance.

`./xR -d data.csv -x x -y y -a variable -b variable -n 1000 -o xr-output -m 300`


