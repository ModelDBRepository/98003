# Example function for plotting output from a spike recorder.

# Copyright 2007 John L Baker. All rights reserved.
# This software is provided AS IS under the terms of the Open Source
# MIT License. See http://www.opensource.org/licenses/mit-license.php.

# Example invocations (see test cases test_baker_100 and test_baker_110)
# plot_spikes("test-baker-spikes.txt")

plot_spikes<-function(
	file,				# File containing data to plot
	idToPlot=NULL,		# id(s) to plot (NULL matches any id)
	tlim=c(0,Inf),		# time limits to plot (in ms)
	drive="C:",			# Drive or file system
	dirPath=file.path(drive,"BNSF_1_0","Data","Default")) # Data directory
{
	# Read the data to be plotted. Note that the first
	# record is a header that provides column names.
	# The first two columns are id (cell of component id) and t (time in ms).
	S<-read.csv(file.path(dirPath,file));

	# If a time or id range specified, make the selection
	if (tlim[1]!=0 || is.finite(tlim[2])) {
		S<-S[S$spike_time>=tlim[1] & S$spike_time<=tlim[2],,drop=FALSE]
	}
	if (!is.null(idToPlot)) {
		S<-S[S$id %in% idToPlot,,drop=FALSE];
	}

	if (nrow(S)==0) {
		stop("No matching data to plot");
	}

	# Do the plot
	plot(S$spike_time,S$id,xlab="Time (ms)",ylab="Spiking cell identifier");
}
