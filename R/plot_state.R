# Example function for plotting output from a state recorder.

# Copyright 2007 John L Baker. All rights reserved.
# This software is provided AS IS under the terms of the Open Source
# MIT License. See http://www.opensource.org/licenses/mit-license.php.

# Example invocations (see test cases test_baker_100 and test_baker_110)
# plot_state("test-baker-voltages.txt",varToPlot="Soma",tlim=c(1000,1200))
# plot_state("test-baker-calcium.txt",varToPlot=c("Soma_CaPool","Soma_CaMDP"),tlim=c(1000,1200))
# plot_state("test-baker-currents.txt",varToPlot="D0060",tlim=c(1000,1040))

plot_state<-function(
	file,				# File containing data to plot
	idToPlot=NULL,		# id(s) to plot (NULL matches any id)
	varToPlot=NULL,		# Column name(s) to plot (NULL for all)
	tlim=c(0,Inf),		# time limits to plot (in ms)
	drive="C:",			# Drive or file system
	dirPath=file.path(drive,"BNSF_1_0","Data","Default")) # Data directory
{
	# Read the data to be plotted. Note that the first
	# record is a header that provides column names.
	# The first two columns are id (cell of component id) and t (time in ms).
	S<-read.csv(file.path(dirPath,file));

	# If a time or id range specified, make the selection
	if (!is.null(varToPlot)) {
		colsToPlot<-names(S) %in% varToPlot;
		colsToPlot[1:2]<-TRUE;
		if (sum(colsToPlot)<=2) {
			stop("No matching variables selected for plot");
		}
		S<-S[,colsToPlot,drop=FALSE];
	}
	if (tlim[1]!=0 || is.finite(tlim[2])) {
		S<-S[S$t>=tlim[1] & S$t<=tlim[2],,drop=FALSE]
	}
	if (!is.null(idToPlot)) {
		S<-S[S$id %in% idToPlot,,drop=FALSE];
	}

	if (nrow(S)==0) {
		stop("No matching data to plot");
	}

	# Do the plot
	x11()
	matplot(S$t,S[,3:ncol(S)],type='l',xlab="Time (ms)",
		ylab=ifelse (length(varToPlot)==1,varToPlot,"Value") );
}
