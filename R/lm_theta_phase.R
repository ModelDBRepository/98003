# Plot theta phase as a function of location using
# a circular-linear regression fit. 

# Copyright 2007 John L Baker. All rights reserved.
# This software is provided AS IS under the terms of the Open Source
# MIT License. See http://www.opensource.org/licenses/mit-license.php.

# See testcase test_baker_670.cpp for details of how the input 
# data files are created. This form of regression is implemented
# in the package "circular", which is available from a CRAN mirror site.

# The dirPath value used below is only an example. 
# The name encodes key parameters used in the simulation 
# so that different runs results can be kept straight 
# (i.e. quick and dirty metadata).

# Example invocations of the function:
# lm_theta_phase(cellId=32546,spikesFile="test-baker-ca3-spikes.txt")
# lm_theta_phase(tlim=c(540,720),dirPath="J:\\GMU\\Data\\Data0306-670T540-00000-50-20")

# Ensure required packages are present
if (!require(circular)) {
	stop("Package circular was not loaded")
}

# Define a function to do the circular regression of theta phase with x coordinate
lm_theta_phase <- function (
	cellId=0,						# Id of cell to plot
	tlim=c(0,Inf),					# Time range to process (in sec)
	xlim=c(-25,25),					# Range of X-coord processed
	ylim=c(-20,-15),					# Range of Y-coord processed
	spikesFile="test-baker-spikes.txt",		# File containing recorded spikes
	mouseFile="test-baker-mouse-states.txt",	# File containing mouse states
	drive="C:",						# Drive of file system
	dirPath=file.path(drive,"BNSF_1_0","Data","Default")) # High-level path to data dir	
{
	# X,Y coordinate values to be processed. The default
	# is a rectangular region along a linear path.
	xlim1 <- xlim[1];
	xlim2 <- xlim[2];
	ylim1 <- ylim[1];
	ylim2 <- ylim[2];

	# Convert time ranges to msec to go with data read
	tlim1 <- tlim[1]*1000;
	tlim2 <- tlim[2]*1000;
	
	# Read target cell spikes and mouse states
	spikes<-read.csv(file.path(dirPath,spikesFile),
		col.names=c("id","t"));

	mouseStates<-read.csv(file.path(dirPath,mouseFile),
		col.names=c("id","t","x","y","hx","hy"));

	# Select just the cell id specified
	spikes<-spikes[spikes$id==cellId,,drop=FALSE];
	if (nrow(spikes)==0) {
		stop("No spikes found for the cell indicated");
	}

	# Convert spikes to theta phase assuming 8 Hz theta
	spikes$theta<-360*((spikes$t+62.5)%%125-62.5)/125;

	# Merge mouse states and spikes to get locations and heading
	ks<-1;
	spikes$x<-0;
	spikes$y<-0;
	spikes$hx<-0;
	spikes$hy<-0;
	for (km in 1:(nrow(mouseStates))) {
		while (ks<=nrow(spikes) && spikes$t[ks]<=mouseStates$t[km]) {
			spikes$x[ks]<-mouseStates$x[km];
			spikes$y[ks]<-mouseStates$y[km];
			spikes$hx[ks]<-mouseStates$hx[km];
			spikes$hy[ks]<-mouseStates$hy[km];
			ks<-ks+1;
		}
	}

	# Take only spikes with the right spatial locations and times
	sel<-spikes$x>=xlim1 & spikes$x<=xlim2 & 
		spikes$y>=ylim1 & spikes$y<=ylim2 &
		spikes$t>=tlim1 & spikes$t<=tlim2;
	tp<-spikes$theta[sel];
	x<-spikes$x[sel];
	
	# Do a circular regression fit
	# Note that x is adjusted to be zero mean
	xmean<-mean(x);
	tpfit<-lm.circular(as.circular(tp,units="degrees"),
		x-xmean,0,type="c-l",verbose=T);
	print(tpfit);

	beta<-tpfit$coefficients;
	mu<-tpfit$mu;

	cat("mu=",mu,"\nxmean=",xmean,"\nbeta=",beta,"\n");

	# Unwrap the points for plotting.
	tp.fitted = mu+atan(beta*(x-xmean))*360/pi;
	tp.unwrapped = tp-360*round((tp-tp.fitted)/360);

	plot(x,tp.unwrapped,type="p",
		ylab="Theta phase (deg)",xlab="X (cm)");
	points(x[order(x)], tp.fitted[order(x)], type='l')
	
}
