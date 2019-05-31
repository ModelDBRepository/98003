# Plot some results related to theta phase precession.

# Copyright 2007 John L Baker. All rights reserved.
# This software is provided AS IS under the terms of the Open Source
# MIT License. See http://www.opensource.org/licenses/mit-license.php.

# This is a small subset of the functions available in the MATLAB
# plot_mouse.m file, but gives the general idea of how something similar
# can be done in R. See also lm_theta-phase for circular regression
# of theta phase along a track.

# Functions implemented here are:
# plot_mouse_path, plot_firing_map, plot_linear_rates, 
# plot_theta_phase, read_mouse_spikes

# Here are some example invocations (customize as needed):
# plot_mouse_path(dirPath="J:\\GMU\\Data\\Data0306-670T540-00000-50-20");
# plot_firing_map(dirPath="J:\\GMU\\Data\\Data0306-670T540-00000-50-20");
# plot_linear_rates(tlim=c(540,720),dirPath="J:\\GMU\\Data\\Data0306-670T540-00000-50-20");
# plot_theta_phase(tlim=c(540,720),dirPath="J:\\GMU\\Data\\Data0306-670T540-00000-50-20");
# plot_theta_phase(cellId=32546,spikesFile="test-baker-ca3-spikes.txt")

# Plot the path taken by the mouse over time. For the theta phase precession
# case, the path used is pretty boring, but here it is.
plot_mouse_path<-function(
	tlim=c(0,Inf),					# Time range to process (in sec)
	mouseFile="test-baker-mouse-states.txt",	# File containing mouse states
	dirPath=file.path("C:","BNSF_1_0","Data","Default")) # Data directory	
{
	# Read mouse location and plot path.
	# Here we assume a rectangular maze.
	mouseStates<-read.csv(file.path(dirPath,mouseFile),
		col.names=c("id","t","x","y","hx","hy"));

	mouseStates<-mouseStates[
		mouseStates$t>=tlim[1]*1000 &
		mouseStates$t<=tlim[2]*1000, ];

	x11();
	plot(mouseStates$x,mouseStates$y,type='l',
		xlab="X (cm)",ylab="Y (cm)",
		main="Mouse path");
}

# Plot a firing map. This is also not so interesting for the
# rectangular path used in theta phase precession.
plot_firing_map<-function(
	cellId=0,						# Id of cell to plot
	tlim=c(0,Inf),					# Time range to process (in sec)
	spikesFile="test-baker-spikes.txt",		# File containing recorded spikes
	mouseFile="test-baker-mouse-states.txt",	# File containing mouse states
	dirPath=file.path("C:","BNSF_1_0","Data","Default")) # Data directory	
{
	# Read the firing data
	rms<-read_mouse_spikes(dirPath=dirPath,
		spikesFile=spikesFile,
		mouseFile=mouseFile,
		cellId=cellId,
		tlim=tlim);
	mouseStates<-rms$mouse;
	spikes<-rms$spikes;

	# Accumulate counts by spatial cell. The size of the
	# maze is assumed to be over -20 to 25 cm in x and y
	# and the cell size is assumed to be 5. All spikes
	# read must fall within this size in x and y values.

	xctr<-seq(-22.5,22.501,5);
	yctr<-seq(-22.5,22.501,5);
	spkCnt<-matrix(0,nrow=10,ncol=10);
	occCnt<-matrix(0,nrow=10,ncol=10)

	# Count spikes within a spatial bin
	nspike<-nrow(spikes);
	for (k in 1:nspike) {
		ix<-1+floor((spikes$x[k]+25)/5);
		iy<-1+floor((spikes$y[k]+25)/5);
		spkCnt[ix,iy]<-spkCnt[ix,iy]+1;
	}

	# Count mouse occupancy in a cell
	for (k in 1:nrow(mouseStates)) {
		ix<-1+floor((mouseStates$x[k]+25)/5);
		iy<-1+floor((mouseStates$y[k]+25)/5);
		occCnt[ix,iy]<-occCnt[ix,iy]+1;
	}

	# Get rates adjusted for the mouse sample rate
	mouseSampleRate<-(nrow(mouseStates)-1)/diff(range(mouseStates$t))*1000;
	rates<-spkCnt/(occCnt/mouseSampleRate);
	maxRate<-max(as.vector(rates),na.rm=TRUE);

	# Plot a firing rate map with actual points of firing added.
	# Obviously this can be refined based on the need.	
	x11()
	image(xctr,yctr,rates[nrow(rates):1,],col=topo.colors(32),
		xlab="X (cm)",ylab="Y (cm)",
		main=paste("Maximum rate =",format(maxRate,digits=3),"Hz"))
	points(spikes$x,spikes$y,pch=19)
		
}

# Plot firing rates along a linear path
plot_linear_rates<-function(
	cellId=0,						# Id of cell to plot
	tlim=c(0,Inf),					# Time range to process (in sec)
	xlim=c(-25,25),					# X-coord range to plot
	ylim=c(-20,-15),					# Y-coord range to plot
	spikesFile="test-baker-spikes.txt",		# File containing recorded spikes
	mouseFile="test-baker-mouse-states.txt",	# File containing mouse states
	dirPath=file.path("C:","BNSF_1_0","Data","Default")) # Data directory
{
	# Read the firing data restricting the mouse location and times
	rms<-read_mouse_spikes(dirPath=dirPath,
		spikesFile=spikesFile,
		mouseFile=mouseFile,
		cellId=cellId,tlim=tlim,
		xlim=xlim,ylim=ylim);
	mouseStates<-rms$mouse;
	spikes<-rms$spikes;

	# Accumulate counts by spatial cell. Only the x coord
	# is used in the spatial location. The maze is assumed
	# to extend from x=-25 to x=25. The bin size is 2.5 cm.

	xctr<-seq(-23.75,23.7501,2.5);
	spkCnt<-numeric(20);
	occCnt<-numeric(20)

	# Count spikes within a spatial bin
	nspike<-nrow(spikes);
	for (k in 1:nspike) {
		ix<-1+floor((spikes$x[k]+25)/2.5);
		spkCnt[ix]<-spkCnt[ix]+1;
	}

	# Count mouse occupancy in a cell
	for (k in 1:nrow(mouseStates)) {
		ix<-1+floor((mouseStates$x[k]+25)/2.5);
		occCnt[ix]<-occCnt[ix]+1;
	}

	# Get rates adjusted for the mouse sample rate
	mouseSampleRate<-(nrow(mouseStates)-1)/diff(range(mouseStates$t))*1000;
	rates<-spkCnt/(occCnt/mouseSampleRate);
	maxRate<-max(rates,na.rm=TRUE);

	# Plot the results
	x11()
	plot(xctr,rates,type='b',xlab="X (cm)",ylab="Firing rate (Hz)",
		main=paste("Firing rates for",tlim[1],"to",tlim[2],"sec"));
}

# Plot the theta phase of target place cell firing as a function of
# location (x coordinate). This is a simple version in which there 
# is no attempt to fit a phase precession rate or unwrap theta phases.
# See lm_theta_phase for a more sophisticate treatment.
plot_theta_phase<- function(
	cellId=0,						# Id of cell to plot
	tlim=c(0,Inf),					# Time range to process (in sec)
	xlim=c(-25,25),					# X-coord range to plot
	ylim=c(-20,-15),					# Y-coord range to plot
	thetaFreq=8,					# assumed theta frequency
	spikesFile="test-baker-spikes.txt",		# File containing recorded spikes
	mouseFile="test-baker-mouse-states.txt",	# File containing mouse states
	dirPath=file.path("C:","BNSF_1_0","Data","Default")) # Data directory
{
	# Read the firing data restricting the mouse location and times
	rms<-read_mouse_spikes(
		dirPath=dirPath,
		spikesFile=spikesFile,
		mouseFile=mouseFile,
		cellId=cellId,tlim=tlim,
		xlim=xlim,ylim=ylim);
	mouseStates<-rms$mouse;
	spikes<-rms$spikes;

	# Plot the theta phase (based on 8 Hz theta). Spikes are shown
	# twice, once with a theta phase in the range -180 to 180 and the
	# other in the range 180 to 540 degrees.
	plot(spikes$x,spikes$theta,col='black',ylim=c(-180,540),
		xlab="X (cm)",ylab="Theta phase (deg)",
		main=paste("Place cell theta phase for",tlim[1],"to",tlim[2],"sec"));
	points(spikes$x,spikes$theta+360,col='blue');
}

# Define a helper function to merge mouse location states 
# and target place cell firing information by time. Value 
# returned is a list containing the mouse states and place
# cell spike times augmented with location, heading, and phase.
read_mouse_spikes<-function(dirPath,spikesFile,mouseFile,
	cellId=0,tlim=c(0,Inf),xlim=c(-Inf,Inf),ylim=c(-Inf,Inf),thetaFreq=8)
{
	# Read target cell spikes and mouse states subsetting
	# spikes to the cell specified.
	spikes<-read.csv(file.path(dirPath,spikesFile),
		col.names=c("id","t"));
	spikes<-spikes[spikes$id==cellId,,drop=FALSE];
	if (nrow(spikes)==0) {
		warning("No spikes read");
	}

	mouseStates<-read.csv(file.path(dirPath,mouseFile),
		col.names=c("id","t","x","y","hx","hy"));

	# Convert spikes to theta phase assuming 8 Hz theta (125 ms/cycle)
	thetaCycle<-1000/thetaFreq;
	spikes$theta<-360/thetaCycle*((spikes$t+thetaCycle/2)%%thetaCycle-thetaCycle/2);

	# Merge mouse states and spikes to get locations and heading
	ks<-1;
	spikes$x<-0;
	spikes$y<-0;
	spikes$hx<-0;
	spikes$hy<-0;
	for (km in 1:nrow(mouseStates)) {
		while (ks<=nrow(spikes) && spikes$t[ks]<=mouseStates$t[km]) {
			spikes$x[ks]<-mouseStates$x[km];
			spikes$y[ks]<-mouseStates$y[km];
			spikes$hx[ks]<-mouseStates$hx[km];
			spikes$hy[ks]<-mouseStates$hy[km];
			ks<-ks+1;
		}
	}

	# Return spikes within the limits specified. Note that time limits
	# are specified in seconds and must be converted to msec.
	sel1<-mouseStates$x>=xlim[1] & mouseStates$x<=xlim[2] & 
		mouseStates$y>=ylim[1] & mouseStates$y<=ylim[2] &
		mouseStates$t>=tlim[1]*1000 & mouseStates$t<=tlim[2]*1000;
	sel2<-spikes$x>=xlim[1] & spikes$x<=xlim[2] & 
		spikes$y>=ylim[1] & spikes$y<=ylim[2] &
		spikes$t>=tlim[1]*1000 & spikes$t<=tlim[2]*1000;

	return(list(mouse=mouseStates[sel1,],spikes=spikes[sel2,]));
}
	