# Example R function for reading and plotting rallpack results.

# Copyright 2007 John L Baker. All rights reserved.
# This software is provided AS IS under the terms of the Open Source
# MIT License. See http://www.opensource.org/licenses/mit-license.php.

# Example invocations
# plot_rallpack_voltages()
# plot_rallpack_states()

plot_rallpack_voltages<-function(
	tlim=c(0,Inf),	# time limits to plot (in ms)
	drive="C:",		# drive (or file system) containing data	
	dirPath=file.path(drive,"BNSF_1_0","Data","Default")) # data directory
{
	# Plot voltages from the near end and far end compartments

	# Set file names to read
	nearVmPath<-file.path(dirPath,"test-rallpack-near-vm.txt");
	farVmPath<-file.path(dirPath,"test-rallpack-far-vm.txt");

	# Read files. Note that column names are picked up from header line.
	nearVm<-read.csv(nearVmPath);
	farVm<-read.csv(farVmPath);

	# Select within time limits
	nearVm<-nearVm[nearVm$t>=tlim[1] & nearVm$t<=tlim[2],,drop=FALSE];
	farVm<-farVm[farVm$t>=tlim[1] & farVm$t<=tlim[2],,drop=FALSE];

	# Do a plot of near and far voltages together
	ylim<-c(min(nearVm$C0001,farVm$C1000),
		max(nearVm$C0001,farVm$C1000));
	plot(nearVm$t,nearVm$C0001,ylim=ylim,
		type='l',col='blue',
		xlab="Time (ms)",ylab="Vm (mV)",
		main="Rallpack test results");
	lines(farVm$t,farVm$C1000);
}

plot_rallpack_states<-function(
	tlim=c(0,Inf),	# time limits to plot (in ms)
	drive="C:",		# drive (or file system) containing data	
	dirPath=file.path(drive,"BNSF_1_0","Data","Default")) # data directory
{
	# Plot states from the near and far end compartments

	# Set file names to read
	nearStatesPath<-file.path(dirPath,"test-rallpack-near-states.txt");
	farStatesPath<-file.path(dirPath,"test-rallpack-far-states.txt");

	# Read files. Note that column names are picked up from header line.
	nearStates<-read.csv(nearStatesPath);
	farStates<-read.csv(farStatesPath);

	# Select within time limits
	nearStates<-nearStates[nearStates$t>=tlim[1] & nearStates$t<=tlim[2],,drop=FALSE];
	farStates<-farStates[farStates$t>=tlim[1] & farStates$t<=tlim[2],,drop=FALSE];

	# Plot both sets of states in one screen
	layout(matrix(c(1,2),2,1));

	ns<-ncol(nearStates);
	matplot(nearStates$t,nearStates[,4:ns],type='l',
		xlab="Time (ms)",ylab="State value",
		main="Near compartment states");

	ns<-ncol(farStates);
	matplot(farStates$t,nearStates[,4:ns],type='l',
		xlab="Time (ms)",ylab="State value",
		main="Far compartment states");

	# Reset layout
	layout(matrix(1,1,1))
}

