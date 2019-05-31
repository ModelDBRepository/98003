function plot_test_rallpack

% Copyright 2007 John L Baker. All rights reserved.
% This software is provided AS IS under the terms of the Open Source
% MIT License. See http://www.opensource.org/licenses/mit-license.php.

% NOTE: files and directories referenced will need to be customized
% to go with the local values used during installation and testing.
% Plots can be further refined as needed, of course.

DATADIR='C:\BNSF_1_0\Data\Default\';

B=dlmread([DATADIR,'test-rallpack-near-vm.txt'],',',1,1);
NT=B(:,1);
NV=B(:,2);

B=dlmread([DATADIR,'test-rallpack-far-vm.txt'],',',1,1);
FT=B(:,1);
FV=B(:,2);

figure
plot(NT,NV,'.-',FT,FV,'-')
title('Rallpack3 axon far end waveform')
xlabel('Time (msec)')
ylabel('Vm (mV)')
legend('Near end Vm','Far end Vm')

spidx=find(NV(1:end-1)<-30 & NV(2:end)>=-30);
num_spikes = length(spidx)

