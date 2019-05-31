// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_300.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Test long term synaptic plasticity for GLU synapses.
// Spike timings here are off slightly from the ideal values
// and small differences are found between results here and
// those of the MATLAB equivalent implementation of the
// plasticity algorithm.

#include "test_baker_300.h"
#include "plasticity_glu_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_300() 
{
	cout<<"Test case 300 - Paired spike plasticity"<<endl;

	Number				ach=0*microM;
	SimTime				deltat = 15*msec;
	Number				freq = 1*Hz;
	Number				w0=1.0f;

	SimTime				tinj,t;
	Number				Iinj;
	int					n;

	// For reasons not explained, Debanne et al. used
	// different numbers of pairings in LTP and LTD induction, ie
	// 60 for LTP and 100 for LTD. If replicating their results
	// this must be taken into account in setting npair.
	int					npair = 60; 

	// Set a temperature to report on
	IonChannel::defaultTempC(37);

	Controller			cont;
	SSRMNeuron			preNeuron;
	
	Test300Neuron		postNeuron;

	AxonProcess*		axon = preNeuron.axonProcess();
	SSRMSoma*			soma = postNeuron.soma();

	NR2A_SynapticResp	dummyNMDAR;	// needed to get NMDAR params

	NMDARDepPlasticityRule* rule;

	ExternalSpikeRecorder	spikeRecorder ("test-baker-spikes.txt");

	// Disable presynaptic plasticity for this test
	soma->ampa()->presynapticRule( NULL );

	// Set the postsynaptic rule to test
	soma->ampa()->postsynapticRule(rule=new CA3_AC_STDPRule(&dummyNMDAR)); 

	// Set the AChLevel for the synapse
	dummyNMDAR.AChLevel(ach);
	soma->ampa()->setModParam("AChModulator",ach);

	postNeuron.numericIdentifier(1);
	preNeuron.numericIdentifier(2);

	preNeuron.addProbe(&spikeRecorder);
	postNeuron.addProbe(&spikeRecorder);

	preNeuron.addToController(&cont);
	postNeuron.addToController(&cont);

	cont.start();

	Synapse* syn = soma->ampa()->createSynapse(axon,w0);
	
	Iinj = 2*nanoA;
	tinj = 0.5*msec;

	cout<<"ACh (microM) = "<<ach/microM<<endl;
	cout<<"Number of pairings = "<<npair<<endl;
	cout<<"Frequency (Hz) = "<<freq<<endl;
	cout<<"Pairing delta T (msec) = "<<deltat/msec<<endl;
	cout<<"Starting synapse weight = "<<soma->ampa()->synapseWeight(syn)<<endl;

	for (t=0,n=0;n<npair;t+=1/freq,n++) {
		cont.runUpTo(t);

		cout<<"t = "<<t/sec
			<<" (sec), synapse weight = "<<soma->ampa()->synapseWeight(syn)
			<<" target weight = "<<rule->targetWeight(syn)<<endl;

		if (deltat>=0) {

			preNeuron.soma()->Iinjected(Iinj);
			cont.runForDuration(tinj);
			preNeuron.soma()->Iinjected(0);

			cont.runUpTo(t+deltat);

			postNeuron.soma()->Iinjected(Iinj);
			cont.runForDuration(tinj);
			postNeuron.soma()->Iinjected(0);
		}
		else {

			postNeuron.soma()->Iinjected(Iinj);
			cont.runForDuration(tinj);
			postNeuron.soma()->Iinjected(0);

			cont.runUpTo(t-deltat);

			preNeuron.soma()->Iinjected(Iinj);
			cont.runForDuration(tinj);
			preNeuron.soma()->Iinjected(0);
		}
	}

	cont.runUpTo(t);
	cout<<"t = "<<t/sec
		<<" (sec), synapse weight = "<<soma->ampa()->synapseWeight(syn)
		<<" target weight = "<<rule->targetWeight(syn)<<endl;

	cont.runForDuration(10*sec);
	cout<<"Ending synapse weight (+10sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;
	cont.runForDuration(20*sec);
	cout<<"Ending synapse weight (+30sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;

	cont.finish();
}
