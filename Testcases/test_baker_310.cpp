// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_310.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Test homosynaptic plasticity for GLU synapses.

#include "test_baker_300.h"
#include "plasticity_glu_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


void test_baker_310() 
{
	cout<<"Test case 310 - homosynaptic LTD"<<endl;

	SimTime				tmax,tinj,t;
	Number				Iinj,freq;

	Controller			cont;
	SSRMNeuron			preNeuron;
	Test300Neuron		postNeuron;
	AxonProcess*		axon = preNeuron.axonProcess();
	SSRMSoma*			soma = postNeuron.soma();

	NR2A_SynapticResp	dummyNMDAR;	// needed to get ndmar params

	NMDARDepPlasticityRule* rule;

	ExternalSpikeRecorder	spikeRecorder ("test-baker-spikes.txt");

	// Disable presynaptic plasticity for this test
	postNeuron.soma()->ampa()->presynapticRule( NULL );

	// Set the postsynaptic rule to test
	postNeuron.soma()->ampa()->postsynapticRule(rule=new CA3_AC_STDPRule(&dummyNMDAR)); 

	postNeuron.numericIdentifier(1);
	preNeuron.numericIdentifier(2);

	preNeuron.addProbe(&spikeRecorder);
	postNeuron.addProbe(&spikeRecorder);

	preNeuron.addToController(&cont);
	postNeuron.addToController(&cont);

	cont.start();

	Synapse* syn = soma->ampa()->createSynapse(axon,1.0f);

	Iinj = 1*nanoA;
	tinj = 1*msec;
	tmax = 3*minute;
	freq = 3*Hz;

	cout<<"Frequency (Hz) = "<<freq/Hz<<endl;
	cout<<"Starting synapse weight = "<<soma->ampa()->synapseWeight(syn)<<endl;

	for (t=0;t<tmax;t+=1/freq) {
		cont.runUpTo(t);

		cout<<"t = "<<t/sec
			<<" (sec), synapse weight = "<<soma->ampa()->synapseWeight(syn)
			<<" target weight = "<<rule->targetWeight(syn)<<endl;

		preNeuron.soma()->Iinjected(Iinj);
		cont.runForDuration(tinj);
		preNeuron.soma()->Iinjected(0);

	}

	cout<<"Ending synapse weight (+ 0sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;
	cont.runForDuration(10*sec);
	cout<<"Ending synapse weight (+10sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;
	cont.runForDuration(20*sec);
	cout<<"Ending synapse weight (+30sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;

	cont.finish();
}
