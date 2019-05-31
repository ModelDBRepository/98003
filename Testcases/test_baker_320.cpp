// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_320.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Test synaptic plasticity for different types of synapses.
// This tests frequency dependencies in long term plasticity.

#include "test_baker_300.h"
#include "plasticity_glu_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_320() 
{
	cout<<"Test case 320 - Frequency dependent plasticity"<<endl;

	// Test parameters (change values as needed)
	SimTime				deltat = -10*msec;
	Number				freq = 50*Hz;
	Number				ach = 0*microM;

	SimTime				tinj,t;
	Number				Iinj;
	int					k,n,nspike;

	Controller			cont;
	SSRMNeuron			preNeuron;
	Test300Neuron		postNeuron;
	AxonProcess*		axon = preNeuron.axonProcess();
	SSRMSoma*			soma = postNeuron.soma();

	NR2A_SynapticResp	dummyNMDAR;	// needed to get ndmar params

	NMDARDepPlasticityRule* rule;

	ExternalSpikeRecorder	spikeRecorder ("test-baker-spikes.txt");

	// Set a temperature to report on
	IonChannel::defaultTempC(37);

	// Disable presynaptic plasticity for this test
	soma->ampa()->presynapticRule( NULL );

	// Set the postsynaptic rule to test
	soma->ampa()->postsynapticRule(rule=new CA3_AC_STDPRule(&dummyNMDAR)); 

	// Set the AChLevel for the synapse (do this after rules set)
	dummyNMDAR.AChLevel(ach);
	soma->ampa()->setModParam("AChModulator",ach);

	postNeuron.numericIdentifier(1);
	preNeuron.numericIdentifier(2);

	preNeuron.addProbe(&spikeRecorder);
	postNeuron.addProbe(&spikeRecorder);

	preNeuron.addToController(&cont);
	postNeuron.addToController(&cont);

	cont.start();

	Synapse* syn = soma->ampa()->createSynapse(axon,1.0f);

	Iinj = 0.5*nanoA;
	tinj = 0.5*msec;
	nspike = 5*20;

	cout<<"ACh (microM) = "<<ach/microM<<endl;
	cout<<"Delta T (msec) = "<<deltat/msec<<endl;
	cout<<"Frequency (Hz) = "<<freq<<endl;
	cout<<"Starting synapse weight = "<<soma->ampa()->synapseWeight(syn)<<endl;

	for (t=0,n=0; n<nspike;t+=3*sec) {

		for (k=0;k<5;k++,n++) {
			cont.runUpTo(t+k*sec/freq);
			if (deltat>0) {
				preNeuron.soma()->Iinjected(Iinj);
				cont.runForDuration(tinj);
				preNeuron.soma()->Iinjected(0);

				cont.runUpTo(t+k*sec/freq+deltat);

				postNeuron.soma()->Iinjected(Iinj);
				cont.runForDuration(tinj);
				postNeuron.soma()->Iinjected(0);
			}
			else {
				postNeuron.soma()->Iinjected(Iinj);
				cont.runForDuration(tinj);
				postNeuron.soma()->Iinjected(0);

				cont.runUpTo(t+k*sec/freq-deltat);

				preNeuron.soma()->Iinjected(Iinj);
				cont.runForDuration(tinj);
				preNeuron.soma()->Iinjected(0);
			}
		}

		cout<<"t = "<<t/sec
			<<" (sec), synapse weight = "<<soma->ampa()->synapseWeight(syn)
			<<" target weight = "<<rule->targetWeight(syn)<<endl;

	}

	cout<<"Ending synapse weight (+ 0sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;
	cont.runForDuration(10*sec);
	cout<<"Ending synapse weight (+10sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;
	cont.runForDuration(10*sec);
	cout<<"Ending synapse weight (+20sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;
	cont.runForDuration(10*sec);
	cout<<"Ending synapse weight (+30sec) = "<<soma->ampa()->synapseWeight(syn)<<endl;


	cont.finish();
}
