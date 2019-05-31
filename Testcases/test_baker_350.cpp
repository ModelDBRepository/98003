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
// Test presynaptic plasticity for GLU synapses.

#include "test_baker_300.h"
#include "plasticity_glu_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_350() 
{
	cout<<"Test case 350 - Presynaptic Glu plasticity"<<endl;

	Number				ach=0*microM;
	Number				freq = 20*Hz;
	int					nspike = 10;
	char				synapseType='A'; // A=AC synapse, M=MF synapse

	SimTime				tinj,t;
	Number				Iinj;
	int					n;

	Controller			cont;
	SSRMNeuron			preNeuron;
	
	Test300Neuron		postNeuron;

	AxonProcess*		axon = preNeuron.axonProcess();
	SSRMSoma*			soma = postNeuron.soma();

	AMPA_PresynapticRule* rule;

	// Disable any postsynaptic plasticity
	soma->ampa()->postsynapticRule( NULL ); 

	// Set the rule to test (AC or MF)
	if (synapseType=='A') {
		soma->ampa()->presynapticRule(rule = new CA3_AC_PairedPulseRule);
	}
	else if (synapseType=='M') {
		soma->ampa()->presynapticRule(rule = new CA3_MF_PairedPulseRule);
	}
	else {
		cerr<<"Unknown synapse type"<<endl;
		exit(1);
	}

	// Set the AChLevel for the synapse (must do this after rules set)
	soma->ampa()->setModParam("AChModulator",ach);

	postNeuron.numericIdentifier(1);
	preNeuron.numericIdentifier(2);

	preNeuron.addToController(&cont);
	postNeuron.addToController(&cont);

	cont.start();

	Synapse* syn = soma->ampa()->createSynapse(axon);
	
	Iinj = 2*nanoA;
	tinj = 0.5*msec;

	cout<<"ACh (microM) = "<<ach/microM<<endl;
	cout<<"Frequency (Hz) = "<<freq<<endl;
	cout<<"Number of spikes = "<<nspike<<endl;
	cout<<"Synapse type = "<<synapseType<<endl;

	for (t=0,n=0;n<nspike;t+=1/freq,n++) {
		cont.runUpTo(t);

		cout<<" n = "<<n
			<<" rel prob or quantity = "<<rule->releaseProbability(syn)
			<<endl;

			preNeuron.soma()->Iinjected(Iinj);
			cont.runForDuration(tinj);
			preNeuron.soma()->Iinjected(0);
	}

	cont.finish();
}
