// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_300.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Test synaptic plasticity for different types of synapses.
// This file defines a class used for plasticity tests.


#include <iostream>
#include "bnsf.h"
#include "bnsf_liaf.h"
#include "synapse_glu_baker_2003.h"
#include "synapse_gaba_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


// Define classes for a neuron with synaptic conductances
class Test300Neuron : public SSRMNeuron {
public:
	Test300Neuron() : SSRMNeuron(NULL,false) { initialize(); }

	// Ensure that implicit methods are tested
	// since those are the ones used in the larger model.
	ODESolver* defaultSolver() { return new NeuronSolver; }

protected:

	virtual void addSynapticConductances()
	{		
		// Set up excitatory and inhibitory conductances
		soma()->ampa(new AMPA_SynapticResp(0));
		soma()->gaba(new GABAa_SynapticResp(0));
	}
};

