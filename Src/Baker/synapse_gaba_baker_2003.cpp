// Common Synapse Dynamics for GABAergic Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: synapse_gaba_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// GABAergic synaptic conductances.
//
// See header file for references.

#include "synapse_gaba_baker_2003.h"
#include "plasticity_gaba_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


// ====================================================================
// GABAa_SynapticResp Class Body
// ====================================================================


const Number GABAa_SynapticResp::_Vrev = -70*UOM::mV;
DualExpSynapticCondClassCache GABAa_SynapticResp::_CC;

// Constructors and destructor
GABAa_SynapticResp::GABAa_SynapticResp(Number gVal)
{ 
	// Set maximum conductance
	gMax(gVal);
	
	// Set the presynaptic rule.
	presynapticRule(new GABAa_PresynapticRule);
}

GABAa_SynapticResp::~GABAa_SynapticResp() {}


// ====================================================================
// GABAas_SynapticResp Class Body
// ====================================================================


const Number GABAas_SynapticResp::_Vrev = -70*UOM::mV;
DualExpSynapticCondClassCache GABAas_SynapticResp::_CC;

// Constructors and destructor
GABAas_SynapticResp::GABAas_SynapticResp(Number gVal)
{ 
	// Set maximum conductance
	gMax(gVal);
	
	// Presynaptic rules are not known. Use same rule as GABAa.
	presynapticRule(new GABAa_PresynapticRule);
}

GABAas_SynapticResp::~GABAas_SynapticResp() {}


// ====================================================================
// GABAb_SynapticCond Class Body
// ====================================================================


const Number GABAb_SynapticResp::_Vrev = -95*UOM::mV;
TripleExpSynapticCondClassCache GABAb_SynapticResp::_CC;

// Constructors and destructor
GABAb_SynapticResp::GABAb_SynapticResp(Number gVal)
{ 
	using namespace UOM;

	// Set maximum conductance
	gMax(gVal);
	
	// Set the presynaptic rule
	presynapticRule(new GABAb_PresynapticRule);

	// Set a postsynaptic rule for AMPAR-dependent weight
	// adjustments (Huang et al. 2005). Experimental results
	// only provide a rough guidance in setting parameters 
	// values. Tetanic LTP induction increased GABAb currents 
	// in synapses by a factor of four with slow expression.
	postsynapticRule(new GABAb_PostsynapticRule(1,4,60*sec)); 
}

GABAb_SynapticResp::~GABAb_SynapticResp() {}

// Return an adjusted conductance.
Number GABAb_SynapticResp::conductance()
{
	Number g=TripleExpSynapticCond::conductance();
	
	GABAb_PostsynapticRule* rule = 
		(GABAb_PostsynapticRule*) postsynapticRule();

	// Allow for the case when the postsynaptic rule is disabled
	return rule==NULL ? g : g * rule->weight();
}
