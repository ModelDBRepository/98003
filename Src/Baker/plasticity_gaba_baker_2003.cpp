// Synaptic Plasticity for GABAergic Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: plasticity_gaba_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file contains the classes used to implement
// GABA synapse plasticity. Two version were considered: one from
// data involving regular spike trains (Scanziani) and the
// other from single pulse pairings (Otis et al.). Of these,
// the spike train data is probably more representative.
// 
// References: see header file

#include "plasticity_gaba_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


// ================================================================
// GABA Presynaptic Rule Class for paired pulse depression
// ================================================================


// Static values
const TokenId GABA_PresynapticRule::_GABA_PresynapticId 
	= token("GABA_PresynapticRule");

// Set event quantity based on the presynaptic rule (Scanziani, 2000)
void GABA_PresynapticRule::finalizeAPEvent(ActionPotentialEvent* apEvent)
{
	// Handle initial spike explicitly since exp
	// does unhelpful things with infinite values.
	// Since this is a one-time event, we assume
	// release occurs unconditionally just this once.
	if (apEvent->isi()==InfiniteFuture) {
		apEvent->quantity(1);
	}

	// Otherwise, compute the ongoing PPF or PPD value
	else {
		Number	isi = apEvent->isi();
		Number	ppfd;

		// Get the paired pulse adjustment and set the event quantity
		ppfd = 1+a0()*qdexp(-isi/tau());
		ppfd *= 1+AChA()*AChLevel()/( AChLevel() + AChKd() );

		// Set quantity to either 0 or 1 if random release or pp if not.
		if (useRandomRelease() ) {
			apEvent->quantity( synapticCond()->runif()<ppfd ? 1 : 0 );
		}
		else {
			apEvent->quantity(ppfd);
		}
	}
	
	// Indicate that this event has been handled
	apEvent->isFinal(true);
}

// Apply an ACh neuromodulation rule
void GABA_PresynapticRule::setModParams(TokenId id, int nv, Number* values)
{
	using namespace UOM;

	static const TokenId AChMod = token("AChModulator");

	// Skip any other forms of modulation and check num of params
	if (id!=AChMod) return;
	if (nv<1) {
		FatalError("(GABA_PresynapticRule::setModParams) "
			"Too few parameter values supplied");
	}

	// Set the new concentration value
	AChLevel( values[0] );
}


// ================================================================
// GABAb Postsynaptic Rule Class
// ================================================================


// Constructors and destructor
GABAb_PostsynapticRule::GABAb_PostsynapticRule(
	Number baseWght, 
	Number amparInc,
	SimTime tau)
{
	// Save parameters
	_baseWeight = baseWght;
	_amparIncrement = amparInc;
	_tauW = tau;

	// Initialize other values
	_amparInitialized = false;
	_weight = _baseWeight;
}

GABAb_PostsynapticRule::~GABAb_PostsynapticRule() {}


// Locate associated AMPAR
void GABAb_PostsynapticRule::locateAMPAR()
{
	static const int		numIds = 2;
	static const TokenId	respToFind[numIds] = {
		token("AC_GluR"), 
		token("PP_GluR") 
	};

	int					k;
	Compartment*		comp = synapticCond()->container();		

	// Discard any AMPAR previously located
	if (_amparInitialized) {
		_ampar.resize(0);
	}

	// Try to find each synaptic response type and save the address.
	// There is no danger of any of these going away unexpectedly.
	for (k=0; k<numIds; k++) {

		SynapticResponse*	resp;

		// Locate the associated response and if found, remember it
		resp = comp->findSynapticResponse(respToFind[k],false);
		if (resp!=NULL) {
			_ampar.push_back(resp);
		}
	}
}


// Take action at the end of the time step (before AP purge)
// This should be done after the glutamate plasticity updates are done,
// but weights change slowly enough it really does not matter.
void GABAb_PostsynapticRule::applyEndOfStep(ActionPotentialEventQueue& apQueue)
{
	SynapticResponseVectorIt	it;
	SimTime						h = timeStepSize();
	SimTime						tau = tauW();
	Number						totalAMPARWeight=0;
	Number						targetWeight;
	int							synapseCount=0;

	// Make sure AMPAR have been located in this compartment.
	// This is postponed until first use to simplify start-up.
	if (!_amparInitialized) {
		locateAMPAR();
		_amparInitialized=true;
	}

	// Get an average synaptic weight by going through each 
	// associated AMPAR response type and accumulating.
	for (it=_ampar.begin(); it!=_ampar.end(); it++) {

		SynapticResponse* resp = *it;
		int n;

		totalAMPARWeight += resp->totalSynapticWeight(&n);
		synapseCount += n;
	}

	// Now get the resulting GABAb synapse weight. This value
	// is accessed by the associated response and used to adjust
	// conductance.
	targetWeight = baseWeight();
	if (synapseCount>0) {
		targetWeight += totalAMPARWeight/synapseCount * amparIncrement();
	}

	// Move current weight towards target weight using a time
	// constant of tau. A rough approximation is used assuming
	// that tauW is much larger than the current time step.
	// Otherwise, qdexp might be appropriate here.
	_weight += h/(tau+h/2)*(targetWeight - _weight);
}

