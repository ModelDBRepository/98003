// Na Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_na_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement variations
// on the sodium (Na) channel. See header file for references and notes.


#include "ionchan_na_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;

using namespace BAKER_2003;


// ==============================================
// Na channel class bodies
// ==============================================


// ----------------------------------------------
// Class static values
// ----------------------------------------------

AlphaBetaEntry* Na_h_gate::_abTable = NULL;

AlphaBetaEntry* Soma_Na_s_gate::_abTable = NULL;
AlphaBetaEntry* Proximal_Na_s_gate::_abTable = NULL;
AlphaBetaEntry* Distal_Na_s_gate::_abTable = NULL;

AlphaBetaEntry* Proximal_Na_m_gate::_abTable = NULL;
AlphaBetaEntry* Distal_Na_m_gate::_abTable = NULL;

AlphaBetaEntry* Soma_Na_m_gate::_abTable = NULL;

AlphaBetaEntry* Axon_Na_m_gate::_abTable = NULL;
AlphaBetaEntry* Axon_Na_h_gate::_abTable = NULL;

AlphaBetaEntry* Persistent_Na_m_gate::_abTable = NULL;
AlphaBetaEntry* Persistent_Na_h_gate::_abTable = NULL;

// Reversal potentials from Gasparini & Magee
// Blended reversal is average of distal and proximal.
const Number Proximal_Na_channel::_Vrev = 55.8*mV; 
const Number Distal_Na_channel::_Vrev = 54.4*mV; 
const Number Blended_Na_channel::_Vrev = 54.9*mV; 

// Reversal pontential for other Na channels
const Number Soma_Na_channel::_Vrev = 55*mV; 
const Number Axon_Na_channel::_Vrev = 55*mV; 
const Number Persistent_Na_channel::_Vrev = 55*mV;


// ----------------------------------------------
// Abstract_Na_h_gate body
// ----------------------------------------------


// Constructor
Abstract_Na_h_gate::Abstract_Na_h_gate(IonChannel* mg)
{
	_mgate = mg;

	// Set initial parameter values.
	// Subclasses can override via constructor
	// or wait for simulationStarted event below.
	_inactRate = 0;
	_tauZero = 0;
}

Abstract_Na_h_gate::~Abstract_Na_h_gate() {}

// Save values for fast path code
void Abstract_Na_h_gate::simulationStarted()
{
	// Set parameters using subclass overrides if any.
	_inactRate = inactRate();
	_tauZero = tauZero();

	// Let superclass take over the rest of the logic
	EnergyBarrierTabGate::simulationStarted();
}

// Fast path computation of derivatives.
void Abstract_Na_h_gate::computeDerivatives()
{
	// Get alpha and beta values
	Number a=alpha();
	Number b=beta();

	// Get the current activation state from the mgate
	Number m = mgate()->value();

	// Now get the adjusted xinf and tau values
	Number tau = 1/(a+b+_inactRate*m*m*m);
	Number xinf = tau*a;

	// Apply xinf and tau to get derivative
	// loadAlphaBetaTable prevents tau=0 case from occurring
	derivValue(0) = (xinf-stateValue(0))/(tau+_tauZero);
}

// Fast path computation of local state update via implicit rule.
void Abstract_Na_h_gate::localStateUpdate(SimTime h, CNStepType stepType)
{
	// Get alpha and beta values
	Number a=alpha();
	Number b=beta();

	// Get the current activation state from the mgate
	Number m = mgate()->value();

	// Now get the adjusted xinf and tau values
	Number tau = 1/(a+b+_inactRate*m*m*m);
	Number xinf = tau*a;
	Number x = stateValue(0);

	// Perform the update (semi-implicit trapezoid rule)
	if (stepType==CNStartingHalfStep) {
		stateValue(0) += h/(tau+_tauZero+h)*(xinf-x);
	}
	else {
		stateValue(0) += h/(tau+_tauZero+h/2)*(xinf-x);
	}
}


// ----------------------------------------------
// Abstract_Na_s_gate body
// ----------------------------------------------


// Constructor and destructor
Abstract_Na_s_gate::Abstract_Na_s_gate(bool isDisabled, Number dval)
{
	using namespace UOM;

	// Initialize values
	_gateDisabled = isDisabled;
	_disabledValue = dval;

	// Set a default value (value is not critical)
	_disabledTau = 2*msec;
}

Abstract_Na_s_gate::~Abstract_Na_s_gate() {}

// Override compute derivatives to test gate status.
// Fall back to the Hodgkin-Huxley alpha-beta formulas.
void Abstract_Na_s_gate::computeDerivatives()
{
	HHIonChannel::computeDerivatives();
}
	
// Override localStateUpdate to handle gate status
// Fall back to the Hodgkin-Huxley alpha-beta formulas.
void Abstract_Na_s_gate::localStateUpdate(SimTime h, CNStepType stepType)
{
	HHIonChannel::localStateUpdate(h, stepType);
}

// Return an HH alpha value, adjusted if needed by gate status
Number Abstract_Na_s_gate::alpha()
{
	return gateDisabled()
		? disabledValue()/disabledTau()
		: EnergyBarrierTabGate::alpha();
}

// Return an HH beta value, adjusted if needed by gate status
Number Abstract_Na_s_gate::beta()
{
	return gateDisabled()
		? (1-disabledValue())/disabledTau()
		: EnergyBarrierTabGate::beta();
}


// ----------------------------------------------
// Blended_Na_s_gate body
// ----------------------------------------------


Blended_Na_s_gate::Blended_Na_s_gate(
	Number bratio,					// blend ratio
	bool isDisabled,				// set gateDisabled on or off
	Number dval)					// disabled value 
: BlendedIonGate(
	new Proximal_Na_s_gate,
	new Distal_Na_s_gate,
	bratio)
{
	_gateDisabled = isDisabled;
	_disabledValue = dval;
}

Blended_Na_s_gate::~Blended_Na_s_gate() {}

void Blended_Na_s_gate::gateDisabled(bool flag)
{
	// Pass-thru to underlying gates
	sgate1()->gateDisabled(flag);
	sgate2()->gateDisabled(flag);
}

void Blended_Na_s_gate::disabledValue(Number s)
{
	// Pass-thru to underlying gates
	sgate1()->disabledValue(s);
	sgate2()->disabledValue(s);
}


// ----------------------------------------------
// Proximal_Na_channel body
// ----------------------------------------------


Proximal_Na_channel::Proximal_Na_channel(
	Number		gSpVal,					// specific conductance
	bool		sGateDisabled,			// s gate disabled flag
	Number		sGateValue)				// s gate value if disabled
{
	IonChannel* mgate;

	gSpecific(gSpVal);
	add( mgate = new Proximal_Na_m_gate );
	add( new Na_h_gate(mgate) );
	add( _sGate = new Proximal_Na_s_gate(sGateDisabled) );
}


// ----------------------------------------------
// Distal_Na_channel body
// ----------------------------------------------


Distal_Na_channel::Distal_Na_channel(
	Number		gSpVal,					// specific conductance
	bool		sGateDisabled,			// s gate disabled flag
	Number		sGateValue)				// s gate value if disabled
{
	IonChannel* mgate;

	gSpecific(gSpVal);
	add( mgate = new Distal_Na_m_gate );
	add( new Na_h_gate(mgate) );
	add( _sGate = new Distal_Na_s_gate(sGateDisabled,sGateValue) );
}


// ----------------------------------------------
// Blended_Na_channel (Proximal blended with Distal)
// ----------------------------------------------

		
Blended_Na_channel::Blended_Na_channel(
	Number		gSpVal,					// specific conductance
	Number		bratio,					// blend ratio
	bool		sGateDisabled,			// s gate disabled flag
	Number		sGateValue)				// s gate value if disabled
{
	IonChannel* mgate;

	gSpecific(gSpVal);
	add( mgate = new Blended_Na_m_gate(bratio));
	add( new Na_h_gate(mgate));
	add( _blendedSGate = new Blended_Na_s_gate(bratio,sGateDisabled,sGateValue) );
}


// ----------------------------------------------
// Soma_Na_channel body
// ----------------------------------------------


Soma_Na_channel::Soma_Na_channel(
	Number		gSpVal,					// specific conductance
	bool		sGateDisabled,			// s gate disabled flag
	Number		sGateValue)				// s gate value if disabled
{
	IonChannel* mgate;

	gSpecific(gSpVal);
	add( mgate = new Soma_Na_m_gate );
	add( new Na_h_gate(mgate) );
	add( _sGate = new Proximal_Na_s_gate(sGateDisabled, sGateValue) );
}


// ----------------------------------------------
// Axon_Na_channel body
// ----------------------------------------------


Axon_Na_channel::Axon_Na_channel(Number gSpVal)
{
	gSpecific(gSpVal);
	add( new Axon_Na_m_gate );
	add( new Axon_Na_h_gate );
}


// ----------------------------------------------
// Persistent_Na_channel body
// ----------------------------------------------


Persistent_Na_channel::Persistent_Na_channel(Number gSpVal)
{
	gSpecific(gSpVal);
	add( new Persistent_Na_m_gate );
	add( new Persistent_Na_h_gate );
}

