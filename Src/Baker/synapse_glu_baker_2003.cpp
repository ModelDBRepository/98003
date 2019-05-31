// Synapse Dynamics for Glutamate Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: synapse_glu_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// glutamate synapse conductances.
//
// See header file for references.

#include "synapse_glu_baker_2003.h"
#include "plasticity_glu_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;



// ================================================================
// AMPA synaptic conductance static values
// ================================================================

const Number AMPA_SynapticResp::_Vrev = 0*mV;
DualExpSynapticCondClassCache AMPA_SynapticResp::_CC;

// Constructors and destructor
AMPA_SynapticResp ::AMPA_SynapticResp (Number gVal) 
: DualExpSynapticCond(gVal) {}

AMPA_SynapticResp ::~AMPA_SynapticResp () {}

// =================================================
// NMDA Synaptic Conductance Class Body
// =================================================


// Static variables
const Number NMDA_SynapticResp::_Vrev = 0*mV;
NMDAVmTableEntry* NMDA_SynapticResp::_NMDAVmTable = NULL;

// Constructors and destructor
NMDA_SynapticResp::NMDA_SynapticResp(Number gVal) 
: DualExpSynapticCond(gVal) {}

NMDA_SynapticResp::~NMDA_SynapticResp() {}

// Load the Mg gate table when the simulation starts
void NMDA_SynapticResp::simulationStarted()
{
	// Let superclass do any of its initializations
	DualExpSynapticCond::simulationStarted();

	// Load the effective potential table on the first call
	if (*pNMDAVmTable()==NULL)
		loadNMDAVmTable();
}

// Delete the effective potential table when the simulation ends
void NMDA_SynapticResp::simulationEnded()
{
	// Delete the table
	delete[] *pNMDAVmTable();
	*pNMDAVmTable() = NULL;

	// Let superclass do its own cleanup
	DualExpSynapticCond::simulationEnded();
}

// Load the effective potential array
void NMDA_SynapticResp::loadNMDAVmTable()
{
	NMDAVmTableEntry*	tbl;

	int				k;
	Number			v;

	Number			CaXin = defaultPeakCaXin();
	Number			CaXout = defaultCaXout();

	Number			dVm = 1e-2*VStepForIndex;

	// Delete any previously allocated table
	if (*pNMDAVmTable()!=NULL) {
		delete[] *pNMDAVmTable();
	}

	// Allocate a new table and remember where it is
	*pNMDAVmTable() = tbl = new NMDAVmTableEntry[VTableSize];

	// Go through the voltages one step at a time
	for (k=0,v=VMinForIndex; k<VTableSize; k++, v+=VStepForIndex ) {

		// Get the (instantaneous) Mg++ gate value
		tbl[k].MgGate = MgGateValueForTable(v);
	}
}

// Set ACh concentration level and modulation
void NMDA_SynapticResp::AChLevel(Number ach)
{
	// Save param value
	_AChLevel = ach;

	// Get a new modulation value
	_AChMod  = 1+AChA1()*ach/( ach+AChKd1() );
	_AChMod *= 1+AChA2()*ach/( ach+AChKd2() );

	// Set the new conductance modulation value
	gModulator( _AChMod );
}

// Apply an ACh neuromodulation rule
void NMDA_SynapticResp::setModParams(TokenId id, int nv, Number* values)
{
	using namespace UOM;

	static const TokenId AChMod = token("AChModulator");

	// Skip any other forms of modulation and check num of params
	if (id!=AChMod) return;
	if (nv<1) {
		FatalError("(NMDA_SynapticResp::setModParams::setModParams) "
			"Too few params");
	}

	// Set the new concentration value
	AChLevel( values[0] );
}

// Get conductance including the Mg++ gate value and Ca++ effective conductance
Number NMDA_SynapticResp::conductance()
{
	NMDAVmTableEntry*	ent = *pNMDAVmTable()+container()->VmIndex();

	Number				mgGate = VTableInterp(container()->VmRem()/VStepForIndex,
								(ent-1)->MgGate,ent->MgGate,
								(ent+1)->MgGate,(ent+2)->MgGate);

	return	g()*s2()*mgGate; 
}

// Get current flow assuming Ohm's law behavior
Number NMDA_SynapticResp::Iion()
{
	NMDAVmTableEntry*	ent = *pNMDAVmTable()+container()->VmIndex();

	Number				mgGate = VTableInterp(container()->VmRem()/VStepForIndex,
								(ent-1)->MgGate,ent->MgGate,
								(ent+1)->MgGate,(ent+2)->MgGate);

	// As a fastpath, include conductance calculation inline.
	return	g()*s2()*mgGate*(Vm()-Vrev()); 
}

// Get the decay time constant based on rated Vm value.
SimTime NMDA_SynapticResp::tau2()
{
	Number tauAtV0 = ratedTau2() /
		maxval(nhdMin(),nhdAtV0()+nhdSlope()*ratedVm());
	Number tauAtVm = tauAtV0 *
		maxval(nhdMin(),nhdAtV0()+nhdSlope()*nominalVm());
	return tauAtVm;
}

// Return the Mg++ plug state for a given voltage
Number NMDA_SynapticResp::MgGateValueForTable(Number vm) 
{
	return 1/(1+MgExtConc()/MgKdAtV0()*exp(-MgVmMult()*vm));
}


// ================================================================
// Combined NR2A and NR2B static values
// ================================================================


DualExpSynapticCondClassCache NR2A_SynapticResp::_CC;
DualExpSynapticCondClassCache NR2B_SynapticResp::_CC;


// ================================================================
// AC_Glu_SynapticResp class body
// ================================================================


// Constructor
AC_Glu_SynapticResp::AC_Glu_SynapticResp (
	Number gAMPAR,		// AMPAR conductance 
	Number gNMDAR )		// NMDAR conductance
{
	add(_ampa = new AMPA_SynapticResp(gAMPAR));
	add(_nmda = new NR2A_SynapticResp(gNMDAR));

	// Set the presynaptic rule
	ampa()->presynapticRule(new CA3_AC_PairedPulseRule);

	// Set the postsynaptic rule passing associated NMDAR
	ampa()->postsynapticRule(new CA3_AC_STDPRule(nmda()) );
}


// ================================================================
// PP_Glu_SynapticResp class body
// ================================================================


// Constructor
PP_Glu_SynapticResp::PP_Glu_SynapticResp (
	Number gAMPAR,		// AMPAR conductance 
	Number gNMDAR )		// NMDAR conductance
{
	add(_ampa = new AMPA_SynapticResp(gAMPAR));
	add(_nmda = new NR2A_SynapticResp(gNMDAR));

	// Presynaptic properties of PP are not well
	// determined. MEC-PP shows PPD while LEC-PP
	// shows PPF. For simplicity, a random
	// release rule is used. 
	ampa()->presynapticRule(new RandomReleaseRule);

	// Set the postsynaptic rule. There is no
	// evidence for this in PP synapses, but it
	// it is the best available data and something
	// like this probably applies for MEC PP.
	// LEC PP is more of an unknown and opiod
	// dependencies have been found experimentally.
	ampa()->postsynapticRule(new CA3_PP_STDPRule(nmda()) );
}


// ================================================================
// MF_AMPA_SynapticResp class body
// ================================================================


const Number MF_AMPA_SynapticResp::_Vrev = 0*mV;
DualExpSynapticCondClassCache MF_AMPA_SynapticResp::_CC;

MF_AMPA_SynapticResp::MF_AMPA_SynapticResp (Number gVal) 
{ 
	// Set the peak achievable conductance per synapse
	gMax(gVal); 

	// Set the presynaptic rule for MF PPF
	presynapticRule(new CA3_MF_PairedPulseRule);
}

MF_AMPA_SynapticResp::~MF_AMPA_SynapticResp () {}


// ================================================================
// MF_Glu_SynapticResp class body
// ================================================================


// Constructor
MF_Glu_SynapticResp::MF_Glu_SynapticResp (
	Number gAMPAR,		// AMPAR conductance 
	Number gNMDAR )		// NMDAR conductance
{
	add(_ampa = new MF_AMPA_SynapticResp(gAMPAR));
	add(_nmda = new NR2A_SynapticResp(gNMDAR));
}


