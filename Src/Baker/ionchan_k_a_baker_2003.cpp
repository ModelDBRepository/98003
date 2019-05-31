// K-A Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_a_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K-A channel based on data from the Hoffman 1997 article.
//
// The K+ reversal potential is chosen to fit behavior. See Bekkers
// for a measurement of -82mV in neocortical cells. Martina et al.
// claimed a measurement of -96mV for CA1. It seems that Hoffman et al.
// used -80mV. though the article is not specific. Klee et al. appear
// to have assumed (or possibly measured) a value of -78mV.
// If we assume a Pna/Pk ratio of 0.01 (Aidley & Stanfield table 5.1),
// then -85mV is a reasonable reversal potential for K+ channels. 
//
// References:
//
// Aidley DJ and Stanfield PR, 1996. Ion Channels: molecules in action.
// New York: Cambridge University Press.
//
// Bekkers JM, 2000. Properties of voltage-gated potassium currents in
// nucleated patches from large laryer 5 cortical pyramidal neurons
// of the rat. J. Physiology 525.3, 593-609.
//
// Borg-Graham LJ, 1998. Interpretations of Data and Mechanisms for 
// Hippocampal Cell Models, in Cerebral Cortex vol 13. New York: Plenum Press.
// Also available via Surf-Hippo web site.
//
// Hoffman DA, Magee JC, Colbert CM, Johnston D, 1997.
// K+ channel regulation of signal propagation in dendrites
// of hippocampal cells. Nature 387(6636), 869-875.
//
// Klee R, Ficker E, Heinemann U, 1995. Comparison of voltage-dependent 
// potassium currents in rat pyramidal neurons acutely isolated from 
// hippocampal regions CA1 and CA3.
//
// Martina M, Schultz JH, Ehmke H, Monyer H, Jonas P 1998. Functional
// and molecular differences between voltage-gated K+ channels of
// fast-spiking interneourns and pyramidal neurons of rat hippocampus.
//
// Pan E and Colbert CM, 2001. Subthreshold incactivation of Na+ and
// K+ channels supports activity-dependent enhancement of back-propagating
// action potentials in hippocampal CA1. J. Neurophysiology 85, 1013-1016.


#include "ionchan_k_a_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;

using namespace BAKER_2003;


// ==============================================
// K_A class bodies
// ==============================================

// ----------------------------------------------
// Gate bodies
// ----------------------------------------------

AlphaBetaEntry* Proximal_K_A_a_gate::_abTable = NULL;
AlphaBetaEntry* Distal_K_A_a_gate::_abTable = NULL;
AlphaBetaEntry* K_A_b_gate::_abTable = NULL;

Number K_A_b_gate::tauForTable(Number v)
{
	// Data from Hoffman shows that for -25mV < v < 55mV, the
	// relationship between tau and v is tau = 0.26*(v+56) at 35 deg C. 
	// Pan showed that at 32 deg C tau is around 6 msec 
	// at -50mV and is roughly 5 msec at -60 mV. Obviously tau
	// must assume some non-zero minimum value.

	// The cutoff at the minimum tau is implemented with a linoid rate
	// function to provide a transition that is smooth near the minimum.

	const Number tauSlope	= 0.26*msec/mV;
	const Number vmin		= -56*mV+tauMin()/tauSlope;
	const Number k			= 5*mV;

	Number			tau;

	tau = linoidRate(tauSlope,v-vmin,k)/Q10Factor();
	tau += tauMin()/Q10FactorForTauMin();
	return tau;
}

/// ----------------------------------------------
// Proximal_K_A_channel body
// ----------------------------------------------

const Number Proximal_K_A_channel::_Vrev = -85*mV; 

Proximal_K_A_channel::Proximal_K_A_channel(Number gSpVal)
{
	gSpecific(gSpVal);
	add( new Proximal_K_A_a_gate );
	add( new K_A_b_gate );
}

// ----------------------------------------------
// Distal_K_A_channel body
// ----------------------------------------------

const Number Distal_K_A_channel::_Vrev = -85*mV; 

Distal_K_A_channel::Distal_K_A_channel(Number gSpVal)
{
	gSpecific(gSpVal);
	add( new Distal_K_A_a_gate);
	add( new K_A_b_gate );
}

// ----------------------------------------------
// Blended_K_A_channel body
// ----------------------------------------------

const Number Blended_K_A_channel::_Vrev = -85*mV; 

Blended_K_A_a_gate::Blended_K_A_a_gate(Number bratio) : BlendedIonGate(
	new Proximal_K_A_a_gate,
	new Distal_K_A_a_gate,
	bratio) {}

Blended_K_A_channel::Blended_K_A_channel(Number gSpVal, Number bratio)
{
	gSpecific(gSpVal);
	add( new Blended_K_A_a_gate(bratio));
	add( new K_A_b_gate );
}
