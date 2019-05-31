// K-M Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_m_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// See header file for references.
//
// Reversal potential is in the mid range of various values derived
// experimentally for potassium channels.
//
// Muscarinic modulation is to be added.


#include "ionchan_k_m_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


// ==============================================
// K_M_channel class bodies
// ==============================================

AlphaBetaEntry* K_M_channel::_abTable = NULL;

const Number K_M_channel::_Vrev = -78*mV; // From Halliwell & Adams


// Special tau computation
Number K_M_channel::tauForTable(Number v)
{

	// Get parameters in shorter notation
	Number	tmax = tauMax();
	Number	tmin = tauMin();
	Number	qt = Q10Factor();
	Number	vh = Vhalf();
	Number	ka = kalpha();
	Number	kb = kbeta();

	// Get voltage (delta from vh) at which maximum tau is reached
	// and also associated rate to use in adjusting to hit tauMax.
	Number	dvm = log(ka/kb)*ka*kb/(ka+kb);
	Number	rmax = exp(dvm/ka)+exp(-dvm/kb);

	// Compute tau = 1/(alpha(v)+beta(v)) where alpha and beta are
	// for time constant only.
	Number	a = exp((v-vh)/ka);
	Number	b = exp(-(v-vh)/kb);

	return rmax/qt*(tmax-tmin)/(a+b)+tmin;
}	
