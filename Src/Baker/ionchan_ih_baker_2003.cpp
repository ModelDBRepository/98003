// Ih Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_ih_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// a K-h channel modeled on data from Magee (1998).
//
// Properties modeled are for the case in which [Na] is typical of
// physiological values. Proximal parameters are an extrapolation from
// the measurements where [Na]=0. Deactivation time constants are increased
// as per description by Magee for the [Na]>0 case.
//
// References:
// 
// Maccaferri G, Mangoni M, Lazzari A, Difrancesco D (1993)
// Properties of the hyperpolarization-activated current in rat
// CA1 hippocampal pyramidal cells. J Neurophysiol. 69, 2129-2136.
// 
// Magee JC, 1998. Dendritic hyperpolarization-activated currents modify
// the integrative properties of hippocampal CA1 pyramidal neurons.
// J. Neuroscience 18(19), 7613-7624.


#include "ionchan_ih_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


// ==============================================
// Ih class bodies
// ==============================================

// ----------------------------------------------
// K_h_channel bodies
// ----------------------------------------------


AlphaBetaEntry* Proximal_Ih_channel::_abTable = NULL;
AlphaBetaEntry* Distal_Ih_channel::_abTable = NULL;


// Reversal potential is estimated using GHK with Na/K perm = .35
const Number Proximal_Ih_channel::_Vrev = -25*mV;
const Number Distal_Ih_channel::_Vrev = -25*mV;
const Number Blended_Ih_channel::_Vrev = -25*mV;


// Constructors ----------------------------

Proximal_Ih_channel::Proximal_Ih_channel(Number gSpVal)
{
	gSpecific(gSpVal);
}

Distal_Ih_channel::Distal_Ih_channel(Number gSpVal)
{
	gSpecific(gSpVal);
}

Blended_Ih_channel::Blended_Ih_channel(Number gSpVal, Number bratio)
: Order1BlendedIonChannel(
	gSpVal,
	new Proximal_Ih_channel,
	new Distal_Ih_channel,
	bratio) {}
