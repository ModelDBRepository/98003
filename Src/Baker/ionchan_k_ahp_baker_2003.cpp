// K-AHP Ion Channel Dynamics
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_ahp_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K-AHP channel definitions. This follows the formalism
// used in Migliore et al. but adjusted for a half activation
// at 400 nanoM [Ca++]-in (see Hirschberg et al.)
//
// Note that Lancaster & Adams suggest a temperature sensivity
// that would imply a Q10 for the current of around 16.
// Temperature effects on the calcium current may be significant. 
// A Q10 for this channel alone is not known and thus no
// Q10 is included in this implementation pending further data.
//
// References:
//
// Hirschberg B, Maulie J, Adelman JP, Marrion NV (1999). Gating
// properties of single SK channels in hippocampal CA1 pyramidal
// neurons. Biophysical Journal 77, 1905-1913.
//
// Lancaster B & Adams PR (1986). Calcium-dependent current
// generating the afterhyperpolatization of hippocampal neurons.
// J Neurophysiology 55, 1268-1282.
//
// Migliore M, Cook EP, Jaffe DB, Turner DA, and Johnston D. (1995).
// Computer simulations of morphologically reconstructed CA3
// hippocampal neurons. J. Neurophysiol. 73(3), 1157-1168.


#include "ionchan_k_ahp_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;

using namespace BAKER_2003;


// ==============================================
// K_AHP_channel class bodies
// ==============================================

const Number K_AHP_channel::_Vrev = -91*mV; 

// Constructor and destructor
K_AHP_channel::K_AHP_channel(Number gSpVal, CalciumPool* pool)
: HHIonChannel(gSpVal) 
{
	calciumPool(pool);
}

K_AHP_channel::~K_AHP_channel() {}

// q gate alpha function
Number K_AHP_channel::alpha()
{
	static const Number a0		= 1.0/(200*msec);
	static const Number CaXhalf	= 600*nanoM;

	Number x = CaX() / CaXhalf;

	return a0*x*x*x*x;
}

// q gate beta function
Number K_AHP_channel::beta()
{
	static const Number b0		= 1.0/(200*msec);
	
	return b0;
}

// compute conductance using the state variable (q)
Number K_AHP_channel::conductance() 
{
	Number q = value();
	return g()*q;
}
