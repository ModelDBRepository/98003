// K-C Ion Channel Dynamics
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_c_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K-C channel. Parameters are adapted from Moczydlowski 
// and Latorre (1982). A similar model is also distributed 
// as an example with Neuron. See also Migliori et al. 1995.
//
// See header file for references.


#include "ionchan_k_c_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;

using namespace BAKER_2003;

// ==============================================
// K_C_channel body
// ==============================================

const Number K_C_channel::_Vrev = -91*mV; 

K_C_channel::K_C_channel(Number gScaled, CalciumPool* pool) 
: HHIonChannel(gScaled)
{
	calciumPool(pool);
}

K_C_channel::~K_C_channel() {}

Number K_C_channel::alpha()
{
	// Parameters from Moczydlowski and Latorre rescaled
	// to fit calcium and voltage sensitivity from 
	// Gong et al. at 2 microM in adult CA1 pyramidal cells.
	
	static const Number calciumScale = 16;
	static const Number voltageScale = 1.95f;							

	static const Number a		= 480/sec; // beta in M&L
	static const Number z		= 2;
	static const Number zd1		= 0.84f * z / voltageScale;
	static const Number k1		= 1.8e-4*molar / calciumScale;
	
	Number c = CaX();

	// Note that this is 1/tau-closed in M&L
	// Rewrite the expression in M&L to allow c=0
	return a*c/(c+k1*exp(-zd1*Vm()*FoverRT() )) * Q10Factor();
}

Number K_C_channel::beta()
{
	// Parameters from Moczydlowski and Latorre rescaled
	// to fit calcium and voltage sensitivity from 
	// Gong et al. at 2 microM in adult CA1 pyramidal cells.

	static const Number calciumScale = 16;
	static const Number voltageScale = 1.95f;							

	static const Number b		= 280/sec;	// alpha in M&L
	static const Number z		= 2;
	static const Number zd2		= 1.0f * z / voltageScale;
	static const Number k2		= 1.1e-5*molar / calciumScale;

	Number c = CaX();

	// Note that this is 1/tau-open in M&L
	return b/(1+c/k2*exp(zd2*Vm()*FoverRT() )) * Q10Factor();
}

Number K_C_channel::conductance()
{
	return g()*value();
}
