// Sodium Ion Channel
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: ionchan_na_rallpack.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the sodium (Na) channel definitions for the Rallpack benchmark.


#include "ionchan_na_rallpack.h"

using namespace std;
using namespace BNSF;

using namespace RALLPACK;
using namespace UOM;


// ==============================================
// Na channel class bodies
// ==============================================

// ----------------------------------------------
// m gate body
// ----------------------------------------------

AlphaBetaEntry* Na_m_gate::_abTable = NULL;

Number Na_m_gate::alphaForTable(Number v)
{
	return linoidRate(0.1/(msec*mV),v+40*mV,10*mV);
}

Number Na_m_gate::betaForTable(Number v)
{
	return 4.0/msec * exp(-(v+65*mV)/(18*mV));
}


// ----------------------------------------------
// h gate body
// ----------------------------------------------

AlphaBetaEntry* Na_h_gate::_abTable = NULL;

Number Na_h_gate::alphaForTable(Number v)
{
	return 0.07/msec * exp(-(v+65*mV)/(20*mV));
}

Number Na_h_gate::betaForTable(Number v)
{
	return 1.0/msec / (1+exp(-(v+35*mV)/(10*mV)));
}

// ----------------------------------------------
// Na_channel body
// ----------------------------------------------

const Number Na_channel::_Vrev = 50*UOM::mV;

Na_channel::Na_channel(Number gSpVal)
: M3HIonChannel(gSpVal) 
{
	add(new Na_m_gate );
	add(new Na_h_gate );
}

Na_channel::~Na_channel() {}
