// Calcium Ion Channel Dynamics
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the Ca++ channel definitions. Only the most relevant
// channel types are included. Ca-N channels are composite
// of medium HVA channels.
//
// See the header file for additional notes and references.

#include "ionchan_ca_baker_2003.h"

using namespace std;
using namespace BNSF;

using namespace BAKER_2003;
using namespace UOM;


// ==============================================
// Ca_T_channel class bodies
// ==============================================


AlphaBetaEntry* Ca_T_m_gate::_abTable = NULL;
AlphaBetaEntry* Ca_T_h_gate::_abTable = NULL;

Ca_T_channel::Ca_T_channel(Number gSpVal)
{
	gSpecific(gSpVal);
	add(new Ca_T_m_gate);
	add(new Ca_T_h_gate);
}


// ==============================================
// Ca_N_channel class bodies
// ==============================================


AlphaBetaEntry* Ca_N_m_gate::_abTable = NULL;
AlphaBetaEntry* Ca_N_h_gate::_abTable = NULL;

Ca_N_channel::Ca_N_channel(Number gSpVal)
{
	gSpecific(gSpVal);
	add(new Ca_N_m_gate);
	add(new Ca_N_h_gate);
}


// ==============================================
// Ca_L_channel class body
// ==============================================


AlphaBetaEntry* Ca_L_channel::_abTable = NULL;

Ca_L_channel::Ca_L_channel(Number gSpVal)
{
	gSpecific(gSpVal);
}
