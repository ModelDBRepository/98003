// K Ion Channel
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: ionchan_k_rallpack.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K channel definition used in Rallpack benchmarks.



#include "ionchan_k_rallpack.h"

using namespace std;
using namespace BNSF;
using namespace UOM;

using namespace RALLPACK;


// ==============================================
// K_DR_channel body
// ==============================================


AlphaBetaEntry* K_channel::_abTable = NULL;
const Number K_channel::_Vrev = -77*UOM::mV; 

Number K_channel::alphaForTable(Number v)
{
	return linoidRate(0.01/(msec*mV),v+55*mV,10*mV);
}

Number K_channel::betaForTable(Number v)
{
	return 0.125/msec * exp(-(v+65*mV)/(80*mV));
}

Number K_channel::conductance()
{
	Number n = value();
	return g()*n*n*n*n;
}
