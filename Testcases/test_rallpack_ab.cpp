// Print Rallpack channel alpha and beta tables
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: test_rallpack_ab.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Print channel alpha-beta tables to verify proper results
// The output files can be graphed with MATLAB or similar tools.
// See plot_alpha_beta.m.

#include "bnsf.h"
#include "ionchan_k_rallpack.h"
#include "ionchan_na_rallpack.h"

#include <iostream>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace RALLPACK;

void test_rallpack_ab()
{
	
	// channel gates
	Na_m_gate		m;
	Na_h_gate		h;
	K_channel		n;

	m.printAlphaBetaTable("test-rallpack-m.txt");
	h.printAlphaBetaTable("test-rallpack-h.txt");
	n.printAlphaBetaTable("test-rallpack-n.txt");

	cout<<"Alpha-beta tables written"<<endl;

}

