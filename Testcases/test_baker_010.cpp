// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_010.cpp
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
#include "ionchan_k_a_baker_2003.h"
#include "ionchan_na_baker_2003.h"
#include "ionchan_k_dr_baker_2003.h"
#include "ionchan_k_m_baker_2003.h"
#include "ionchan_ih_baker_2003.h"

#include <iostream>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_010()
{
	// Set a temperature to report on
	IonChannel::defaultTempC(37);
	
	// Na channel gates
	Proximal_Na_m_gate		na_pmg;
	Distal_Na_m_gate		na_dmg;
	Blended_Na_m_gate		na_bmg;

	Soma_Na_s_gate			na_ssg;
	Proximal_Na_s_gate		na_psg;
	Distal_Na_s_gate		na_dsg;
	Blended_Na_s_gate		na_bsg;

	Soma_Na_m_gate			na_smg;
	Axon_Na_m_gate			na_amg;
	Axon_Na_h_gate			na_ahg;

	Persistent_Na_m_gate	nap_mg;
	Persistent_Na_h_gate	nap_hg;

	Na_h_gate		na_hg;

	na_pmg.printAlphaBetaTable("test-baker-prox-na-m.txt");
	na_dmg.printAlphaBetaTable("test-baker-dist-na-m.txt");
	na_bmg.printAlphaBeta("test-baker-blend-na-m.txt");

	na_ssg.printAlphaBetaTable("test-baker-soma-na-s.txt");
	na_psg.printAlphaBetaTable("test-baker-prox-na-s.txt");
	na_dsg.printAlphaBetaTable("test-baker-dist-na-s.txt");
	na_bsg.printAlphaBeta("test-baker-blend-na-s.txt");

	na_smg.printAlphaBetaTable("test-baker-soma-na-m.txt");
	na_amg.printAlphaBetaTable("test-baker-axon-na-m.txt");
	na_ahg.printAlphaBetaTable("test-baker-axon-na-h.txt");

	nap_mg.printAlphaBetaTable("test-baker-nap-m.txt");
	nap_hg.printAlphaBetaTable("test-baker-nap-h.txt");

	na_hg.printAlphaBetaTable("test-baker-na-h.txt");

	// K-A channel gates
	Proximal_K_A_a_gate		ka_pag;
	Distal_K_A_a_gate		ka_dag;
	K_A_b_gate				ka_bg;

	ka_pag.printAlphaBetaTable("test-baker-prox-ka-a.txt");
	ka_dag.printAlphaBetaTable("test-baker-dist-ka-a.txt");
	ka_dag.printAlphaBeta("test-baker-blend-ka-a.txt");
	ka_bg.printAlphaBeta("test-baker-ka-b.txt");

	// K-DR  and K-M channels
	K_DR_channel			kdr_chan;
	kdr_chan.printAlphaBetaTable("test-baker-kdr-n.txt");

	K_M_channel				km_chan;
	km_chan.printAlphaBetaTable("test-baker-km.txt");

	// Ih channel
	Proximal_Ih_channel		ih_pchan;
	Distal_Ih_channel		ih_dchan;
	Blended_Ih_channel		ih_bchan;

	ih_pchan.printAlphaBetaTable("test-baker-prox-ih.txt");
	ih_dchan.printAlphaBetaTable("test-baker-dist-ih.txt");
	ih_bchan.printAlphaBeta("test-baker-blend-ih.txt");

	cout<<"Alpha-beta tables written"<<endl;

}

