// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_020.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Print channel alpha-beta tables to verify proper results
// The output files can be graphed with MATLAB or similar tools.
// See plot_alpha_beta.m

#include "bnsf.h"
#include "ionchan_ca_baker_2003.h"
#include "ionchan_k_c_baker_2003.h"
#include "ionchan_k_ahp_baker_2003.h"
#include <iostream>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_020()
{

	// Set a temperature to report on
	IonChannel::defaultTempC(37);

	// Ca channel gates
	Ca_T_m_gate			cat_mg;
	Ca_T_h_gate			cat_hg;
	Ca_N_m_gate			can_mg;
	Ca_N_h_gate			can_hg;
	Ca_L_channel		cal_m;

	Number		v,Vghk,Gghk;
	Number		caxIn,caxOut,temp;

	// Print GHK voltage conversion tables
	if (1) {
		caxIn=cat_mg.defaultPeakCaXin();
		caxOut=cat_mg.defaultCaXout();
		temp=cat_mg.defaultTempC();

		cout<<"GHK voltage conversions for CaXin="<<caxIn/nanoM<<" nM"<<endl;
		for (v=-150;v<150;v+=10) {
			Vghk=IonChannel::ghkCaEffectivePotential(v*mV,caxIn,caxOut,temp);
			Gghk=IonChannel::ghkCaEffectiveCond(v*mV,caxIn,caxOut,temp);
			cout<<"v="<<v<<"(mV)  \tVghk="<<Vghk/mV<<"(mV)  \tGghk ="<<Gghk<<endl;
		}
	}


	char				kcfile[256];
	Number				a,b;
	EulerMethodSolver	kca_solver;	
	Model				kca_model;
	Compartment			kca_comp;

	// The following are deleted automatically
	SimpleCalciumPool*	kca_pool	= new SimpleCalciumPool;
	K_C_channel*		kc_c		= new K_C_channel(1*mS/cm_2,kca_pool);
	K_AHP_channel*		kahp_q		= new K_AHP_channel(1*mS/cm_2,kca_pool);

	// Print alpha beta values for Ca dependent channels
	cat_mg.printAlphaBetaTable("test-baker-cat-m.txt");
	cat_hg.printAlphaBetaTable("test-baker-cat-h.txt");
	can_mg.printAlphaBetaTable("test-baker-can-m.txt");
	can_hg.printAlphaBetaTable("test-baker-can-h.txt");
	cal_m.printAlphaBetaTable("test-baker-cal-s.txt");

	// OK it's the hard way, but connect up with a model
	// and solver just enough to get things to point where
	// all initializations are done.

	kca_model.solver(&kca_solver);
	kca_comp.model(&kca_model);
	kca_comp.add(kca_pool);
	kca_comp.add(kc_c);
	kca_comp.add(kahp_q);

	kca_solver.start();

	cout<<endl;
	for (caxIn=0*nanoM;caxIn<=2001*nanoM;caxIn+=100*nanoM) {
		kca_pool->stateValue(0)=caxIn;

		// Dump K-C alpha-beta values to a file
		sprintf(kcfile,"test-baker-kc-c-%d.txt",int(0.5+caxIn/nanoM));
		FILE* kcf = fopen(kcfile,"w");
		for (v=-150;v<=150;v+=.25) {
			kca_comp.Vm(v*mV);
			fprintf(kcf,"%g,%g,%g\n",v*mV,kc_c->alpha(),kc_c->beta());
		}
		fclose(kcf);

		a=kahp_q->alpha();
		b=kahp_q->beta();
		cout<<"CaX="<<int(0.5+caxIn/nanoM)<<" nanoM\t";
		cout<<"ahp xinf="<<a/(a+b)<<" tau="<<1/msec/(a+b)<<" ms";
		cout<<endl;
	}
	cout<<endl;
	cout<<"Calcium related alpha-beta tables written"<<endl;

	kca_solver.finish();
}
