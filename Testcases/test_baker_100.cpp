// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_100.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This is a test case using the pyramidal cell model.
// The objective of this test case is to measure whole
// resistance and I-V relationships.
//
// The primary aspect tested here is response to continuous
// current injections to get input resistance and time constant.

#include "bnsf.h"
#include "neuron_baker_2003.h"
#include <iostream>
#include <vector>
#include <ctime>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

double estimatedTau(double h, double v1, double v2, double v3);

void test_baker_100()
{
	cout<<"Test case 100 - measure cell input resistance"<<endl;

	Number	Iinj;
	Controller* cont = new Controller;
	CA3PyramidalCell* pyr1 = new L56aPyramidalCell;

	// Set test parameters
	IonChannel::defaultTempC(37);
	SimTime settleTime=1*sec; // 1 sec is enough but 3 are better
	Iinj = -10*picoA;
	pyr1->AChLevel(0*microM);

	time_t startTime,endTime;
	Number	v0,v1,v2,v3,v4;
	double  h1,h2,tau;

	Compartment* pcomp;
	
	pyr1->numericIdentifier(1);
	pyr1->addToController(cont);

	ExternalRecorder* voltageRecorder = new ExternalVoltageRecorder(
		"test-baker-voltages.txt");
	voltageRecorder->minInterval(1*msec); // no need for high sample rate here


	pyr1->addProbe(voltageRecorder);


	// Choose current injection site (soma or dendrite)
	pcomp = pyr1->somaComp();

	// L56a is assumed cell morphology for dendrite selection
	// Note that different spatial resolutions can be used.

	// pcomp = pyr1->dendriteCompByBranch(1781,331*micron);		// obliqueDend
	// pcomp = pyr1->dendriteCompByBranch(2102,550*micron);		// distalDend
	// pcomp = pyr1->dendriteCompByBranch(1956,251*micron);		// mediumDend
	// pcomp = pyr1->dendriteCompByBranch(1467,51*micron);		// proximalDend
	// pcomp = pyr1->dendriteCompByBranch(410,239*micron);		// basalDend

	cout<<"ACh level (microM) = "<<pyr1->AChLevel()/microM<<endl;
	cout<<"Current injection (picoA) = " <<Iinj/picoA<<endl;
	cout<<"Injection compartment = "<<pcomp->componentName()<<endl;

	// Run the simulation
	cout<<endl<<"Simulation starting"<<endl;

	cont->start();
	time(&startTime);

	// Allow time to settle to an equilibrium state
	// Can use a large time step for this simulation since
	// spikes are not relevant and only equilibrium voltages matter.
	pyr1->solver()->timeStep(5*msec);
	pyr1->solver()->debugPerformance(false);
	cont->runForDuration(settleTime);

	cout<<"Resting voltage (mV) = ";
	cout<<pcomp->Vm()/UOM::mV<<endl;
	v0 = pcomp->Vm();

	// Start current injection
	pcomp->Iinjected(Iinj);

	// Use smaller time step to get time constant
	pyr1->solver()->timeStep(1*msec);
	
	h1=5*msec;
	cont->runForDuration(h1);
	v1 = pcomp->Vm();
	cont->runForDuration(h1);	
	v2 = pcomp->Vm();
	cont->runForDuration(h1);	
	v3 = pcomp->Vm();

	tau = estimatedTau(h1,v1,v2,v3);
	cout<<"Early time constant h (ms) = "<<h1/msec
		<<" tau (ms) = "<<tau/msec<<endl;

	h2=20*msec;
	cont->runForDuration(h2-2*h1);
	v2 = pcomp->Vm();
	cont->runForDuration(h2);	
	v3 = pcomp->Vm();

	tau = estimatedTau(h2,v1,v2,v3);
	cout<<"Late time constant h (ms) = "<<h2/msec
		<<" tau (ms) = "<<tau/msec<<endl;

	// Get whole cell resistance after Ih sag settles out
	// This corresponds with typical experimental protocols
	// even though full equilibrium may not be reached.
	cont->runForDuration(1000*msec-h1-2*h2);	
	v4 = pcomp->Vm();
	pcomp->Iinjected(0*picoA);

	cout<<"Voltage change = "<< (v4-v0)/mV <<" mV"<<endl;
	cout<<"Input resistance  (megaohm) = "
		<<(v4-v0)/Iinj/megaohm <<endl
		<<endl;

	cont->runForDuration(100*msec);

	time(&endTime);

	cout<<"Time to execute = ";
	cout<< (endTime-startTime);
	cout<<" sec";
	cout<<endl;

	cout<<"Simulated time = "<<(cont->evalTime()/UOM::msec);
	cout<<" msec"<<endl;

	cout<<"Derivative evaluations = "<<pyr1->solver()->derivativeEvals()<<endl;
	cout<<"Time steps done = "<<pyr1->solver()->timeStepsDone()<<endl;

	cout<<"Ending voltage (mV) = ";
	cout<<pcomp->Vm()/UOM::mV<<endl;

	cont->finish();

	// Do deletes one at a time for debugging
	delete pyr1;
	delete voltageRecorder;
	delete cont;

	cout<<"Done with deletes"<<endl;
}

// Estimate a time constant based on voltage change.
// Model is v(t1+h)=v_inf+(v1-v_inf)*exp(-h/tau).
// This gives (v3-v1)/(v2-v1)=(exp(-2h/tau)-1)/(exp(-h/tau)-1).
// A matching value for tau can be found by binary search.
// t1 is chosen to avoid an initial transient.
// Ih sag is ignored for this calculation and is included
// as part of the time constant. Obviously a variety of
// alternative fitting procedures are possible.

double estimatedTau(double h, double v1, double v2, double v3)
{
	double  r,tau,tauMin,tauMax;

	r = (v3-v1)/(v2-v1);
	tauMin=1*microsec;
	tauMax=500*msec;
	do {
		tau=(tauMin+tauMax)/2;
		if ( r>(exp(-2*h/tau)-1)/(exp(-h/tau)-1) ) {
			tauMin = tau;
		}
		else {
			tauMax = tau;
		}
	} while (tauMax-tauMin>1e-4*msec);

	return tau;
}

