// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_110.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
// This is a test case using the pyramidal cell model.
// This case involves current injections to the soma
// to create spike reponses. 

#include "bnsf.h"
#include "neuron_baker_2003.h"

#include <iostream>
#include <vector>
#include <ctime>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_110()
{
	cout<<"Test case 110"<<endl;

	time_t startTime,endTime;
	Controller* cont = new Controller;
	CA3PyramidalCell* pyr1 = new L56aPyramidalCell;

	// Set simulation parameters
	IonChannel::defaultTempC(37);
	SimTime settleTime=1*sec;
	pyr1->AChLevel(0*microM);

	pyr1->numericIdentifier(1);
	pyr1->addToController(cont);

	ExternalRecorder* calciumRecorder = new ExternalCalciumRecorder(
		"test-baker-calcium.txt");
	calciumRecorder->minInterval(1*msec);

	ExternalRecorder* voltageRecorder = new ExternalVoltageRecorder(
		"test-baker-voltages.txt");
	voltageRecorder->minInterval(0*msec);

	// Add probes
	pyr1->addProbe(voltageRecorder);
	pyr1->addProbe(calciumRecorder);
	
	cout<<"Simulation starting"<<endl;
	cout<<"ACh level (microM) = "<<pyr1->AChLevel()/microM<<endl;

	cont->start();
	time(&startTime);

	// Allow time to settle to an equilibrium state
	// Use a longer time step for this phase
	pyr1->solver()->timeStep(5*msec);
	pyr1->solver()->debugPerformance(false);
	cont->runForDuration(settleTime);

	time(&endTime);
	
	cout<<"Settling period ended"<<endl;
	cout<<"Resting soma voltage (mV) = ";
	cout<<pyr1->somaComp()->Vm()/UOM::mV<<endl;
	cout<<"Settling real time (sec) = "<<endTime-startTime<<endl;

	time(&startTime);

	// Allow different frequencies of sampling voltages

	// pyr1->solver()->timeStep(250*microsec);	// 4 kHz
	pyr1->solver()->timeStep(50*microsec);	// 20 kHz

	pyr1->solver()->errTol(0.001f);
	pyr1->solver()->debugPerformance(false);

	// Somatic current injections
	pyr1->somaComp()->Iinjected(1*nanoA);
	cont->runForDuration(5*msec);

	pyr1->somaComp()->Iinjected(0*picoA);
	cont->runForDuration(97*msec);

	time(&endTime);
	cout<<"Derivative evaluations = "<<pyr1->solver()->derivativeEvals()<<endl;
	cout<<"Time steps done = "<<pyr1->solver()->timeStepsDone()<<endl;

	cout<<"Ending soma voltage (mV) = ";
	cout<<pyr1->somaComp()->Vm()/UOM::mV<<endl;

	cout<<"Simulation real time (sec) = ";
	cout<< (endTime-startTime);
	cout<<endl;

	cout<<"Simulated time = "<<(cont->evalTime()/UOM::msec);
	cout<<" msec"<<endl;

	cont->finish();

	// Do deletes one at a time for debugging
	delete pyr1;
	delete cont;
	delete voltageRecorder;
	delete calciumRecorder;

	cout<<"Done with deletes"<<endl;
}

