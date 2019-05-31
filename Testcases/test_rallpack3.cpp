// Rallpack3 test implementation
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: test_rallpack3.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This implements the rallpack3 test case.

#include "bnsf.h"
#include "neuron_rallpack.h"
#include <iostream>
#include <ctime>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace RALLPACK;

void test_rallpack3()
{
	cout<<"Rallpack3 - cable with Na and K channels"<<endl;

	int nComp = 1000;
	RPNeuron3* cell = new RPNeuron3(nComp);

	Number	Iinj;
	Controller* cont = new Controller;
	time_t startTime,endTime;

	ExternalRecorder* nearVoltage = new ExternalVoltageRecorder(
		"test-rallpack-near-vm.txt");
	ExternalRecorder* farVoltage = new ExternalVoltageRecorder(
		"test-rallpack-far-vm.txt");

	ExternalRecorder* nearStates = new ExternalStateRecorder(
		"test-rallpack-near-states.txt");
	ExternalRecorder* farStates = new ExternalStateRecorder(
		"test-rallpack-far-states.txt");

	// Set test parameters
	Iinj = 0.1*nanoA;

	// Connect cell with recorders and controller.
	cell->numericIdentifier(1);
	cell->nearComp()->addProbe(nearVoltage);
	cell->farComp()->addProbe(farVoltage);

	// Ion channels states are not part of Rallpack
	// as such, but may be interesting. Connecting
	// the recorder with the neuron would capture
	// the states of ion channels in all compartments,
	// but this would be of unworkable size for display.
	cell->nearComp()->addProbe(nearStates);
	cell->farComp()->addProbe(farStates);
	cell->addToController(cont);

	// Run the simulation
	cout<<endl<<"Simulation starting"<<endl;
	cont->start();
	time(&startTime);	// Start current injection
	cell->nearComp()->Iinjected(Iinj);

	// Set a time step size to compare with genesis results (or not)
	// Note that accuracy is minimally affected by step size.

	cell->solver()->timeStep(250*microsec);
	// cell->solver()->timeStep(50*microsec);
	// cell->solver()->timeStep(1*msec);
	cell->solver()->errTol(0.001f);

	// Turn off (default) or on debug output from the solver.
	cell->solver()->debugPerformance(false);
	
	// Simulate for the required duration
	cont->runForDuration(250*msec);

	time(&endTime);

	cout<<"Time to execute = ";
	cout<< (endTime-startTime);
	cout<<" sec";
	cout<<endl;

	cout<<"Simulated time = "<<(cont->evalTime()/UOM::msec);
	cout<<" msec"<<endl;

	cout<<"Number of compartments = "<<nComp<<endl;
	cout<<"Injection current (nanoA) = "<<Iinj/nanoA<<endl;

	cout<<"Derivative evaluations = "<<cell->solver()->derivativeEvals()<<endl;
	cout<<"Time steps done = "<<cell->solver()->timeStepsDone()<<endl;

	cont->finish();

	delete cell;
	delete nearVoltage;
	delete farVoltage;
	delete nearStates;
	delete farStates;
	delete cont;

	cout<<"Done with deletes"<<endl;
}
