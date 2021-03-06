// Rallpack2 test implementation
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: test_rallpack2.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This implements the rallpack2 test case.
// There is a small error arising from the injection point
// being in the center of the first compartment.

#include "bnsf.h"
#include "neuron_rallpack.h"
#include <iostream>
#include <ctime>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace RALLPACK;


void test_rallpack2()
{
	cout<<"Rallpack2 - passive branching cable"<<endl;

	Number	Iinj;
	Controller* cont = new Controller;
	RPNeuron2* cell = new RPNeuron2;

	// Set test parameters
	Iinj = 0.1*nanoA;

	time_t startTime,endTime;

	ExternalRecorder* nearRecorder = new ExternalVoltageRecorder(
		"test-rallpack-near-vm.txt");
	ExternalRecorder* farRecorder = new ExternalVoltageRecorder(
		"test-rallpack-far-vm.txt");

	// Connect cell with recorders and controller
	cell->numericIdentifier(1);
	cell->nearComp()->addProbe(nearRecorder);
	cell->farComp()->addProbe(farRecorder);
	cell->addToController(cont);

	// Run the simulation
	cout<<endl<<"Simulation starting"<<endl;
	cont->start();
	time(&startTime);

	// Start current injection
	cell->nearComp()->Iinjected(Iinj);

	// Set a time step size to compare with genesis results (or not)
	// cell->solver()->timeStep(50*microsec);
	cell->solver()->timeStep(250*microsec);
	// cell->solver()->timeStep(1*msec);
	cell->solver()->errTol(0.001f);
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

	cout<<"Derivative evaluations = "<<cell->solver()->derivativeEvals()<<endl;
	cout<<"Time steps done = "<<cell->solver()->timeStepsDone()<<endl;

	cont->finish();

	delete cell;
	delete nearRecorder;
	delete farRecorder;
	delete cont;

	cout<<"Done with deletes"<<endl;
}
