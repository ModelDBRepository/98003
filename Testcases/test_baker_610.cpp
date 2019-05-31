// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_610.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Simulate random mouse movement and generated inputs.
// Results are written to files for visualization.
// The network is not used to speed up the simulation.

#include <iostream>
#include <ctime>

#include "bnsf.h"
#include "mouse_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_610() 
{
	cout<<"Test case 610 - mouse movement and input generation"<<endl;


	double				duration = 1*minute;	// how long to simulate
	int					seedIndex=0;			// which random seeds to use
	bool				useCircularMaze=false;	// use circular vs square maze
	bool				useNovelArea=false;		// use a novelty region

	Controller*			cont = new Controller;

	Maze*				mwm;
	Mouse*				mouse; 

	SpatialRegion*		novelArea = new RectangularRegion(-15*cm,-10*cm,5*cm,10*cm);

	Number				feat1[] = {1,0,0};
	Number				feat2[] = {0,1,0};
	Number				feat3[] = {0,0,1};

	PlaceCell*			pc;
	int					k;

	time_t				startTime,endTime;
	double				diffTime;

	double				t;

	time(&startTime);

	// Create the maze and mouse objects without a network hookup
	if (useCircularMaze) {
		mwm = new CircularMaze;
		mouse = new Mouse(mwm, Mouse::randomWalk, 
					Mouse::defaultSeeds(seedIndex), false);
	}
	else {
		mwm = new RectangularMaze;
		mouse = new Mouse(mwm, Mouse::wallFollower,
					Mouse::defaultSeeds(seedIndex), false);
	}

	// Define recorders
	ExternalRecorder*	stateRecorder = new ExternalStateRecorder("test-baker-mouse-states.txt");
	ExternalRecorder*	ECSpikeRecorder = new ExternalSpikeRecorder("test-baker-ec-spikes.txt");
	ExternalRecorder*	DGSpikeRecorder = new ExternalSpikeRecorder("test-baker-dg-spikes.txt");
	ExternalRecorder*	CA3SpikeRecorder = new ExternalSpikeRecorder("test-baker-ca3-spikes.txt");
	ExternalRecorder*	INSpikeRecorder = new ExternalSpikeRecorder("test-baker-in-spikes.txt");
	ExternalRecorder*	subsetSpikeRecorder = new ExternalSpikeRecorder("test-baker-subset-spikes.txt");

	// Hook up the mouse for simulation
	mouse->addToController(cont);
	mouse->addProbe(stateRecorder);

	// Save place cell centers
	cout<<"Saving place field centers"<<endl;
	mouse->ECPlaceCells()->printPlaceFieldCenters("test-baker-ec-pcloc.txt");
	mouse->DGPlaceCells()->printPlaceFieldCenters("test-baker-dg-pcloc.txt");
	mouse->CA3PlaceCells()->printPlaceFieldCenters("test-baker-ca3-pcloc.txt");

	// Add the spike recorder probe to all layers
	mouse->ECPlaceCells()->addProbeToAll(ECSpikeRecorder);
	mouse->DGPlaceCells()->addProbeToAll(DGSpikeRecorder);
	mouse->CA3PlaceCells()->addProbeToAll(CA3SpikeRecorder);

	mouse->CA3AxoAxonicCells()->addProbeToAll(INSpikeRecorder);
	mouse->CA3BasketCells()->addProbeToAll(INSpikeRecorder);
	mouse->CA3BistratifiedCells()->addProbeToAll(INSpikeRecorder);
	mouse->CA3OLMCells()->addProbeToAll(INSpikeRecorder);

	// Add selected cells to subset recorder for ease of plotting
	for (k=1;k<=200;k++) {
		pc=mouse->ECPlaceCells()->placeCell(k);
		pc->addProbe(subsetSpikeRecorder);
	}
	for (k=1;k<=200;k++) {
		pc=mouse->DGPlaceCells()->placeCell(k);
		pc->addProbe(subsetSpikeRecorder);
	}
	for (k=1;k<=200;k++) {
		pc=mouse->CA3PlaceCells()->placeCell(k);
		pc->addProbe(subsetSpikeRecorder);
	}

	// Add some landmarks at different locations
	// These are not really used at present, so this
	// is just a test of the interface.
	mwm->add(new Landmark(  0, 50, 1, 3, feat1) );
	mwm->add(new Landmark( 50,  0, 2, 3, feat2) );
	mwm->add(new Landmark(-50,  0, 3, 3, feat3) );

	// Optionally, specify novelty region
	if (useNovelArea) {
		mouse->noveltyRegion(novelArea);
	}

	// Pick one DG cell as the training input and make it active
	// Start by making possible afferent DG cells inactive.
	mouse->DGPlaceCells()->setActivity(false,1,50);	
	mouse->DGPlaceCells()->placeCell(21)->setCenter(0*cm,-17.5*cm);
	mouse->DGPlaceCells()->placeCell(21)->isActive( true );

	// Progress report
	time(&endTime);
	diffTime = difftime(endTime,startTime);

	cout<<endl;
	cout<<"Time to build layers (sec) = "<<diffTime<<endl;
	cout<<"Simulated duration (sec) = "<<duration/sec<<endl;

	// Kick things off
	cout<<"Simulation starting"<<endl;
	cont->start();
	
	// Run the simulation
	for (t=0;t<duration;t+=1*sec) {
		cout<<"time (sec)="<<t<<endl;
		cont->runForDuration(1*sec);
	}

	cout<<"Simulation ending"<<endl;
	cont->finish();

	time(&endTime);
	diffTime = difftime(endTime,startTime);
	cout<<"Time to execute (sec) = "<<diffTime<<endl;

	cout<<"Deleting data structures"<<endl;
	delete mouse;
	delete cont;
	delete mwm;
	delete novelArea;
	delete stateRecorder;
	delete ECSpikeRecorder;
	delete DGSpikeRecorder;
	delete CA3SpikeRecorder;
	delete INSpikeRecorder;
	delete subsetSpikeRecorder;
	cout<<"Done"<<endl;
}
