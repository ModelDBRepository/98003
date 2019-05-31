// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_670.cpp
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		6 March 2007
//
// Description:
//
// Simulate regular mouse movement and generated inputs.
// Results are written to files for visualization.
//
// A traing DG cells supplies the initial place cell
// location. During the test, ACh is set to an initial
// level and then changed over time (1 sec intervals).
// ACh levels do not go below the final value specified.
//
// For this test the maze is rectangular and the path 
// followed is parallel with maze walls. This allows 
// place field characteristics in a linear track to 
// be simulated.
//
// A warning is in order if whole cell voltage probing
// is done. Over 1 gigabyte is generated per 5 minutes 
// of simulated time. Obviously this can add up quickly 
// for longer runs.

#include <iostream>
#include <ctime>

#include "bnsf.h"
#include "mouse_baker_2003.h"

using namespace std; 
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_670() 
{
	cout<<"Test case 670 - linear track place cells with decreasing ACh"<<endl;

	int					seedIndex=0; // random number seed index (0 to 4)

	double				trainingDuration = 9*minute;
	double				noveltyDuration = 6*minute;
	double				phaseProbeDuration = 3*minute;
	double				completionProbeDuration = 3*minute;

	Number				initialACh = 50*microM;
	Number				finalACh = 20*microM;

	Number				deltaACh = (initialACh-finalACh)/(noveltyDuration);
	Number				currentACh;

	// Selective knockout switches
	bool				koNMDARCurrent = false;
	bool				koNMDARPlasticity = false;
	bool				koNMDARAChMod = false;
	bool				koNMDARCaDepSupp = false;
	bool				koThetaPhasePrecession = false;

	// Other control switches and params
	bool				useTeacherCell = true;
	bool				disablePlasticityForProbe = true;
	bool				chgDGThetaPhaseToMatchCA3 = false;
	bool				disableTheta = false;
	bool				probeVm = false;
	bool				probeCa = false;

	bool				varySpeedForProbe = false;
	Number				secondaryProbeSpeed = 5*cm/sec;

	// Simulation objects
	Controller*			cont = new Controller;
	Maze*				maze = new RectangularMaze;
	Mouse*				mouse = new Mouse(maze, 
							Mouse::wallFollower,
							Mouse::defaultSeeds(seedIndex));
	
	PlaceCell*			pc;  // current place cell (for debug access)
	int					k;

	time_t				startTime,endTime;
	double				diffTime;
	double				t;

	time(&startTime);

	// Set up external recorders to capture results
	ExternalRecorder*	mouseStateRecorder = new ExternalStateRecorder("test-baker-mouse-states.txt");

	ExternalRecorder*	cellVoltageRecorder = new ExternalVoltageRecorder("test-baker-voltages.txt");
	ExternalRecorder*	cellCalciumRecorder = new ExternalCalciumRecorder("test-baker-calcium.txt");

	ExternalRecorder*	cellSpikeRecorder = new ExternalSpikeRecorder("test-baker-spikes.txt");
	ExternalRecorder*	dendSpikeRecorder = new ExternalDendriteSpikeRecorder("test-baker-dendrite-spikes.txt");
	ExternalRecorder*	ECSpikeRecorder = new ExternalSpikeRecorder("test-baker-ec-spikes.txt");
	ExternalRecorder*	DGSpikeRecorder = new ExternalSpikeRecorder("test-baker-dg-spikes.txt");
	ExternalRecorder*	CA3SpikeRecorder = new ExternalSpikeRecorder("test-baker-ca3-spikes.txt");
	ExternalRecorder*	INSpikeRecorder = new ExternalSpikeRecorder("test-baker-in-spikes.txt");
	ExternalRecorder*	subsetSpikeRecorder = new ExternalSpikeRecorder("test-baker-subset-spikes.txt");

	ExternalRecorder*	cellACSynapseRecorder = 
		new ExternalSynapseRecorder("test-baker-ac-synapse-weights.txt","AC_GluR");
	cellACSynapseRecorder->minInterval(10*sec);
	
	ExternalRecorder*	cellPPSynapseRecorder = 
		new ExternalSynapseRecorder("test-baker-pp-synapse-weights.txt","PP_GluR");
	cellPPSynapseRecorder->minInterval(10*sec);

	// Show statistics
	mouse->printSynapseCounts();

	// Save place cell centers
	mouse->ECPlaceCells()->printPlaceFieldCenters("test-baker-ec-pcloc.txt");
	mouse->DGPlaceCells()->printPlaceFieldCenters("test-baker-dg-pcloc.txt");
	mouse->CA3PlaceCells()->printPlaceFieldCenters("test-baker-ca3-pcloc.txt");

	// Hook up the mouse
	mouse->addToController(cont);
	mouse->addProbe(mouseStateRecorder);

	// Add a spike recorder probe to all layers
	mouse->ECPlaceCells()->addProbeToAll(ECSpikeRecorder);
	mouse->DGPlaceCells()->addProbeToAll(DGSpikeRecorder);
	mouse->CA3PlaceCells()->addProbeToAll(CA3SpikeRecorder);

	mouse->CA3AxoAxonicCells()->addProbeToAll(INSpikeRecorder);
	mouse->CA3BasketCells()->addProbeToAll(INSpikeRecorder);
	mouse->CA3BistratifiedCells()->addProbeToAll(INSpikeRecorder);
	mouse->CA3OLMCells()->addProbeToAll(INSpikeRecorder);

	// Debug ODE solver performance (or not)
	mouse->targetCell()->solver()->debugPerformance(false);

	// Add probes to target cell
	mouse->targetCell()->addProbe(cellSpikeRecorder);
	mouse->targetCell()->addProbe(dendSpikeRecorder);
	mouse->targetCell()->addProbe(subsetSpikeRecorder);
	mouse->targetCell()->addProbe(cellACSynapseRecorder);
	mouse->targetCell()->addProbe(cellPPSynapseRecorder);

	// Optional probe of voltage and calcium
	if (probeVm) {
		mouse->targetCell()->addProbe(cellVoltageRecorder);
	}
	if (probeCa) {
		mouse->targetCell()->addProbe(cellCalciumRecorder);
	}

	// Add probes to selected cells for ease of plotting
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

	// Perform any selective NMDAR knockouts
	mouse->targetCell()->knockoutNMDAR(
		koNMDARCurrent, 
		koNMDARPlasticity,
		koNMDARAChMod, 
		koNMDARCaDepSupp );

	// If requested, disable theta phase precession
	if (koThetaPhasePrecession) {
		mouse->ECPlaceCells()->disableThetaPhasePrecession();
		mouse->DGPlaceCells()->disableThetaPhasePrecession();
		mouse->CA3PlaceCells()->disableThetaPhasePrecession();
	}

	// If requested, change DG theta phase to match CA3
	if (chgDGThetaPhaseToMatchCA3) {

		Number CA3ThetaPhase = mouse->CA3PlaceCells()->placeCell(1)->thetaPhase();

		for (k=1;k<=mouse->DGPlaceCells()->numCells();k++) {
			pc=mouse->DGPlaceCells()->placeCell(k);
			pc->thetaPhase(CA3ThetaPhase);
		}
	}

	// If requested, set theta amplitude to zero suppressing theta.
	if (disableTheta) {
		cout<<"Setting theta amplitudes to zero"<<endl;
		mouse->ECPlaceCells()->thetaAmplitude(0);
		mouse->DGPlaceCells()->thetaAmplitude(0);
		mouse->CA3PlaceCells()->thetaAmplitude(0);
		mouse->CA3AxoAxonicCells()->thetaAmplitude(0);
		mouse->CA3BasketCells()->thetaAmplitude(0);
		mouse->CA3BistratifiedCells()->thetaAmplitude(0);
		mouse->CA3OLMCells()->thetaAmplitude(0);
	}

	// Progress report
	time(&endTime);
	diffTime = difftime(endTime,startTime);

	cout<<endl;
	cout<<"Time to build layers (sec) = "<<diffTime<<endl;
	cout<<"Simulated trainging duration (sec) = "<<trainingDuration/sec<<endl;
	cout<<"Simulated phase probe duration (sec) = "<<phaseProbeDuration/sec<<endl;
	cout<<"Simulated completion probe duration (sec) = "<<completionProbeDuration/sec<<endl;
	cout<<"Knockout flags = "
		<< koNMDARCurrent 
		<< koNMDARPlasticity
		<< koNMDARAChMod
		<< koNMDARCaDepSupp
		<< koThetaPhasePrecession
		<< endl;
	cout<<"Vary speed for a 2nd phase probe pass = "<<varySpeedForProbe<<endl;
	cout<<"Random number seed index = "<<seedIndex<<endl;

	// Run the simulation
	cout<<endl;
	cout<<"Simulation starting"<<endl;
	cont->start();

	if (useTeacherCell) {
		// Pick one DG cell as the training input and make it active
		// Start by making possible afferent DG cells inactive.
		mouse->DGPlaceCells()->setActivity(false,1,50);	
		mouse->DGPlaceCells()->placeCell(21)->setCenter(0*cm,-17.5*cm);
		mouse->DGPlaceCells()->placeCell(21)->isActive( true );
	}

	// Set training pass ACh Level
	currentACh = initialACh;
	mouse->targetCell()->AChLevel(currentACh);
	cout<<"Initial ACh Level (microM) = "<<mouse->targetCell()->AChLevel()/microM<<endl;
	cout<<"Final ACh Level (microM) = "<<finalACh/microM<<endl;

	// Run the training interval simulation
	cout<<"time (sec)=";
	for (t=0;t<trainingDuration;t+=1*sec) {

		// Run for 1 sec
		cout<<t<<" ";
		cont->runForDuration(1*sec);

		// Set new ACh level
		currentACh = maxval(currentACh-deltaACh,finalACh);
		mouse->targetCell()->AChLevel(currentACh);
		if (int(t)%10==9) {
			cout<<" ACh="<<mouse->targetCell()->AChLevel()/microM<<endl;
		}
	}
	cout<<endl;

	// Run the probe pass simulations using current ACh levels
	cout<<"Simulation continuing - phase probe pass"<<endl;
	
	// Perform any selective NMDAR knockouts.
	mouse->targetCell()->knockoutNMDAR(
		koNMDARCurrent, 
		disablePlasticityForProbe,
		koNMDARAChMod, 
		koNMDARCaDepSupp );

	// Run for more time
	cout<<"time (sec)=";
	for (t=0;t<phaseProbeDuration;t+=1*sec) {
		cout<<t<<" ";
		if (int(t)%10==9) cout<<endl;
		cont->runForDuration(1*sec);
	}
	cout<<endl;

	if (varySpeedForProbe) {

		Number oldSpeed=mouse->speed();

		// Set a new speed
		mouse->speed(secondaryProbeSpeed);

		// Run the probe pass simulations using current ACh levels
		cout<<"Phase probe continuing - movement speed = "
			<<secondaryProbeSpeed / (cm/sec)<<" cm/sec"
			<<endl;

		// Run for more time
		cout<<"time (sec)=";
		for (t=0;t<phaseProbeDuration;t+=1*sec) {
			cout<<t<<" ";
			if (int(t)%10==9) cout<<endl;
			cont->runForDuration(1*sec);
		}
		cout<<endl;

		// Restore the original speed
		mouse->speed(oldSpeed);
	}


	// Run the probe pass simulations using current ACh levels
	cout<<"Simulation continuing - pattern completion probe pass"<<endl;
	cout<<"Plasticity disabled for probe = " <<disablePlasticityForProbe<<endl;
	
	// Perform any selective NMDAR knockouts.
	mouse->targetCell()->knockoutNMDAR(
		koNMDARCurrent, 
		disablePlasticityForProbe,
		koNMDARAChMod, 
		koNMDARCaDepSupp );

	// Make the training inputs inactive for the probe
	mouse->DGPlaceCells()->setActivity(false,1,50);	

	// Run for more time
	cout<<"time (sec)=";
	for (t=0;t<completionProbeDuration;t+=1*sec) {
		cout<<t<<" ";
		if (int(t)%10==9) cout<<endl;
		cont->runForDuration(1*sec);
	}
	cout<<endl;

	// Wrap up
	cont->runForDuration(1*msec);	// ensure final synapse reports
	cout<<"Simulation ending"<<endl;

	cout<<"Derivative evaluations = "
		<<mouse->targetCell()->solver()->derivativeEvals()<<endl;

	cout<<"Time steps done = "
		<<mouse->targetCell()->solver()->timeStepsDone()<<endl;

	cont->finish();

	time(&endTime);
	diffTime = difftime(endTime,startTime);
	cout<<"Time to execute (sec) = "<<diffTime<<endl;

	// Delete allocated objects. Doing it this way
	// allows debugging if something goes wrong.
	cout<<"Deleting data structures"<<endl;
	delete mouse;
	delete cont;
	delete maze;
	delete mouseStateRecorder;
	delete cellVoltageRecorder;
	delete cellCalciumRecorder;
	delete cellSpikeRecorder;
	delete dendSpikeRecorder;
	delete ECSpikeRecorder;
	delete DGSpikeRecorder;
	delete CA3SpikeRecorder;
	delete INSpikeRecorder;
	delete subsetSpikeRecorder;
	delete cellACSynapseRecorder;
	delete cellPPSynapseRecorder;
	cout<<"Done"<<endl;
}
