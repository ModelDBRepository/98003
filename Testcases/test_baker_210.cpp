// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_210.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This is a test case for the pyramidal cell model.
// Synaptic stimulation by various types of synapses are tested.

#include "bnsf.h"
#include "bnsf_liaf.h"
#include "neuron_baker_2003.h"

#include "synapse_glu_baker_2003.h"

#include <iostream>
#include <vector>
#include <ctime>

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_210()
{
	cout<<"Test case 210"<<endl;

	time_t			startTime,endTime;
	Number			v1;
	Number			weight;

	Controller*		cont = new Controller;
	CA3PyramidalCell* pyr1 = new L56aPyramidalCell;
	LIAFNeuron*		liaf = new LIAFNeuron;
	AxonProcess*	axon = liaf->axonProcess();

	Compartment*	distalDend = NULL;
	Compartment*	mediumDend = NULL;
	Compartment*	proximalDend = NULL;
	Compartment*	basalDend = NULL;
	Compartment*	obliqueDend = NULL;

	char*			synapseType;
	Synapse*		syn;

	pyr1->numericIdentifier(1);
	pyr1->addToController(cont);

	liaf->numericIdentifier(2);
	liaf->addToController(cont);

	// Make LIAF more like IAF to get a wide range of frequencies
	liaf->soma()->Rm_specific(500*kohm/cm_2);

	ExternalRecorder* voltageRecorder = new ExternalVoltageRecorder(
		"test-baker-voltages.txt");
	voltageRecorder->minInterval(0*msec);

	ExternalRecorder* currentRecorder = new ExternalCurrentRecorder(
		"test-baker-currents.txt");
	voltageRecorder->minInterval(0*msec);

	ExternalRecorder* spikeRecorder = new ExternalSpikeRecorder(
		"test-baker-spikes.txt");

	pyr1->addProbe(voltageRecorder);
	pyr1->addProbe(currentRecorder);
	liaf->addProbe(spikeRecorder);

	// Set a temperature to simulate
	IonChannel::defaultTempC(37);

	// Set ACh level for pyramidal cell
	pyr1->AChLevel(0*microM);
	
	// Locate exemplar compartments (assumes L56a cell)
	obliqueDend = pyr1->dendriteCompByBranch(1781,331*micron);
	distalDend = pyr1->dendriteCompByBranch(2102,550*micron);
	mediumDend = pyr1->dendriteCompByBranch(1956,251*micron);
	proximalDend = pyr1->dendriteCompByBranch(1467,51*micron);
	basalDend = pyr1->dendriteCompByBranch(410,239*micron);

	// Stimulate via single synaptic connnections
	weight=1.0f;

	synapseType = "AC_GluR";
	// synapseType = "PP_GluR";
	// synapseType = "MF_GluR";
	// synapseType = "GABAa";
	// synapseType = "GABAas";
	// synapseType = "GABAb";

	syn=basalDend->createSynapse(synapseType,axon,weight);
	// syn=proximalDend->createSynapse(synapseType,axon,weight);
	// syn=mediumDend->createSynapse(synapseType,axon,weight);
	// syn=obliqueDend->createSynapse(synapseType,axon,weight);
	// syn=distalDend->createSynapse(synapseType,axon,weight);

	// Disable any plasticity rule in effect for this response.
	syn->postsynapticProcess()->disablePlasticity();

	// Add a recorder to probe synapse states
	// ExternalRecorder* synapseRecorder = new ExternalSynapseRecorder(
	// 	"test-baker-synapses.txt","AC_GluR");
	// pyr1->addProbe(synapseRecorder);

	// Optionally, voltage clamp soma to measure EPSC or IPSC
	// pyr1->somaComp()->setVoltageClamp(-70*mV);	// for EPSC
	// pyr1->somaComp()->setVoltageClamp(-55*mV);	// for IPSC

	// Optionally suppress various ionic conductances
	// pyr1->setGModulator("Na", 0);
	// pyr1->setGModulator("NaP", 0);
	// pyr1->setGModulator("CaT", 0);

	// Optionally set gModulator for AMPAR or NR2A by accessing synapse. 
	// We cannot use setGModulator as above because AMPAR and NR2A are 
	// not unique ids.

	// PP_Glu_SynapticResp* ppresp= (PP_Glu_SynapticResp*) (syn->postsynapticProcess() );
	// ppresp->nmda()->gModulator(0);

	// AC_Glu_SynapticResp* acresp= (AC_Glu_SynapticResp*) (syn->postsynapticProcess() );
	// acresp->nmda()->gModulator(0);

	cout<<"Simulation starting"<<endl;
	cout<<"AChLevel (microM) = "<<pyr1->AChLevel() / microM <<endl;
	if (obliqueDend!=NULL) 
		cout<<"Oblique compartment = "<<obliqueDend->componentName()<<endl;
	if (distalDend!=NULL) 
		cout<<"Distal compartment = "<<distalDend->componentName()<<endl;
	if (mediumDend!=NULL) 
		cout<<"Medium compartment = "<<mediumDend->componentName()<<endl;
	if (proximalDend!=NULL) 
		cout<<"Proximal compartment = "<<proximalDend->componentName()<<endl;
	if (basalDend!=NULL) 
		cout<<"Basal compartment = "<<basalDend->componentName()<<endl;
	cout<<endl;

	cont->start();
	time(&startTime);

	// Allow time to settle to an equilibrium.
	// Use a longer time step for this phase
	// When voltage clamp is used, settling is slow
	pyr1->solver()->debugPerformance(false);
	pyr1->solver()->timeStep(5*msec);

	// cont->runForDuration(10*sec); // get really, really settled
	cont->runForDuration(1*sec); // just get calmed down some

	cout<<"Settling period ended"<<endl;

	time(&endTime);
	cout<<"Execution time so far = ";
	cout<< (endTime-startTime);
	cout<<" sec";
	cout<<endl;

	v1=pyr1->somaComp()->Vm();
	cout<<"Soma voltage (mV) at rest = ";
	cout<<v1/UOM::mV<<endl;

	// Reset time step to allow tracking of voltages and currents
	pyr1->solver()->errTol(0.0001f);
	pyr1->solver()->timeStep(250*microsec);	// 4 kHz sample rate

	// Force spikes to trigger synaptic responses
	liaf->soma()->Iinjected(10*picoA); // 10 pA gives about 40Hz
	cont->runForDuration(25*msec);

	// Stop triggering input spikes
	liaf->soma()->Iinjected(0*picoA);

	// Optionally, create a spike after the synaptic input
	if (0) {
		cont->runForDuration(1*msec);
		pyr1->somaComp()->Iinjected(800*picoA);
		cont->runForDuration(3*msec);
		pyr1->somaComp()->Iinjected(0*picoA);
	}

	cont->runForDuration(100*msec);	
	// cont->runForDuration(10*sec); // useful for GABA-b	

	time(&endTime);

	cout<<"Time to execute = ";
	cout<< (endTime-startTime);
	cout<<" sec";
	cout<<endl;

	cout<<"Simulated time = "<<(cont->evalTime()/UOM::msec);
	cout<<" msec"<<endl;

	cout<<"Derivative evaluations = "<<pyr1->solver()->derivativeEvals()<<endl;
	cout<<"Time steps done = "<<pyr1->solver()->timeStepsDone()<<endl;

	cout<<"Ending soma voltage (mV) = ";
	cout<<pyr1->somaComp()->Vm()/UOM::mV<<endl;

	cont->finish();

	delete pyr1;
	delete liaf;
	delete voltageRecorder;
	delete currentRecorder;
	delete spikeRecorder;
	delete cont;

	cout<<"Done with deletes"<<endl;
}

