// Simple Standard Neuron Models
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: neuron_bnsf_liaf.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// a leaky integrate and fire neuron model.
//
// See header file for references.


#include "bnsf_liaf.h"

using namespace std;
using namespace BNSF;


// ====================================================================
// SSRM Soma class body
// ====================================================================


// Create a new instance and set initial defaults
SSRMSoma::SSRMSoma()
{
	using namespace UOM;

	// Set a nominal size for the compartment
	// and set Rm,Cm to get a 30 msec time constant
	radius(10*micron);
	Rm_specific(30*kohm*cm_2);
	Cm_specific(1*microF/cm_2);

	// Set Vleak to a commonly used value
	Vleak(-65*mV);

	// Set values related to spiking
	// These are further customized by the
	// associated neuron class.
	_Vspike = +40*mV;
	_Vafter = -70*mV;
	_spikeWidth = 0;
	_refractoryInterval = 0;
	_RmRefractoryRatio = 1;

	// Initialize conductances to NULL here.
	// Actual value is provided by the neuron.
	_ampa = NULL;
	_gaba = NULL;
}

// Destroy this instance
SSRMSoma::~SSRMSoma() 
{
	// Delete allocated objects.

	// Note that this would also be cleaned up
	// by the Compartment destructor but it is
	// cleaner to take care of things explicitly.

	delete _ampa;
	delete _gaba;
}

// Set the AMPA conductance
void SSRMSoma::ampa(SynapticConductance* synCond)
{
	// Dispose of any old instance
	delete _ampa;

	// Add a new AMPA channel
	add(_ampa = synCond);
}

// Set the GABA conductance
void SSRMSoma::gaba(SynapticConductance* synCond)
{
	// Dispose of any old instance
	delete _gaba;

	// Add a new GABA channel
	add(_gaba = synCond);
}

// Add a new excitatory synapse
Synapse* SSRMSoma::addExcSynapseFrom(SpikingNeuron* aff, Number weight, Number dist)
{
	return ampa()->createSynapse(aff->axonProcess(), weight, dist);
}

// Add a new inhibitory synapse
Synapse* SSRMSoma::addInhSynapseFrom(SpikingNeuron* aff, Number weight, Number dist)
{
	return gaba()->createSynapse(aff->axonProcess(), weight, dist);
}

// Compute the leakage current. During the period immediately 
// follwing a spike, the leakage reversal potential is set to 
// Vafter and the membrane resistance is reduced.
double SSRMSoma::Ileak()
{
	if (currentTime()<neuron()->spikeTime()+spikeWidth() ) {
		return (Vm()-Vafter()) / (Rm()*RmRefractoryRatio());
	}
	else {
		return (Vm()-Vleak())/Rm();
	}
}

// Set the membrane potential immediately following a spiking event.
void SSRMSoma::setVmAfterSpike(SimTime tspike)
{
	// If there is an active voltage clamp, nothing is changed
	if (isVoltageClamped() )
		return;

	// If there is a refractory interval, set voltage to its peak
	if (refractoryInterval()>0) {

		// Set to spiking voltage at time of spike.
		// We ignore the effect of changing voltages
		// between the time of the spike and the end
		// of the time step.
		Vm(Vspike());
	}

	// Otherwise, set immediately to Vafter.
	// Similarly, change in Vm between the time
	// of the spike and end of step are ignored.
	else {
		Vm(Vafter());
	}

	// Make derivatives match the current state.
	recomputeDerivatives();
}



// ====================================================================
// SSRMNeuron class body
// ====================================================================



// Create a new instance with a default model and solver
SSRMNeuron::SSRMNeuron(bool doInit)
{
	_soma = NULL;
	if (doInit) {
		initialize();
	}
}

// Create a new instance using an existing model
SSRMNeuron::SSRMNeuron(Model* m, bool doInit) : VoltageTriggeredSpikingNeuron(m)
{
	_soma = NULL;
	if (doInit) {
		initialize();
	}
}

// Destroy this instance
SSRMNeuron::~SSRMNeuron() 
{
	delete _soma;
}

// Initialize the new instance
void SSRMNeuron::initialize()
{
	using namespace UOM;

	// Add the sole compartment
	add( _soma = new SSRMSoma );

	// Now override for SSRM defaults
	soma()->refractoryInterval(3*msec);
	soma()->spikeWidth(1*msec);
	soma()->RmRefractoryRatio(Number(0.01));

	// Set up excitatory and inhibitory conductances
	addSynapticConductances();
}

// Set up for synapses (subclass may override)
void SSRMNeuron::addSynapticConductances()
{
	soma()->ampa(new ExcitatorySynapticConductance);
	soma()->gaba(new InhibitorySynapticConductance);
}

// Return a solver to use if none was created earlier
ODESolver* SSRMNeuron::defaultSolver()
{
	const SimTime		defaultTimeStep = 0.25*UOM::msec;
	ODESolver*			mySolver = new AdaptiveSolver;

	mySolver->maxTimeStep(defaultTimeStep);
	mySolver->timeStep(defaultTimeStep);
	mySolver->errTol(0.001f);

	return mySolver;
}

// Add an excitatory synapse from here to a neuron
Synapse* SSRMNeuron::addExcSynapseFrom(SpikingNeuron* afferentNeuron, Number w, Number dist)
{
	return soma()->addExcSynapseFrom(afferentNeuron,w,dist);
}

// Add an inhibitory synapse from here to a neu=ron
Synapse* SSRMNeuron::addInhSynapseFrom(SpikingNeuron* afferentNeuron, Number w, Number dist)
{
	return soma()->addInhSynapseFrom(afferentNeuron,dist,w);
}

// Add an excitatory synapse to an efferent neuron and return the synapse
Synapse* SSRMNeuron::addExcSynapseTo(SSRMNeuron* efferentNeuron, Number w, Number dist)
{
	return efferentNeuron->addExcSynapseFrom(this,w,dist);
}

// Add an inhibitory synapse to an efferent neuron and return the synapse
Synapse* SSRMNeuron::addInhSynapseTo(SSRMNeuron* efferentNeuron, Number w, Number dist)
{
	return efferentNeuron->addInhSynapseFrom(this,w,dist);
}

// Handle a firing event by adjusting membrane voltage
void SSRMNeuron::postFiringActions(SimTime tspike)
{
	// Set the new membrane voltage
	soma()->setVmAfterSpike(tspike);
}



// ====================================================================
// LIAFNeuron class body
// ====================================================================



// Create a new instance with a default model and solver
LIAFNeuron::LIAFNeuron(bool doInit) : SSRMNeuron(false)
{
	// Do initialization during construction if requested
	if (doInit) {
		initialize();
	}
}

// Create a new instance using an existing model
LIAFNeuron::LIAFNeuron(Model* m, bool doInit) : SSRMNeuron(m,false)
{
	// Do initialization during construction if requested
	if (doInit) {
		initialize();
	}
}

// Destroy this instance
LIAFNeuron::~LIAFNeuron() {}

// Initialize the instance as a LIAF neuron
// Note that SSRMNeuron initialize is not done during
// construction of this object.
void LIAFNeuron::initialize()
{
	// Add the sole compartment
	add( _soma = new SSRMSoma );

	// Set up excitatory and inhibitory conductances
	soma()->ampa(new ExcitatorySynapticCurrent);
	soma()->gaba(new InhibitorySynapticCurrent);	

	// Override for LIAF defaults
	soma()->refractoryInterval(0);
	soma()->RmRefractoryRatio(1);
}



// ====================================================================
// Poisson Neuron Class
// ====================================================================



// Constructors and destructor
PoissonNeuron::PoissonNeuron(Model* m) : SpikingNeuron(m)
{
	using namespace UOM;

	// Set defaults (order matters here)
	peakFiringRate(100/sec);
	minISI(5*UOM::msec);
	firingRate(0);

	// At this stage, there are no spikes scheduled.
	_nextSpikeTime = InfinitePast;
}

PoissonNeuron::~PoissonNeuron() {}

// Create a clock solver to trigger end of time step
ODESolver* PoissonNeuron::defaultSolver()
{
	// Only a fixed interval solver is supported.
	ODESolver*		mySolver = new ClockSolver;

	// Set a default time step. This limits only the
	// times at which rates can be changed and not
	// the generated spike rate.
	mySolver->timeStep(5*UOM::msec);

	return mySolver;
}

// Set a new peak firing rate and compute a new next
// firing time if one has already been set.
void PoissonNeuron::peakFiringRate(Number rate)
{
	_peakFiringRate = rate;

	if (_nextSpikeTime>InfinitePast) {
		_nextSpikeTime = currentTime()+rexp(1/rate);
	}
}

// Set a firing rate. If necessary, increase peak firing
// to ensure that it is greater than the current firing rate.
void PoissonNeuron::firingRate(Number rate)
{
	// Save the rate provided
	_firingRate = rate;

	// Increase rate of candidate spikes if necessary
	if (rate>peakFiringRate() ) {
		peakFiringRate(rate);
	}
}

// Return the probability of selecting a spike among those 
// candidate spikes generated at peakFiringRate.
Number PoissonNeuron::spikeSelectionProbability()
{
	// Return the adjusted probability. Subclasses
	// may want to do something more elaborate
	// before making the adjustment, hence the
	// extra function call to get the rate.

	return isiAdjustedSelProb( firingRate() );
}

// Return an adjusted selection probability taking into
// account ISI restrictioned times when no firing can occur.
Number PoissonNeuron::isiAdjustedSelProb(Number fr)
{
	// See what fraction of the time is not taken up by
	// ISI restrictions.
	Number a=1-fr*minISI();
	
	// Return a rate with the adjustment. If no adjustment
	// is possible, fire every time just to keep going.
	return a<=0 ? 1: fr/a/peakFiringRate();
}

// At end of a time step, figure out which spikes come next
// and dispatch them in advance.
void PoissonNeuron::timeStepEnded()
{
	SimTime			now = currentTime();
	SimTime			next = now+timeStep();
	SimTime			meanISI = 1.0/peakFiringRate();

	// See if new peak firing spikes are needed. Since
	// candidate spikes are Poisson distributed, new
	// spike times are determined irrespective of history.
	if (_nextSpikeTime < now) {
		_nextSpikeTime = now + rexp(meanISI);
	}

	// Generate any spikes that fit into the next time step
	while (_nextSpikeTime<next) {

		// See if the next candidate time is sufficiently
		// past the previous accepted spike time to meet
		// the minISI constraint. If so, see if this spike 
		// is selected based on selection probability.

		if (_nextSpikeTime-spikeTime()>minISI() &&
			runif() < spikeSelectionProbability() ) {

			// Signal that a spike occurred
			signalSpikeEvent(_nextSpikeTime);
		}

		// Generate next peak firing spike time for candidate spikes
		_nextSpikeTime += rexp(meanISI);
	}

	// Let the superclass notify any probes that the step ended.
	SpikingNeuron::timeStepEnded();
}



// ====================================================================
// Static values for Synaptic Conductances
// ====================================================================



DualExpSynapticCondClassCache ExcitatorySynapticConductance	:: _CC;
DualExpSynapticCondClassCache ExcitatorySynapticCurrent		:: _CC;
DualExpSynapticCondClassCache InhibitorySynapticConductance	:: _CC;
DualExpSynapticCondClassCache InhibitorySynapticCurrent		:: _CC;
