// Leaky Integrate and Fire and Spike Response Neuron Models
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// File: bnsf_liaf.h
//
// Description:
//
// This header file contains the classes used to implement
// a leaky integrate and fire neuron model, a simple spike
// response neuron model and other models of theoretical interest. 
//
// While the classic leaky integrate and fire (LIAF) neuron is included
// here, with refactory behavior and realistic synaptic conductances,
// this implementation is more like a simplified version of Gerstner's
// spike response model (SSRM). The leaky integrate and fire
// implementation is a specialization of the SSRM model. 
//
// Synaptic conductances defined here are fairly elementary. These
// are patterned after AMPA and GABA (fast) synapses. For SSRM
// neurons these are conductances. For LIAF neurons, these are the 
// pure current sources more typically used in such models.
//
// Poisson neurons fire at random with a specified, but changable, firing
// frequency. The implementation is through the use of spike thinning.
// Unlike most neurons, firing events are scheduled to occur at a later
// time than when actually generated. This allows neurons receiving these 
// spikes to have them already queued up well in advanced of the time at 
// which the spike occurs. A consequence is that changes in firing frequency
// take effect one time step after they are made.
//
// Reference:
//
// These models are widely used and there are many possible references.
// A few selected ones are below.
//
// Dayan P and Abbott LF (2001). Theoretical Neuroscience: computational
// and mathematical modeling of neural systems.  Cambridge MA: MIT Press.
//
// Gerstner W. (1997). Spiking Neurons. In: Pulsed Neural
// Networks, ed. by Bishop, CM and Maass, W. Cambridge MA: MIT Press.
//
// Gerstner W and Kistler W (2002). Spiking Neuron Models.
// New York: Cambridge University Press.
//
// Dextexhe A, Mainen ZF, Sejnowski TJ (1998). Kinetic models of synaptic
// transmission, in Methods of Neuronal Modeling 2nd edition, ed. Kock C
// and Segev I. Cambridge MA: MIT Press.
//
// Protopapas AD, Vanier M, Bower J (1998). Simulating Large Networks
// of Neurons, in Methods of Neuronal Modeling 2nd edition, ed. Kock C
// and Segev I. Cambridge MA: MIT Press.
 

// Only include this header once
#ifndef __BNSF_LIAF_H_
#define __BNSF_LIAF_H_


// Required BNSF headers
#include "bnsf_base.h"
#include "bnsf_math.h"
#include "bnsf_sim.h"
#include "bnsf_nmod.h"

// Names spaces used in this definition
using namespace std;

// Add these definitions to the BNSF namespace
namespace BNSF {

	// Identify classes used in this application
	// to allow forward references and for documentation.

	class Compartment;					// defined elsewhere
		class SphericalCompartment;		// defined elsewhere
			class SSRMSoma;

	class Neuron;						// defined elsewhere
		class SpikingNeuron;			// defined elsewhere	
			class SSRMNeuron;
				class LIAFNeuron;
			class PoissonNeuron;	

	class ExcitatorySynapticConductance;
	class ExcitatorySynapticCurrent;

	class InhibitorySynapticConductance;
	class InhibitorySynapticCurrent;

	// ----------------------------------------------------------------
	// CLASS:	SSRMSoma
	// EXTENDS:	SphericalCompartment
	// DESC:	Define the soma compartment for both a
	//			simple spike response model and a leaky
	//			integrate and fire model.
	// RESP:
	//		1.	Provide parameters for spiking.
	//		2.	Allocate excitatory and inhibitory synapses.
	//		3.	Compute leakage currents immediately following
	//			refraction to simulate a fixed spike width
	//			followed by a refractory period when Rm is
	//			decreased by some fixed fraction.
	//
	// NOTES:	This is only one possible spike response model,
	//			but it can be tuned to provide a reasonable match
	//			to more realistic spiking neurons. LIAF neurons use
	//			a degrenerate case of the soma.
	// ----------------------------------------------------------------

	class SSRMSoma : public SphericalCompartment {

	public:
		
		// Constructors and destructor
		SSRMSoma();
		virtual ~SSRMSoma();

		// Accessors for excitatory and inhibitory synaptic conductances
		inline  SynapticConductance*	ampa() { return _ampa; }
		virtual void			ampa(SynapticConductance* synCond);

		inline  SynapticConductance* gaba() { return _gaba; }
		virtual void			gaba(SynapticConductance* synCond);

		// Accessor for the peak voltage following a spike. Default = +40 mV
		inline  Number			Vspike() { return _Vspike; }
		virtual void			Vspike(Number v) { _Vspike = v; }
		
		// Accessor for the potential after the spike. Default = -70 mV
		inline  Number			Vafter() { return _Vafter; }
		virtual void			Vafter(Number v) { _Vafter = v; }

		// Accessor for time interval following a spike during which Rm is
		// subject to adjustment to reduce voltage sensitivity.
		inline  SimTime			refractoryInterval() { return _refractoryInterval;}
		virtual void			refractoryInterval(Number t) { _refractoryInterval = t; }

		// Accessor for the time interval between the spike occurs and when
		// the refractory adjustment is applied.
		inline  SimTime			spikeWidth() { return _spikeWidth; }
		virtual void			spikeWidth(SimTime t) { _spikeWidth = t; }

		// Get/set a ratio used to set the membrane resistance during the
		// refractory period (if any). During the refractory period, Rm is
		// effectively Rm_after_spike = refractoryRatio*Rm_normal.
		inline  Number			RmRefractoryRatio() { return _RmRefractoryRatio; }
		virtual void			RmRefractoryRatio(Number r) { _RmRefractoryRatio = r; }

		// Add excitatory synpatic connection (AMPA) and return the synapse added
		virtual Synapse*		addExcSynapseFrom(
			SpikingNeuron* n,						// afferent neuron
			Number weight=1.0,						// synapse weight
			Number dist=100*UOM::micron);			// axon distance (used for delay)

		// Add inhibitory synpatic connection (GABA) and return the synapse added
		virtual Synapse*		addInhSynapseFrom(
			SpikingNeuron* n,						// afferent neuron
			Number weight=1.0,						// synapse weight
			Number dist=100*UOM::micron);			// axon distance (used for delay)

		// Framework interfaces -----------------------------------------------

		// Compute the leakage current. This is adjusted to simulate
		// a spike event and reduce responsiveness afterwards.
		virtual double			Ileak();

		// Adjust the current voltage to reflect a spiking event.
		virtual void			setVmAfterSpike(SimTime tspike);

		// State vector label functions
		virtual const char* componentName() {return "soma"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "v" }; 
			return sl; }

	protected:
		SynapticConductance*	_ampa;					// Excitatory AMPA conductances
		SynapticConductance*	_gaba;					// Inhibitory GABA conductances
		SimTime					_spikeWidth;			// Time during which spike occurs
		SimTime					_refractoryInterval;	// Period of limited response after spike
		Number					_RmRefractoryRatio;		// Ratio of Rm change (SSRM)
		Number					_Vspike;				// Spike peak potential (SSRM)
		Number					_Vafter;				// Potential following spike
	};

	// ----------------------------------------------------------------
	// CLASS:	SSRMNeuron
	// EXTENDS:	VoltageTriggeredSpikingNeuron
	// DESC:	Defines a spiking neuron containing a
	//			single compartment for the soma. This
	//			implements a restricted for of a spike
	//			response model (see SSRMSoma).
	// RESP:
	//		1.	Provide parameters for spiking.
	//		2.	Initialize the soma compartment
	//		3.	Allocate synapses (with soma compartment)
	//
	// NOTES:	Excitatory and inhibitory synaptic conductance
	//			are patterned after AMPA and fast GABA.
	// ----------------------------------------------------------------

	class SSRMNeuron : public VoltageTriggeredSpikingNeuron {

	public:

		// Constructors and destructor
		SSRMNeuron(bool doInit=true);
		SSRMNeuron(Model* m, bool doInit=true );

		virtual ~SSRMNeuron();	
		
		// Accessors
		virtual SSRMSoma*				soma() { return _soma; }
		virtual Compartment*			somaComp(int n=1) { return soma(); }

		virtual SynapticConductance*	ampa() { return soma()->ampa(); }
		virtual SynapticConductance*	gaba() { return soma()->gaba(); }

		// Add an excitatory synapse from an afferent neuron and return the synapse
		virtual Synapse*		addExcSynapseFrom(
			SpikingNeuron*			afferentNeuron, 
			Number					weight=1.0,
			Number					dist=100*UOM::micron); 

		// Add an inhibitory synapse from an afferent neuron and return the synapse
		virtual Synapse*		addInhSynapseFrom(
			SpikingNeuron*			afferentNeuron, 
			Number					weight=1.0,
			Number					dist=100*UOM::micron); 

		// Add an excitatory synapse to an efferent neuron and return the synapse
		virtual Synapse*		addExcSynapseTo(
			SSRMNeuron*				efferentNeuron, 
			Number					weight=1.0,
			Number					dist=100*UOM::micron); 

		// Add an inhibitory synapse to an efferent neuron and return the synapse
		virtual Synapse*		addInhSynapseTo(
			SSRMNeuron*				efferentNeuron, 
			Number					weight=1.0,
			Number					dist=100*UOM::micron); 

	protected:
		SSRMSoma*				_soma;
		virtual void			initialize();		// initialize with SSRM defaults
		virtual void			addSynapticConductances();	// set up for synapses

		// Spiking Neuron interface functions -------------------------

		// Default solver for this neuron
		virtual ODESolver*		defaultSolver();

		// Voltage to use for determining whether firing occurred
		virtual Number			firingVoltage() { return soma()->Vm(); }

		// Time derivative of firing voltage
		virtual Number			firingDeriv() { return soma()->VmDot(); }

		// Voltage threshold which must be exceeded to be counted as a spike
		virtual Number			firingThreshold() { return -50*UOM::mV; }

		// Minimum interspike interval (default)
		virtual SimTime			minISI() { return 3*UOM::msec; }

		// Handle reset after firing
		virtual void			postFiringActions(SimTime tspike);
	};

	// ----------------------------------------------------------------
	// CLASS:	LIAFNeuron
	// EXTENDS:	SpikingNeuron
	// DESC:	Defines a spiking neuron containing a
	//			single compartment for the soma. This
	//			implements a leaky integrate-and-fire
	//			neuron model.
	// RESP:
	//		1.	Provide parameters for spiking.
	//		2.	Initialize the soma compartment
	//		3.	Allocate synapses (with soma compartment)
	//
	// NOTES:	Excitatory and inhibitory synaptic conductances
	//			are patterned after AMPA and fast GABA.
	//			Synapses are modelled as pure current sources
	//			rather than conductances with a reversal potential.
	// ----------------------------------------------------------------

	class LIAFNeuron : public SSRMNeuron {

	public:

		// Constructors and destructor
		LIAFNeuron(bool doInit=true);
		LIAFNeuron(Model* m, bool doInit=true );
		virtual ~LIAFNeuron();

	protected:
		virtual void			initialize();		// initialize with LIAF defaults
	};

	// ----------------------------------------------------------------
	// CLASS:	PoissonNeuron
	// EXTENDS:	SpikingNeuron
	// DESC:	Defines a spiking neuron in which spikes
	//			occur at random such that the number of spikes
	//			in an interval follow Poisson statistics.
	// RESP:
	//		1.	Allocate ClockSolver as a source of time steps.
	//		2.	Generate candidate spike times for thinning
	//		3.	Generate random spikes for the upcoming time step.
	//		4.	Enforce minimum ISI value.
	//
	// NOTES:	Spike thinning is used to efficiently generate
	//			random spike trains. True Poisson statistics
	//			would require minimum ISI = 0, but this allows
	//			spikes separated by unrealisticly short intervals.
	//			A refractory period is introduced by adjusting spike
	//			selection probabilities to give the desired rate.
	// ----------------------------------------------------------------

	class PoissonNeuron: public SpikingNeuron {

	public:

		// Constructors and destructor
		PoissonNeuron(Model* m=NULL);
		virtual ~PoissonNeuron();

		// Accessor for the fixed time step. Spikes are generated at
		// the end of one time step to occur in the next so that
		// receiving neurons will not have late arriving events.
		inline  SimTime			timeStep() { return solver()->timeStep(); }
		virtual void			timeStep(SimTime step) { solver()->timeStep(step); }

		// Accessor for peak firing rate. This is the rate from which
		// spikes are thinned by the current firing rate. The solver
		// time step is adjusted to match the peak firing rate.
		inline  Number			peakFiringRate() { return _peakFiringRate; }
		virtual void			peakFiringRate(Number rate);

		// Accessor for current firing rate. Peak rate is increased if needed. 
		// Subclass may set firingRate if rate is determined dynamically. 
		virtual Number			firingRate() { return _firingRate; }
		virtual void			firingRate(Number rate);

		// Minimum interspike interval - default value here is 1 msec.
		// For true poisson statistics, minISI can be set to 0.
		virtual SimTime			minISI() { return _minISI; }
		virtual void			minISI(SimTime isi) { _minISI = isi; }

		// Framework interfaces ---------------------------------------
		virtual void			timeStepEnded();

	protected:
		SimTime					_nextSpikeTime;
		SimTime					_minISI;			// minimum ISI constraint
		Number					_peakFiringRate;	// rate for candidate spikes
		Number					_firingRate;		// effective firing rate

		// Apply a formula for adjusting spike selection probability
		// to compensate for imposed minISI restriction. fr is the
		// desired firing rate after allowing for minISI.
		virtual  Number			isiAdjustedSelProb(Number fr);

		// Framework interfaces ---------------------------------------
		virtual ODESolver*		defaultSolver();

		// Voltage to use for determining whether firing occurred (unused)
		virtual Number			firingVoltage() { return 0; }

		// Subclass responsibilities ----------------------------------

		// Return the probability of selecting a spike among those generated
		// at peakFiringRate and at the time of the next time. This is
		// adjusted slightly to compensate for ISI rules.
		virtual Number			spikeSelectionProbability();
	};

	// ----------------------------------------------------------------
	// CLASS:	ExcitatorySynapticConductance
	// EXTENDS:	DualExpSynapticCond
	// DESC:	Defines a synaptic conductance with timing
	//			typical of excitatory AMPA synapses.
	// RESP:
	//		1.	Provide time constants for synaptic response.
	//		2.	Provide other associated parameter values.
	//
	// ----------------------------------------------------------------

	class ExcitatorySynapticConductance : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		ExcitatorySynapticConductance(Number gMaxValue = 1*UOM::nanoS) 
			: DualExpSynapticCond(gMaxValue) {}
		
		virtual ~ExcitatorySynapticConductance () {}

		// Parameter accessors
		virtual SimTime		tau1() { return 2*UOM::msec; }
		virtual SimTime		tau2() { return 2*UOM::msec; }
		virtual Number		Vrev() { return 0*UOM::mV; }

		// State vector label functions
		virtual const char* componentName() {return "AMPA"; }

	protected:
		static  DualExpSynapticCondClassCache _CC;	// class cache
		virtual DualExpSynapticCondClassCache*		// accessor
							pDualExpClassCache() { return &_CC; }
	};

	// ----------------------------------------------------------------
	// CLASS:	InhibitorySynapticConductance
	// EXTENDS:	DualExpSynapticCond
	// DESC:	Defines a synaptic conductance with timing
	//			typical of inhibitory GABA synapses.
	// RESP:
	//		1.	Provide time constants for synaptic response.
	//		2.	Provide other associated parameter values.
	// ----------------------------------------------------------------

	class InhibitorySynapticConductance : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		InhibitorySynapticConductance(Number gMaxValue = 1*UOM::nanoS) 
			: DualExpSynapticCond(gMaxValue) {}

		virtual ~InhibitorySynapticConductance() {}

		// Parameter accessors
		virtual SimTime		tau1() { return 1*UOM::msec; }
		virtual SimTime		tau2() { return 7*UOM::msec; }
		virtual Number		Vrev() { return -80*UOM::mV; }

		// State vector label functions
		virtual const char* componentName() {return "GABA"; }

	protected:
		static  DualExpSynapticCondClassCache _CC;	// class cache
		virtual DualExpSynapticCondClassCache*		// accessor
							pDualExpClassCache() { return &_CC; }
	};

	// ----------------------------------------------------------------
	// CLASS:	ExcitatorySynapticCurrent
	// EXTENDS:	DualExpSynapticCond
	// DESC:	Defines a synaptic current with timing
	//			typical of excitatory AMPA synapses.
	// RESP:
	//		1.	Provide time constants for synaptic response.
	//		2.	Provide other associated parameter values.
	//		3.	Compute current as a function of time.
	// ----------------------------------------------------------------

	class ExcitatorySynapticCurrent : public DualExpSynapticCurrent {

	public:

		// Constructors and destructor
		ExcitatorySynapticCurrent(Number Imax=-65*UOM::picoA) 
			: DualExpSynapticCurrent(Imax) {}

		virtual ~ExcitatorySynapticCurrent () {}

		virtual SimTime		tau1() { return 2*UOM::msec; }
		virtual SimTime		tau2() { return 2*UOM::msec; }

		// State vector label functions
		virtual const char* componentName() {return "AMPA"; }

	protected:
		static  DualExpSynapticCondClassCache _CC;	// class cache
		virtual DualExpSynapticCondClassCache*		// accessor
							pDualExpClassCache() { return &_CC; }
	};

	// ----------------------------------------------------------------
	// CLASS:	InhibitorySynapticCurrent
	// EXTENDS:	DualExpSynapticCond
	// DESC:	Defines a synaptic current with timing
	//			typical of inhibitory GABA synapses.
	// RESP:
	//		1.	Provide time constants for synaptic response.
	//		2.	Provide other associated parameter values.
	//		3.	Compute current as a function of time.
	// ----------------------------------------------------------------

	class InhibitorySynapticCurrent : public DualExpSynapticCurrent {

	public:

		// Constructors and destructor
		InhibitorySynapticCurrent(Number Imax=15*UOM::picoA) 
			: DualExpSynapticCurrent(Imax) {}

		virtual ~InhibitorySynapticCurrent() {}

		virtual SimTime		tau1() { return 1*UOM::msec; }
		virtual SimTime		tau2() { return 7*UOM::msec; }

		// State vector label functions
		virtual const char* componentName() {return "GABA"; }

	protected:
		static  DualExpSynapticCondClassCache _CC;	// class cache
		virtual DualExpSynapticCondClassCache*		// accessor
							pDualExpClassCache() { return &_CC; }
	};
};

#endif // #ifndef
