// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_nmod.h
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		6 March 2007
//
// Description:
//
// This header file contains the basic classes used as a
// framework for building biological network simulations.
// As is normal for frameworks, many of these classes are
// abstract and subclassed as needed for the simulation.
//
// This header defines classes used in define neural models.


// Only include this header once per compilation unit
#ifndef __BSNF_NMOD_H_
#define __BSNF_NMOD_H_


// ====================================================================
// MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================
#ifdef WIN32

// Disable warning C4786: symbol greater than 255 character,
#pragma warning( disable: 4786)

#endif
// ====================================================================
// END OF MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================


// ====================================================================
// Header files included here by default
// ====================================================================

// Standard C++ Template Library headers
#include <functional>
#include <queue>
#include <map>
#include <list>
#include <set>

// Incorporate all the names from the std library by reference
using namespace std;

// Required BNSF headers
#include "bnsf_base.h"
#include "bnsf_math.h"
#include "bnsf_sim.h"


// ====================================================================
// Primary namespace for the framework
// ====================================================================


namespace BNSF {


	// ====================================================================
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ====================================================================


	// Extensions to the simulator hierarchy
	class Solver;							// defined elsewhere
		class ODESolver;					// defined elsewhere
			class NeuronSolver;
	
	class Probe;							// defined elsewhere
		class ExternalRecorder;				// defined elsewhere
			class ExternalStateRecorder;	// defined elsewhere
				class ExternalCalciumRecorder;
				class ExternalIonCurrentRecorder;
				class ExternalCompartmentRecorder;
					class ExternalVoltageRecorder;
					class ExternalCurrentRecorder;
			class ExternalSpikeRecorder;
			
	// Public classes defining table entries
	class AlphaBetaEntry;
	class GHKTableEntry;
	class MorphologyEntry;

	// Neural model classes
	class Neuron;
		class SpikingNeuron;
			class VoltageTriggeredSpikingNeuron;
				class MorphologicalNeuron;

	class Compartment;
		class SphericalCompartment;

	class AxonProcess;
	class ElectricalCoupling;
		class ElectricalJunction;

	class CalciumPool;
		class SimpleCalciumPool;

	// Classes for defining ion channels
	class IonChannel;
		class MultiGateIonChannel;
			class MnHIonChannel;
				class M1HIonChannel;
				class M2HIonChannel;
				class M3HIonChannel;
				class M4HIonChannel;
				class MhHCaIonChannel;
					class M1HCaIonChannel;
					class M2HCaIonChannel;
					class M3HCaIonChannel;
			class MnHSIonChannel;
				class M1HSIonChannel;
				class M2HSIonChannel;
				class M3HSIonChannel;
				class M4HSIonChannel;
		class HHIonChannel;
			class HHIonGate;
			class VoltageDepTabChannel;
				class VoltageDepTabGate;
				class EnergyBarrierTabChannel;
					class Order1EnergyBarrierTabChannel;
					class Order2EnergyBarrierTabChannel;
					class Order3EnergyBarrierTabChannel;
					class Order4EnergyBarrierTabChannel;
					class EnergyBarrierTabGate;
		class SynapticResponse;
			class SynapticGroupResponse;
			class SynapticConductance;
				class SingleExpSynapticCond;
					class DualExpSynapticCondWithVarTau;
					class DualExpSynapticCond;
						class DualExpSynapticCurrent;
					class TripleExpSynapticCond;

	// Classes related to synaptic connections
	class Synapse;
	class SynapseState;

	class PlasticityRule;
	class PlasticityState;

	// Network events
	class Event;							// defined elsewhere
		class ActionPotentialEvent;			// eventClassId=0x0001

	class ActionPotentialEventQueue;


	// ====================================================================
	// Common typedefs and globals
	// ====================================================================

	// ----------------------------------------------------------------
	// Vector and iterator typedefs for (most) classes
	// Note that by convention, pointers are stored in
	// vectors so this info is not in the typedef name.
	// ----------------------------------------------------------------

	typedef vector<Neuron*>							NeuronVector;
	typedef NeuronVector::iterator					NeuronVectorIt;

	typedef vector<Compartment*>					CompartmentVector;
	typedef CompartmentVector::iterator				CompartmentVectorIt;

	typedef vector<Synapse*>						SynapseVector;
	typedef SynapseVector::iterator					SynapseVectorIt;

	typedef list<Synapse*>							SynapseList;
	typedef SynapseList::iterator					SynapseListIt;

	typedef vector<ElectricalCoupling*>				ElectricalCouplingVector;
	typedef ElectricalCouplingVector::iterator		ElectricalCouplingVectorIt;

	typedef vector<CalciumPool*>					CalciumPoolVector;
	typedef CalciumPoolVector::iterator				CalciumPoolVectorIt;

	typedef vector<SynapticResponse*>				SynapticResponseVector;
	typedef SynapticResponseVector::iterator		SynapticResponseVectorIt;

	typedef vector<SynapticGroupResponse*>			SynapticGroupResponseVector;
	typedef SynapticGroupResponseVector::iterator	SynapticGroupResponseVectorIt;

	typedef vector<SynapticConductance*>			SynapticConductanceVector;
	typedef SynapticConductanceVector::iterator		SynapticConductanceVectorIt;

	typedef vector<PlasticityRule*>					PlasticityRuleVector;
	typedef PlasticityRuleVector::iterator			PlasticityRuleVectorIt;

	typedef vector<IonChannel*>						IonChannelVector;
	typedef IonChannelVector::iterator				IonChannelVectorIt;

	typedef vector<ActionPotentialEvent*>			ActionPotentialEventVector;
	typedef ActionPotentialEventVector::iterator	ActionPotentialEventVectorIt;

	typedef vector<ActionPotentialEvent>			APEventObjectVector;
	typedef APEventObjectVector::iterator			APEventObjectVectorIt;

	// ====================================================================
	// Global functions used to compute index into lookup
	// tables for alpha beta values etc. These are global
	// because all tabels use the same index system to
	// avoid multiple conversions from floating point
	// voltage to integer index, which can be an expensive
	// operation. The first and last two entries of the
	// table are not directly accessed via index values
	// but may be used by interpolation functions.
	// ====================================================================

	static const Number VMinForIndex		= -120*UOM::mV;
	static const Number VMaxForIndex		= 120*UOM::mV; 
	static const Number VStepForIndex		= 0.125*UOM::mV;
	static const int	VTableIndexMax		= int((VMaxForIndex-VMinForIndex)/VStepForIndex-2.5);
	static const int	VTableSize			= VTableIndexMax+3;

	// Determine the table index from a voltage
	inline int VTableIndex(Number v) 
	{
		// Calculate an integer index from v
		const int index = int((v-VMinForIndex)/VStepForIndex);

		// Ensure that the index is in range. Index values 0, 
		// VIndexMax+1, and VIndexMax+2 are reserved to simplify 
		// interpolation (see VTableInterp below).
		if (index<1)					return 1;
		else if (index>VTableIndexMax)	return VTableIndexMax;
		else							return index;
	}

	// Compute a remainder for use in interpolating from voltage
	inline Number VTableRem(Number v, int vmIdx)
	{
		const Number vrem = v-(VMinForIndex+vmIdx*VStepForIndex);

		// Check for cases where v is out of range of the table
		if (vrem<0)						return 0;
		else if (vrem>VStepForIndex)	return VStepForIndex;
		else							return vrem;
	}

	// Estimate y=f(v) using an interpolation of table values from
	// from voltages in the range VMinForIndex through VMaxForIndex.
	// Let v=v1+s*h where h=v(i+1)-v(i) and y(i)=f(x(i)).
	// Typically 0<=s1<1 for best accuracy.
	
	// Estimate y=f(v) by piecewise linear interpolation.
	inline Number VTableInterpLI(Number s, Number y1, Number y2)
	{
		return y1+s*(y2-y1);
	}

	// Estimate y=f(v) with a piecewise cubic polynomial.
	// Note that the first derivative is not continuous at 
	// knot points. This is based on Neville's method but
	// has been restated for a fixed h value with a bit
	// of factoring thrown in.
	inline Number VTableInterpPC(Number s, 
		Number y0, Number y1, Number y2, Number y3)
	{
		return y1+s/2*(y2-y0+s*(y2-y1+y0-y1)+(1-s*s)*(y2-y1+(y0-y3)/3));
	}

	// Estimate y=f(v) by averaging two quadratic fits based
	// on y0,y1,y2 and y1,y2,y3. This is somewhat less accurate
	// than a cubic polynomial but error terms are less biased.
	// It also uses fewer floating point operations.
	inline Number VTableInterpQA(Number s, 
		Number y0, Number y1, Number y2, Number y3)
	{
		return y1+s*(y2-y1+(s-1)*(y3+y0-y2-y1)/4);
	}

	// Provide a common default interpolation method that 
	// can easily be changed as needed.
	inline Number VTableInterp(Number s, 
		Number y0, Number y1, Number y2, Number y3)
	{
		return VTableInterpQA(s,y0,y1,y2,y3);
	}



	// ====================================================================
	// Framework classes for neural simulation objects
	// ====================================================================



	// ----------------------------------------------------------------
	// CLASS:	NeuronSolver
	// EXTENDS:	ODESolver
	// DESC:	Numerically evaluate the ODE system using
	//			a method specialized for neurons. A Crank-Nicolson
	//			scheme of staggered evaluations is used to achieve
	//			roughly second order accuracy. Then Richardson
	//			extrapolation is used to increase accuracy and
	//			adapt the procedure to maintain an error tolerance.
	//			Error is estimated by comparing the results of
	//			different order exptrapolations.
	//			
	// RESP:
	//		1.	Build Jacobian for compartments
	//		2.	Apply semi-implicit solution to state transitions
	//		3.	Apply Richardson extrapolation to produce final state
	//		4.	Develop a posteriori error estimates
	//
	// NOTE:	This is related to the semi-implicit method
	//			described by Press et al. in Numerical Recipes.
	//			See Bader and Deuflhard's modification of the
	//			method by Bulirsch-Stoer for stiff systems.
	//
	//			If high-order extrapolation is used, changing
	//			Number to type double may be needed. Otherwise
	//			floating point round-off inhibits convergence.
	//
	//			See Liem CB, Lu T, and Shih TM, The Splitting
	//			Extrapolation Method for an overview of Richardson
	//			extrapolation theory. See Mascagni MV and Sherman AS
	//			chapter in Methods in Neuronal Modeling for
	//			a description of the Crank-Nicolson method.
	// ----------------------------------------------------------------

	class NeuronSolver : public ODESolver {

	public:

		// Constructors and destructor
		NeuronSolver();
		virtual ~NeuronSolver();

		// Paramter accessors -----------------------------------------

		// timeStep, errTol, and safetyMargin are inherited.

		// timeStep controls the macro step size, that is, time 
		// interval that is crossed by using variable numbers of 
		// micro steps in separate passes over the interval. 
		// A counter-intuitive property here is that decreasing 
		// timeStep may not necessarily decrease cumulative error.

		// maxPasses controls how many passes are used in the
		// extrapolation process. Additional passes may be made
		// to achieve error tolerance, but only maxPasses are used
		// in the final extrapolation.
		inline  int				maxPasses() { return _maxPasses; }
		virtual void			maxPasses(int n);

		// minSteps controls the minimum number of steps used.
		inline  int				minSteps() { return _minSteps; }
		virtual void			minSteps(int n);

		// maxSteps controls how many steps can be used in
		// attempting to cross the interval specified in timeStep.
		// If this is exceeded, the interval is halved.
		inline	int				maxSteps()	{ return _maxSteps; }
		virtual void			maxSteps(int n) { _maxSteps = n; }

		// Always assume only first order accuracy when doing extrapolations
		// regardless of component specifications. Note that compartment
		// voltages are assumed second order accurate unless this is set.
		// Default is to assume at least second order accuracy unless otherwise
		// specified by the component by returning isSecondOrderAccurate()==false.
		inline  bool			assumeFirstOrder() { return _assumeFirstOrder; }
		virtual void			assumeFirstOrder(bool x) { _assumeFirstOrder=x; }

		// AlwaysRebuildJacobian controls whether the Jacobian
		// matrix for compartments is retained between steps.
		// Retaining the matrix prevents redundant computation at
		// the expense of additional storage.
		inline  bool			alwaysRebuildJacobian() { return _alwaysRebuildJacobian; }
		virtual void			alwaysRebuildJacobian(bool f) { _alwaysRebuildJacobian=f; }

		// Other public interfaces ------------------------------------

		// Return the estimated error in the previous time step.
		// This is a weighted worst case error in any component.
		inline  double			errorEst() {return _errorEst; }

		// Tell all components we are starting now
		virtual void			start();				// initialize

	protected:

		CompartmentVector		_compartments;			// all compartments in model
		ModelComponentVector	_timeAlignedComponents;	// components evaluated in phase with Vm
		ModelComponentVector	_timeOffsetComponents;	// components evaluated out of phase

		SparseMatrix			_jacobian;				// compartmental Jacobian

		int						_compSVOffset;			// state vector offset to first compartment
		int						_compSVSize;			// size of compartment state vector
		int						_maxPasses;				// max passes over interval
		int						_curPasses;				// current number of passes used
		double					_errorEst;				// estimate of error on last step

		// Note: the following arrays should agree with _PerPassArraySize
		// which cannot be used here for allocations.

		int						_nsteps[8];				// number of steps per pass
		int						_maxSteps;				// maximum steps per pass
		int						_minSteps;				// minimum steps per pass
		double					_hvect[8];				// micro time step sizes
		double					_stateExtrap[2][8];		// coeff for state extr
		double					_errorExtrap[2][8];		// coeff for error extr
		bool					_assumeFirstOrder;		// use 1st order extrap. formula
		bool					_alwaysRebuildJacobian; // flag for optional rebuild

		// Size of arrays allocated for per-pass data
		const static int		_PerPassArraySize;	// i.e. 8

		// Look through model components and sort by recipient groups
		virtual void			prepareRecipients();

		// Do one time step of the duration indicated
		virtual void			processStep(SimTime maxDuration);

		// Compute implicit update by solving the linear system.
		// LU decomp is used to solve the system. 
		// compDeriv is both input and output vector.
		virtual void			solveImplicitSystem(SparseMatrix& IhJ, double* compDeriv);

		// Build vector of time steps based on macro step time
		// and number of steps per pass
		virtual void			buildHvect(SimTime H);

		// Build interpolation coefficients given a current 
		// macro step size and the array of steps per pass.
		virtual void			buildExtrapolationCoeff();

		// Assign Jacobian indexes to minimize LU decomp fill
		// in a manner similar to Hines ordering.
		virtual void			assignJacobianIndexes();

		// Build the Jacobian matrix by interrogating compartments
		virtual void			buildJacobian();

		// Update the Jacobian matrix by interrogating compartments
		virtual void			updateJacobian();

		// Compute derivatives based on current state data
		// This applies only to compartments
		virtual void			computeDerivatives();

		// Update non-compartment states via local methods
		virtual void			updateOffsetStates(SimTime h, CNStepType stepType);
		virtual void			updateAlignedStates(SimTime h, CNStepType stepType);

		// Notify compartments that their state vector has been changed
		virtual void			notifyOnStateVectorChanged();
	};

	
	// ----------------------------------------------------------------
	// CLASS:	ExternalSpikeRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records spiking events in an external file.
	// RESP:
	//		1.	Write out id and spike time.
	//
	// NOTE:	Minimum reporting interval does not not apply
	//			since all spikes should be recorded.
	// ----------------------------------------------------------------

	class ExternalSpikeRecorder : public ExternalRecorder {

	public:

		// Constructors and destructor
		ExternalSpikeRecorder(char* fn);
		virtual ~ExternalSpikeRecorder();

		// Signal that an action potential occurred.
		// neuron is assumed to be a subclass of SpikingNeuron.
		virtual void signalEvent(unsigned int classId, ModelComponent* neuron);
	
	protected:

		// Provide a stub for reportOn since the interface is not used
		virtual void reportOn(ModelComponentVector& comps, SimTime t, int id) {}
	};
	
	// ----------------------------------------------------------------
	// CLASS:	ExternalDendriteSpikeRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records spiking events in an external file.
	// RESP:
	//		1.	Test for dendrite voltage exceeding a threshold
	//		2.	Determine if this follow a whole cell spike
	//		3.	Write a record indicating the spike event
	//
	// NOTES:	The neuron reported on is found in the collection of
	//			compartments provided to reportOn. Only the first
	//			found MorphologicalNeuron is reported on.
	//
	//			Voltages are tested only at the end of time steps.
	//			The record written contains neuron id, time, and dendrite
	//			number, and time since the previous cell spike. 
	//			Each occurence of transition from below the voltage threshold
	//			to above it is reported. Dendrite spike times are reported
	//			as of the end of the time step. Times are reported in units 
	//			of milliseconds. The soma is used as the origin point of
	//			backpropagating spikes because those in axons might not
	//			completely propagate back into the soma and dendrites.
	// ----------------------------------------------------------------

	class ExternalDendriteSpikeRecorder : public ExternalRecorder {

	public:

		// Constructors and destructor
		ExternalDendriteSpikeRecorder(
			char*				fn,						// file name to write to
			Number				spikeVm=-30*UOM::mV);	// voltage threshold for spike
		virtual ~ExternalDendriteSpikeRecorder();

		// Accessors
		inline	Number			spikeVm() { return _spikeVm; }
		virtual void			spikeVm(Number v) { _spikeVm=v; }

		// Note when removed from a component
		virtual void			removedFrom(ModelComponent* mc);

	protected:
		Number					_spikeVm;		// voltage threhold
		MorphologicalNeuron*	_neuron;		// neuron reported on				
		int						_nDend;			// number of dendrite comps

		// Work areas to keep track of previous voltages
		Number*					_prevVm;		// previous Vm by dendrite
		Number					_prevSomaVm;	// previous Vm at soma
		SimTime					_somaSpikeTime;	// time of last soma spike

		// Set state to an initial condition
		virtual void			initialize();
		
		// Report on any spike events found
		virtual void			reportOn(
			ModelComponentVector&	comps,		// components to report on
			SimTime					t,			// time of probe
			int						id);		// id of component probed

		// Identify a root node where backpropagating action potentials
		// are found before those in dendrites. Spikes that initiate in 
		// the axon can come after those found in the soma based on 
		// the voltage threshold used here. Subclasses can provide
		// alternative designations as needed.
		virtual Compartment*	somaComp();

		// Write a column header
		virtual void			writeColumnHeader();

		// Check for dendritic spikes and write them out
		virtual void			writeSpikes(SimTime t);
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalSynapseRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records synpatic weights in a comma separate format.
	// RESP:
	//		1.	Write the weight of each synapse to a file
	//
	// NOTES:	Since synapses can be transient, each weight is
	//			written as a separate line. This can result in many
	//			lines per reporting epoch. Since there can be many
	//			synapses, this is not a report that should be written
	//			with great frequency.
	// ----------------------------------------------------------------

	class ExternalSynapseRecorder : public ExternalRecorder {

	public:

		// Constructors and destructor
		ExternalSynapseRecorder(char* fn, TokenId synType=NullTokenId);
		ExternalSynapseRecorder(char* fn, char* synName);
		virtual ~ExternalSynapseRecorder();

		// Get/set the synapse type to report on via token or string
		inline  TokenId				synapseType() { return _synapseType; }
		virtual void				synapseType(TokenId type) { _synapseType = type; }
		virtual void				synapseType(char* name) { _synapseType = token(name); }

		// Report on a collection of components
		virtual void				reportOn(ModelComponentVector& comps, SimTime t, int id);

	protected:
		TokenId						_synapseType;	// component id of synapse to report
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalCompartmentRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records a compartment related value in external file
	//			using comma separated value (CSV) format. This class
	//			is an abstract class.
	// RESP:
	//		1.	Write header line with compartment names.
	//		2.	Write common header for each value line.
	//		3.	Write compartment value (via subclass)
	// ----------------------------------------------------------------

	class ExternalCompartmentRecorder : public ExternalStateRecorder {

	public:

		// Constructors and destructor
		ExternalCompartmentRecorder(char* fn, char* mapFileName = NULL);
		virtual ~ExternalCompartmentRecorder();

	protected:
		virtual void				writeColumnHeader(ModelComponentVector& comp);
		virtual void				writeState(ModelComponentVector& comp, SimTime t, int id);

		// Subclass responsibility
		virtual void				writeValue(Compartment* pcomp) = 0;
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalVoltageRecorder
	// EXTENDS:	ExternalCompartmentRecorder
	// DESC:	Records compartment voltages in external file
	//			using comma separated value (CSV) format.
	// RESP:
	//		1.	Write state vector values for the voltage
	//			for each compartment of the model associated
	//			with this recorder.
	// ----------------------------------------------------------------

	class ExternalVoltageRecorder : public ExternalCompartmentRecorder {

	public:

		// Constructors and destructor
		ExternalVoltageRecorder(char* fn, char* mapFileName = NULL);
		virtual ~ExternalVoltageRecorder();

	protected:
		virtual void				writeValue(Compartment* pcomp);
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalCurrentRecorder
	// EXTENDS:	ExternalCompartmentRecorder
	// DESC:	Records compartment currents in external file
	//			using comma separated value (CSV) format.
	// RESP:
	//		1.	Write total current for each compartment of the 
	//			model associated with this recorder.
	//
	// NOTES:	Total current is computed during computeDerivatives
	//			and is reported using that value. Even though a
	//			compartment may under voltage clamp, the current
	//			is still computed even though dV/dt=0 if clamped.
	// ----------------------------------------------------------------

	class ExternalCurrentRecorder : public ExternalCompartmentRecorder {

	public:

		// Constructors and destructor
		ExternalCurrentRecorder(char* fn, char* mapFileName = NULL);
		virtual ~ExternalCurrentRecorder();

	protected:
		virtual void				writeValue(Compartment* pcomp);
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalIonCurrentRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records state variables in external file
	//			using comma separated value (CSV) format.
	// RESP:
	//		1.	Write state vector values for the current
	//			for each ion channel of the model associated
	//			with this recorder.
	// ----------------------------------------------------------------

	class ExternalIonCurrentRecorder : public ExternalStateRecorder {

	public:

		// Constructors and destructor
		ExternalIonCurrentRecorder(char* fn, char* mapFileName = NULL);
		virtual ~ExternalIonCurrentRecorder();

	protected:
		virtual void				writeColumnHeader(ModelComponentVector& comp);
		virtual void				writeState(ModelComponentVector& comp, SimTime t, int id);
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalCalciumRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records state variables in external file
	//			using comma separated value (CSV) format.
	// RESP:
	//		1.	Write state vector values for the calcium
	//			pools of the model associated with this recorder.
	// ----------------------------------------------------------------

	class ExternalCalciumRecorder : public ExternalStateRecorder {

	public:

		// Constructors and destructor
		ExternalCalciumRecorder(char* fn);
		virtual ~ExternalCalciumRecorder();

	protected:
		virtual void				writeColumnHeader(ModelComponentVector& comp);
		virtual void				writeState(ModelComponentVector& comp, SimTime t, int id);
	};
	
	// ----------------------------------------------------------------
	// CLASS:	Neuron
	// EXTENDS:	ModelComponent
	// DESC:	The primary active brain cell type simulated.
	// RESP:
	//		1.	Know cell compartments
	//		2.	Connect compartments electrotonically
	//		3.	Set compartment topological distances
	//		4.	Notify probes when state changes
	// ----------------------------------------------------------------

	class Neuron : public ModelComponent {

	public:

		// Constructors and destructors
		Neuron(Model* m = NULL);
		virtual ~Neuron();

		// Model accessors
		inline Model*			model() { return _model; }
		virtual void			model(Model* m);

		// ODE solver accessors
		virtual ODESolver*		solver();
		virtual void			solver(ODESolver* newSolver);

		// Accessor for numeric identifier
		virtual int				numericIdentifier() { return _numericIdentifier; }
		virtual void			numericIdentifier(int id) { _numericIdentifier = id; }

		// Add-remove compartments
		virtual void			addCompartment(Compartment* comp);
		virtual void			removeCompartment(Compartment* comp);
		inline  int				numComp() {return _compartments.size(); }

		// Add-remove probes used to trace state changes and events
		virtual void			addProbe(Probe* pr);
		virtual void			removeProbe(Probe* pr);

		// Polymorphic add for convenience
		virtual void			add(Probe* pr) { addProbe(pr); }
		virtual void			add(Compartment* c) { addCompartment(c); }

		// Provide total membrane area by summing of all compartments
		virtual Number			membraneArea();

		// By default the neuron does not keep its own state
		virtual int				numStateVar() { return 0; }
		virtual void			setInitialState() {}
		virtual void			computeDerivatives() {}

		// Tell about this object if asked
		virtual bool			isTopLevelComponent() { return true; }
		virtual bool			isNeuron() { return true; }
		virtual bool			notifyOnTimeStepEnded() { return true; }

		// Subclass responsibilities (some optional) ------------------

		// Identify special categories of neurons
		virtual bool			isSpikingNeuron() { return false; }
		virtual bool			isMorphologicalNeuron() { return false; }

		// Interface accessors for spiking
		virtual SimTime			spikeTime()		{ return InfinitePast; }
		virtual AxonProcess*	axonProcess()	{ return NULL; }

		// Default ODESolver to use when neuron owns the model.
		virtual ODESolver*		defaultSolver();

		// Notifications from ODE solver
		virtual void			simulationStarted();
		virtual void			simulationEnded();
		virtual void			timeStepEnded();

	protected:

		int						_numericIdentifier;		// an arbitrary identifier
		ProbeVector				_probes;				// vector of associated probes
		CompartmentVector		_compartments;			// vector of compartments
		ODESolver*				_allocatedSolver;		// solver allocated by this neuron
		bool					_usesOwnModel;			// allocate own model and solver

		// Delete any held compartment objects
		virtual void			deleteAllCompartments();

		// Utility functions for subclasses ---------------------------

		// Assign compartments their distance from a root compartment
		virtual void			setDistFromSoma(
			Compartment*		root = NULL);			// NULL means compartments[0]
	};

	// ----------------------------------------------------------------
	// CLASS:	SpikingNeuron
	// EXTENDS:	Neuron
	// DESC:	Abstract class representing a neuron
	//			capable of action potential spikes.
	// RESP:
	//		1.	Know associated axon process
	//		2.	Notify associated efferent synapses of spikes
	//		3.	Estimate of the instantaneous firing rate (subclass)
	//		4.	Estimate of the ongoing firing rate (subclass)
	//
	// NOTE:	Axon processes are insulated from the main neuron
	//			objects to allow future expansion for distributed
	//			processing. Otherwise this class just provides
	//			common protocols for subclasses that can be very
	//			different in their implementations.
	// ----------------------------------------------------------------

	class SpikingNeuron : public Neuron {

	public:

		// Constructors and destructor
		SpikingNeuron(Model* m = NULL);
		virtual ~SpikingNeuron();

		// Accessors
		virtual bool			isSpikingNeuron() { return true; }
		inline  SimTime			spikeTime() { return _spikeTime; }
		inline  AxonProcess*	axonProcess() { return _axonProcess; }

		inline  int				numericIdentifier() { return _numericIdentifier; }
		virtual void			numericIdentifier(int id);

		// Subclass responsibilities ----------------------------------

		// Return an estimate of instanteous firing rate.
		virtual Number			firingRate() = 0;

	protected:
		SimTime					_spikeTime;		// time of most recent spike
		AxonProcess*			_axonProcess;	// connection with synapses

		// Signal that a spike occurred, notifying efferent neurons and probes
		virtual void			signalSpikeEvent(SimTime t);
	};

	// ----------------------------------------------------------------
	// CLASS:	VoltageTriggeredSpikingNeuron
	// EXTENDS:	SpikingNeuron
	// DESC:	Abstract class representing a neuron
	//			capable of action potential spikes.
	// RESP:
	//		1.	Detect spiking events via voltage threshold
	//		2.	Make an estimate of the instantaneous firing rate
	//
	// NOTE:	Defining when a given voltage constitutes
	//			a spike is not an absolute. The objective
	//			is to create an event that correlates well
	//			with efferent presynaptic terminals.
	//			Defaults chosen here heuristically detect a
	//			rising voltage above a firing threshold.
	//			Depending on characteristics of the
	//			axon and the presynaptic terminal, different
	//			rules might apply.
	//
	//			Ongoing firing rate is not really different
	//			than instantaneous firing rate except that the
	//			interval over which activity is averaged is
	//			longer. Default algorithms are provided for
	//			estimating these values but they can be
	//			overridden in subclasses as necessary.
	// ----------------------------------------------------------------

	class VoltageTriggeredSpikingNeuron : public SpikingNeuron {

	public:

		// Constructors and destructor
		VoltageTriggeredSpikingNeuron(Model* m = NULL);
		virtual ~VoltageTriggeredSpikingNeuron();

		// Return an estimate of instanteous firing rate based on spikes.
		// Default is an estimate based on a unitary alpha function kernel.
		// Subclasses may override to use a different estimator.
		virtual Number			firingRate();

		// Access time constants used in estimating current firing rate.
		virtual SimTime			firingRateTau()	{ return _firingRateTau; }
		virtual void			firingRateTau(SimTime tau);

		// Handle the end of the time step by checking for spikes
		virtual void			timeStepEnded();			

	protected:
		Number					_previousFiringVm;	// previous voltage
		Number					_previousFiringVmDot; // previous dv/dt

		SimTime					_firingRateTau;	// time constant for estimating rate
		Number					_ExpHFRTau;		// firing rate cache of exp(-h/tau)
		Number					_firingRateS1;	// state 1 for firing rate est
		Number					_firingRateS2;	// state 2 for firing rate est

		// Handle additional actions after a firing -- default is no-op
		virtual void			postFiringActions(SimTime tspike) {}

		// Update internal states when estimating firing rate.
		// spikeNow is true if this step has a spike.
		virtual void			updateFiringRate(bool spikeNow);

		// Subclass responsibilities (some optional)-------------------

		// Voltage to use for determining whether firing occurred
		virtual Number			firingVoltage() = 0;

		// Time derivative of firing voltage
		virtual Number			firingDeriv() = 0;

		// Voltage threshold which must be exceeded to be counted as a spike (default)
		virtual Number			firingThreshold()	{ return -30*UOM::mV; }

		// Absolute minimum time between spikes (default)
		// If firing voltage is persistently above threshold, then
		// generated spikes will be separated by this amount.
		virtual SimTime			minISI()			{ return 5*UOM::msec; }

	private:
		void					initialize();
	};


	// ----------------------------------------------------------------
	// CLASS:	MorphologyEntry	
	// EXTENDS:	none
	// DESC:	Public structure for storing one entry in
	//			an array describing compartents to use in
	//			a cell model. 
	//
	//			The soma is usually the first entry, which
	//			is numbered as 0. Parent refers to the compartment
	//			that is next most proximal to the soma.
	//
	//			Branch is an arbitrary identifier that can be used
	//			to group all compartments between common branch points.
	//
	// RESP:
	//		1.	Know compartment size 
	//		2.	Know distance from soma
	//		3.	Know parent compartment in dendrite tree
	//
	// NOTE:	Types codes are taken from .swc format.
	//			1 = soma, 2 = axon, 3=basal dendrite, and
	//			4 = apical dendrite. -1 is an end marker.
	//			The type is defined separately to simplify
	//			reading and writing table entries.
	//			Entries with parent=0 imply connection with
	//			the tree root (first morphology entry).
	//
	//			idnum is included to make entries self documenting.
	//			In the overall table, idnum should equal array index.
	//			branch is an aribtrary identifier of the branch
	//			on which a segment appears.
	//
	//			Orientation (dx, dy, dz) is the vector from the 
	//			proximal to distal ends of the compartment. Because
	//			a compartment does not exactly correspond with dendrite
	//			locations within the compartment, the orientation
	//			vector may of different length than the compartment.
	// ----------------------------------------------------------------

	class MorphologyEntry {
	public:

		enum swcType {
			somaType=1,
			axonType=2,
			basalDendriteType=3,
			apicalDendriteType=4,
			endType=-1};

		int				type;			// type of entry
		int				idnum;			// identifier of this entry (should = index)
		int				parent;			// identifier of parent entry in tree structure
		int				branch;			// arbitrary identifier of the branch
		double			r;				// radius in microns
		double			len;			// length in micron
		double			dist;			// distance from soma in microns
		double			x;				// X coord of this compartment
		double			y;				// Y coord of this compartment
		double			z;				// Z coord of this compartment
		double			dx;				// X coord of the orientation vector
		double			dy;				// Y coord of the orientation vector
		double			dz;				// Z coord of the orientation vector
	};

	// ----------------------------------------------------------------
	// CLASS:	MorphologicalNeuron
	// EXTENDS:	VoltageTriggeredSpikingNeuron
	// DESC:	Abstract class representing a spiking neuron
	//			defined in terms of a morphology table
	//			supplying compartment sizes and relationships.
	// RESP:
	//		1.	Know cell morphology
	//		2.	Allocate compartments
	//		3.	Electrically connect compartment per morphology
	//		4.	Locate compartments
	//		5.	Apply neuromodulation factors
	//
	// NOTE:	A morphology table is an array, usually static,
	//			of entries (see MorphologyTableEntry).
	//			The first entry is the root of the tree
	//			and is either the whole soma or a middle
	//			compartment of the soma. The table is
	//			assumed to be organized with all dendrite
	//			compartments contiguous. Axon compartments
	//			may or may not be included in the table.
	//			Because morphologies might not include
	//			axons, the only treatment here is to allocate
	//			compartments found in the table.
	// ----------------------------------------------------------------

	class MorphologicalNeuron : public VoltageTriggeredSpikingNeuron {

	public:

		// Constructors and destructor
		MorphologicalNeuron(Model* m = NULL);
		virtual ~MorphologicalNeuron();

		// Access morphology table
		inline MorphologyEntry* morphology() { return _morphology; }
		virtual void			morphology(MorphologyEntry* mt);

		// Flag indicating that this is a morphological neuron
		virtual bool			isMorphologicalNeuron() { return true; }
		
		// Access number of compartments of different types
		inline int				numSomaComp()		{ return _numSomaComp; }
		inline int				numDendriteComp()	{ return _numDendriteComp; }
		inline int				numAxonComp()		{ return _numAxonComp; }
		inline int				numISComp()			{ return _numISComp; }

		// Access a unit normal vector defining the laminar orientation
		inline Number			orientationX()		{ return _orientationX; }
		inline Number			orientationY()		{ return _orientationY; }
		inline Number			orientationZ()		{ return _orientationZ; }

		// Access a soma compartment by number starting with 1.
		virtual Compartment*	somaComp(int n=1);

		// Access axon compartment by number starting with 1
		virtual Compartment*	axonComp(int n=1);

		// Access initial segment compartment by number starting with 1
		virtual Compartment*	ISComp(int n=1);

		// Access dendrite compartment by number starting with 1.
		virtual Compartment*	dendriteComp(int n=1);

		// Get Y coordinate of a dendrite compartment
		virtual Number			dendriteYCoord(int n)
		{ return dendriteMorphologyEntry(n)->y*UOM::micron; }

		// Access a morphology entry by dendrite number starting with 1
		virtual MorphologyEntry* dendriteMorphologyEntry(int n);

		// Access a dendrite compartment by branch id. The dendrite
		// with the soma distance closest to the value provided is returned.
		// If the branch is not found, NULL is returned. Default distance
		// should always locate the most distal compartment in the branch.
		virtual Compartment*	dendriteCompByBranch(
			int					branchId,				// numerical id of branch
			Number				dist=UOM::meter);		// distance from soma

		// Set ion channel conductance multiplier in soma and all dendrites
		// for the identified ion channel to be the value provided.
		// mustMatch esures that either the soma or at least one compartment
		// constains the indicated channel.
		virtual void			setGModulator(
			TokenId				chanId,					// id of channel
			Number				value,					// value to set
			bool				mustMatch=true);		// must be matched somewhere

		virtual void			setGModulator(
			char*				chanName,				// id of channel as string
			Number				value,					// value to set
			bool				mustMatch=true)			// must be matched somewhere
		{	setGModulator(token(chanName),value,mustMatch); }

		// Set ion channel conductance multiplier for an ion channel
		// using a variant of the Michaelis-Menten formula: 
		// g_as_modulated = g_normal * (1+a*X/(X+Kd))
		virtual void			setMichaelisMentenMod(
			TokenId				chanId,					// id of channel
			Number				X,						// modulator concentration
			Number				amax,					// a in above formula
			Number				Kd,						// Kd in above formula
			bool				mustMatch=true);		// must be matched somewhere 

		virtual void			setMichaelisMentenMod(
			char*				chanName,				// id of channel as string
			Number				X,						// modulator concentration
			Number				amax,					// a in above formula
			Number				Kd,						// Kd in above formula
			bool				mustMatch=true)			// must be matched somewhere 
		{	setMichaelisMentenMod(token(chanName),X,amax,Kd,mustMatch); }

		// Test dendrite type based on number.
		virtual bool			isApicalDendrite(int n) { return !isBasalDendrite(n); }
		virtual bool			isBasalDendrite(int n);

	protected:
		MorphologyEntry*		_morphology;		// array of entries
		int						_morphologySize;	// number of entries

		// Define a unit normal vector that can be used to locate
		// different parts of a laminar structure. Default is (0,1,0).
		Number					_orientationX;
		Number					_orientationY;
		Number					_orientationZ;

		// Soma compartments are allocated contiguously.
		// By convention, soma compartments are allocated first.
		int						_numSomaComp;

		// Dendrite compartments are allocated contiguously 
		// and can thus be located among compartments by knowing
		// how many there are and where to start.
		int						_numDendriteComp;
		int						_dendriteCompOffset;

		// Axon compartments are allocated contiguously 
		// and can thus be located among compartments by knowing
		// how many there are and where to start. If there are no
		// axons, the number and offset are zero.
		int						_numAxonComp;
		int						_axonCompOffset;

		// Initial segments compartments are allocated contiguously 
		// and can thus be located among compartments by knowing
		// how many there are and where to start. If there are no
		// axons, the number and offset are zero.
		int						_numISComp;
		int						_ISCompOffset;

		// Return offset into morphology table for first dendrite entry.
		virtual int				dendriteMorphologyOffset();

		// Utility function(s) for use by subclasses  --------------------
		
		// Allocate the soma as a unipotential ball with membrane area the
		// same as the first morphology entry compartment.
		virtual void			createSoma(
			Number				cm_sp=1*UOM::microF/UOM::cm_2,		// Cm specific capacitance
			Number				rm_sp=50000*UOM::ohm*UOM::cm_2);	// Rm specific conductivity

		// Allocate and connect dendrite compartments using common defaults, 
		// i.e. soma is first compartment and also the first morphology entry. 
		// Soma and any axon compartments must be added before invoking this function. 
		// Distance from soma is automatically set after all dendrites are added.
		virtual void			createDendrites(
			Number				cm_sp=1*UOM::microF/UOM::cm_2,		// Cm_specific for all comp
			Number				ri_sp=200*UOM::microF*UOM::cm,		// Ri_specific for all comp
			Number				rm_sp=50000*UOM::ohm*UOM::cm_2,		// Rm_specific for all comp
			Number				areaAdjust=1);						// membrane area adjustment

		// Electrically connect compartments using common defaults. 
		// A single compartment soma is assumed. The relationship 
		// between dendrite morphology entries and comartments must
		// be established prior to invocation and all dendrites must
		// be contiguous in the compartments vector and morphology table.
		virtual void			connectDendrites();

		// Electrically connect compartments with full parameters.
		// The root compartment is assumed to be an equipotential ball.
		// Parentage outside the indicated range is ignored.
		virtual void			connectCompartments(
			int					fromIdx,			// starting morphology index
			int					toIdx,				// ending index
			int					offset,				// starting offset into compartments vector
			Compartment*		root);				// root compartment to connect to

	private:

		// Minimally initialize the instance during construction
		virtual void			initialize();
	};

	// ----------------------------------------------------------------
	// CLASS:	Compartment	
	// EXTENDS:	ModelComponent
	// DESC:	Model one compartment of in a compartmental
	//			neuron simulation. Default geometry is that
	//			of a cylinder.
	// RESP:
	//		1.	Know containing neuron
	//		2.	Know physical parameters of compartment
	//		3.	Know channels
	//		4.	Know membrane voltage
	//		5.	Accumulate ionic currents and compute dV/dt
	//		6.	Set and clear current injection
	//		7.	Set and clear voltage clamp
	//
	//NOTES:	The convention for unit of measure in which
	//			resistance and capacitance are scaled
	//			based on cell geometry is not used here.
	//
	//			Specific parameter values can be provided
	//			and these are scaled based on compartment
	//			size. The default geometry for a compartment
	//			is a cylinder, but subclasses can use other
	//			conventions if specific values are scaled. 
	//
	//			Values for Rm, Ri, Cm are absolute units
	//			(for example, ohms and farads) and should
	//			be converted to resistance or capacitance before 
	//			being set, for example, Ri(2*gigaohm) vs. 
	//			Ri_specific(100*ohm*cm). However, as a practical 
	//			matter, parameters are usually set via specific
	//			unit interfaces and the conversion is automatic.
	//
	//			Membrane area adjustment factor is a useful fudge
	//			factor to account for spines or non-cylindrical
	//			geometry of the compartment.
	// ----------------------------------------------------------------

	class Compartment : public ModelComponent {

	public:

		// Constructors and destructor
		Compartment(bool doInit=true);
		Compartment(
			Number r,								// radius
			Number len,								// length
			Number cm_sp=1*UOM::microF/UOM::cm_2,	// Cm_specific
			Number ri_sp=200*UOM::microF*UOM::cm,	// Ri_specific
			Number rm_sp=50000*UOM::ohm*UOM::cm_2,	// Rm_specific
			Number areaAdj=1);	// Membrane area adjustment factor

		virtual ~Compartment();

		// Model accessors
		inline Model*			model() {return _model; }
		virtual void			model(Model* m);

		// Neuron accessors
		inline Neuron*			neuron() { return _neuron; }
		virtual void			neuron(Neuron* n) { _neuron = n; }

		// Get and set membrane voltage
		inline Number			Vm() { return _Vm; }
		virtual void			Vm(Number v);
		inline int				VmIndex() { return _VmIndex; }
		inline Number			VmRem() { return _VmRem; }

		// Get derivative of the membrance voltage
		inline Number			VmDot() { return derivValue(0); }

		// Get/set this compartment's leak reversal potential.
		inline  Number			Vleak() { return _Vleak; };
		virtual void			Vleak(Number v) { _Vleak=v; }

		// Get/set this compartment's initial resting potential
		inline  Number			Vinit() { return _Vinit; };
		virtual void			Vinit(Number v) { _Vinit=v; }

		// Absolute RC values -- see note about units of measure
		inline Number			Cm() { return _Cm; }
		inline Number			Rm() { return _Rm; }
		inline Number			Ri() { return _Ri; }
		inline Number			gLeak() { return 1/_Rm; }

		// Compartment geometry accessors
		inline  Number			radius() { return _radius; }
		virtual void			radius(Number r);

		inline  Number			length() { return _length; }
		virtual void			length(Number len);

		inline  Number			areaAdjustment() { return _areaAdjustment; }
		virtual void			areaAdjustment(Number x);

		// Access the path distance of this compartment from the soma
		// or equivalent root compartment (soma may be multiple comp).
		inline  Number			distFromSoma() { return _distFromSoma; }
		virtual void			distFromSoma(Number d) { _distFromSoma = d; }

		// Specific capacitivity (absolute units, eg converted from microF/cm^2)
		virtual Number			Cm_specific() { return _Cm_specific; }
		virtual void			Cm_specific(Number c);

		// Specific membrane resistivity (SI units, eg converted from k-ohm*cm^2)
		virtual Number			Rm_specific() { return _Rm_specific; }
		virtual void			Rm_specific(Number r);

		// Specific internal resistivity (SI units, eg converted from k-ohm*cm)
		virtual Number			Ri_specific() { return _Ri_specific; }
		virtual void			Ri_specific(Number r);

		// Specfic leak conductivity (SI units, eg converted from mS/cm^2)
		virtual Number			gLeak_specific() { return 1/Rm_specific(); }
		virtual void			gLeak_specific(Number g) { Rm_specific(1/g); }
		
		// Compute geometric properties from dimensions
		virtual Number			membraneArea() { return _membraneArea; }
		virtual Number			crossSectionArea();
		virtual Number			volume();
		virtual Number			subshellVolume(Number depth);	// shell inside
		virtual Number			extraShellArea(Number dist);	// shell outside

		// Maintain ion channel relationships with the compartment
		virtual void			addIonChannel(IonChannel* ic);
		virtual void			removeIonChannel(IonChannel* ic);

		// Maintain electrical coupling relationships with the compartment
		virtual void			addCoupling(ElectricalCoupling* ec);
		virtual void			removeCoupling(ElectricalCoupling* ec);

		// Return all compartments coupled with this one
		virtual CompartmentVector connectedCompartments();

		// Maintain calcium pool relationships with the compartment
		virtual void			addCalciumPool(CalciumPool* pool);
		virtual void			removeCalciumPool(CalciumPool* pool);

		// Polymorphic add for convenience
		virtual void			add(IonChannel* ic) { addIonChannel(ic); }
		virtual void			add(ElectricalCoupling* ec) { addCoupling(ec); }
		virtual void			add(CalciumPool* pool) {addCalciumPool(pool); }

		// Polymorphic remove for convenience
		virtual void			remove(IonChannel* ic) { removeIonChannel(ic); }
		virtual void			remove(ElectricalCoupling* ec) { removeCoupling(ec); }
		virtual void			remove(CalciumPool* pool) {removeCalciumPool(pool); }

		// Get calcium channels associated with the compartment
		virtual IonChannelVector getCalciumChannels();

		// Get an ion channel by matching on the component name or token. The first 
		// matching channel is returned. No warning is provided if multiple matches 
		// are present (basically a bug but you were warned -- this can be fixed later).
		// If no channels match and mustMatch is true, a fatal error occurs. 
		// Otherwise if none match, NULL is returned.
		virtual IonChannel*		findIonChannel(TokenId channelId, bool mustMatch=true);
		virtual IonChannel*		findIonChannel(char* channelName, bool mustMatch=true);

		// Get a particular synaptic response by component name or token. The first
		// match is returned and must be present if mustMatch=true.
		virtual SynapticResponse* findSynapticResponse(TokenId respId, bool mustMatch=true);
		virtual SynapticResponse* findSynapticResponse(char* respName, bool mustMatch=true);

		// Create a synapse given the component name or token of an existing synaptic response.
		virtual Synapse*		createSynapse(
			TokenId					id,							// synapse component id
			AxonProcess*			axon = NULL,				// presynaptic axon process
			Number					wght = 1.0,					// initial weight value
			Number					dist = 100*UOM::micron);	// axonal distance

		virtual Synapse*		createSynapse(
			char*					name,						// synapse component name
			AxonProcess*			axon = NULL,				// presynaptic axon process
			Number					wght = 1.0,					// initial weight value
			Number					dist = 100*UOM::micron);	// axonal distance

		// Interfaces for experimental manipulations ------------------

		// Get and set injected current (inward negative)
		inline Number			Iinjected() { return _Iinjected; }
		virtual void			Iinjected(Number ic) { _Iinjected = ic; }
		virtual void			IinjectedSpecific(Number ic) { Iinjected(ic*membraneArea()); }

		// Set or clear a clamp at a fixed membrane potential
		virtual void			setVoltageClamp(Number vm);
		virtual void			clearVoltageClamp();
		inline  bool			isVoltageClamped() { return _isVoltageClamped; }

		// Reporting interfaces ---------------------------------------

		// Add to components to probe when reporting on
		// this compartment. By default this includes
		// ion channels and calcium pools as well as
		// the compartment itself.
		virtual void			addToComponentsToProbe(ModelComponentVector& comps);

		// Access the component name for external reporting
		virtual const char*		componentName();
		virtual void			componentName(char* nameStr);

		// Access the numeric identifier for external reporting
		virtual int				numericIdentifier() { return _numericIdentifier; }
		virtual void			numericIdentifier(int n) { _numericIdentifier=n; }

		// Get total current as of last computeDerivatives
		virtual Number			Itotal() { return _Itotal; }

		// Default Itotal unit of measure for reporting
		virtual Number			ItotalUnits() { return UOM::picoA; }

		// Default state label for reporting
		virtual const char** stateLabels() { 
			static const char* sl="Vm"; 
			return &sl; }

		// Default statee value units of measure for reporting
		virtual Number* unitsOfMeasure() {
			static Number units=UOM::mV; 
			return &units; }

		// ODE solver interface functions -----------------------------

		// Return the domain identifier for compartments
		virtual unsigned int	domain() { return 0x0010; }

		// Return number of state variables in this compartment
		virtual int				numStateVar() { return 1; }

		// Do early initializations before components set initial values
		virtual void			simulationStarted();

		// Set initial states
		virtual void			setInitialState();

		// Set default weight value for membrane voltage variable
		virtual void			setWeightValues();

		// Tell about this object if asked
		virtual bool			notifyOnStateVectorChanged() { return true; }
		virtual bool			isCompartment() { return true; }

		// Compartments are usually second order accurate
		virtual bool			isSecondOrderAccurate() { return true; }

		// Update cached values when state vector is changed
		virtual void			stateVectorChanged(); 

		// Compute derivatives of state variables
		virtual void			computeDerivatives();

		// Compute derivatives after the fact
		virtual void			recomputeDerivatives();

		// Access index of this compartment in a Jacobian matrix
		inline  int				jacobianIndex() { return _jacobianIndex; }
		virtual void			jacobianIndex(int n) { _jacobianIndex = n; }

		// Build Jacobian row for this compartment from scratch.
		// Compute derivatives must be done before invoking this.
		virtual void			buildJacobian(SparseMatrix& jacobian);

		// Update Jacobian matrix with current values following
		// the initial build. Unless cell geometry changes, this
		// would only involve changing terms on the matrix diagonal
		// associated with this compartment. Compute derivatives must 
		// be done before invoking this.
		virtual void			updateJacobian(SparseMatrix& jacobian);

		// Compute the leak current from the current potential
		virtual double			Ileak() { return (Vm()-Vleak())/Rm(); }

		// Compute the channel current from the current state
		virtual double			Ichan();

		// Compute the current from all electrical couplings
		virtual double			Icouple();

	protected:

		// Link back to containing neuron
		Neuron*					_neuron;

		// Name and id of this compartment for external reporting
		char*					_name;			// NULL if no value supplied
		int						_numericIdentifier;

		// Current state information cached here
		int						_VmIndex;		// Vm for table lookup
		Number					_Vm;			// membrane voltage
		Number					_VmRem;			// remainder of Vm/VStepForIndex;

		// Values saved as a side-effects during compute derivatives etc.
		// This method is used to avoid some redundant arithmetic.
		double					_gChanTotal;	// total ion channel conductance
		double					_Itotal;		// total save current

		// Structure information
		IonChannelVector			_ionChannels;	// ion channels used here
		ElectricalCouplingVector	_couplings;		// electrotonic connections
		CalciumPoolVector			_calciumPools;	// calcium pools used here

		// Absolute parameters
		Number					_radius;		// cylinder radius
		Number					_length;		// cylinder length
		Number					_Cm;			// membrane capacitance
		Number					_Rm;			// membrane resistance
		Number					_Ri;			// internal resistance
		Number					_areaAdjustment; // adjustment factor for area
		Number					_membraneArea;	// area of membrane surface
		Number					_distFromSoma;	// path distance from soma

		// Specific (scaled by membrane area) values of parameters
		Number					_Cm_specific;	// specific membrane capacitance
		Number					_Rm_specific;	// specific membrane resistivity
		Number					_Ri_specific;	// specific internal resistivity

		// Leak and resting potential for this compartment
		Number					_Vinit;			// initial value only
		Number					_Vleak;			// resting leak reversal
		Number					_Vclamp;		// voltage clamp value

		// Value of injected current taken with "inward negative" sign conventions,
		// that is, positive injected currents tend to depolarize the compartment.
		Number					_Iinjected;

		// Flags relating to voltage clamps
		bool					_isVoltageClamped;		// voltage clamp in effect
		bool					_forceJacobianRebuild;	// rebuild on next update

		// Value used in setting up the Jacobian matrix for ODE
		// methods using a Jacobian. 
		int						_jacobianIndex;

		// Initialize the instance (e.g. during construction)
		void					initialize();

		// Set default values for parameters during construction
		virtual void			setDefaults();

		// Adjust actual RC values based on specific ones
		virtual void			updateParams();

		// Propagate changes to RC values to any dependents
		virtual void			propagateUpdates();

		// Delete all held ion channels
		virtual void			deleteAllIonChannels();

		// Delete all associated couplings
		virtual void			deleteAllCouplings();

		// Delete all held calcium pools
		virtual void			deleteAllCalciumPools();
	};

	// ----------------------------------------------------------------
	// CLASS:	SphericalCompartment	
	// EXTENDS:	Compartment
	// DESC:	Model a single iopotential cell in the
	//			shape of a sphere.
	// RESP:
	//		1.	Know physical parameters of compartment
	//		2.	Scale parameters as appropriate to a ball.
	//
	// NOTES:	See note in Compartment concerning specific
	//			units of measure.
	// ----------------------------------------------------------------

	class SphericalCompartment : public Compartment {

	public:

		// Constructors and destructor
		SphericalCompartment(bool doInit=true);
		SphericalCompartment(
			Number r,								// radius
			Number cm_sp=1*UOM::microF/UOM::cm_2,	// Cm_specific
			Number rm_sp=50000*UOM::ohm*UOM::cm_2,	// Rm_specific
			Number spines=1);						// Spine area adjustment

		virtual ~SphericalCompartment();

		// Return derived geometric values
		virtual Number			volume();			// volume of the compartment
		virtual Number			subshellVolume(Number depth); // shell inside
		virtual Number			extraShellArea(Number dist); // shell outside

		// Adjust actual RC values based on specific ones
		virtual void			updateParams();

		// Ri is always zero for these compartments
		virtual Number			Ri() { return 0; }

		// Interfaces not supported for this object -------------------

		virtual Number			Ri_specific() { return 0; }
		virtual void			Ri_specific(Number r);
		virtual Number			length() { return 0; }
		virtual void			length(Number len);
		virtual Number			crossSectionArea();
	};

	// ----------------------------------------------------------------
	// CLASS:	ElectricalCoupling
	// EXTENDS:	none
	// DESC:	Represents electrotonic coupling between
	//			two compartments
	// RESP:
	//		1.	Know pair of related compartments
	//		2.	Know conductance between compartments
	//		3.	Compute conductance if requested
	//		4.	Compute current flow
	//
	// NOTE:	The default computation of conductance is
	//			based on the assumption that each compartment
	//			contributes half of its internal resistance.
	//			This may be reasonable for cylinders with
	//			connections on the ends but can be wrong in
	//			other cases. Subclasses may provide for
	//			alternative configurations.
	// ----------------------------------------------------------------

	class ElectricalCoupling {

	public:
		// Constructors and destructor
		ElectricalCoupling();
		ElectricalCoupling(Compartment* c1, Compartment* c2);
		ElectricalCoupling(Compartment* c1, Compartment* c2, Number gval);
		virtual ~ElectricalCoupling();

		// Accessors
		virtual Compartment*		comp1();
		virtual void				comp1(Compartment* c1);

		virtual Compartment*		comp2();
		virtual void				comp2(Compartment* c2);

		// Return the number of nodes in the coupling
		virtual int					nodeCount();

		// Return all nodes in the coupling
		virtual CompartmentVector	nodes() { return _nodes; }

		// Get/set coupling conductance. If set this way, g is not
		// recalculated if compartment dimensions change.
		virtual double				g() { return _g; }
		virtual void				g(double gval);

		// Get the conductance between two peer in the coupling.
		// Except for capacitance, this is the Jacobian entry for
		// the row corresponding to c1 and column corresponding to c2.
		virtual double				gForPair(Compartment* c1, Compartment* c2);

		// Compute outward current with respect to a compartment
		virtual double				Iec(Compartment* comp);

		// Given one compartment, answer the other of the pair
		virtual Compartment*		peer(Compartment* comp);

		// Compute conductance based on comparment internal resistance
		virtual void				updateParams();

	protected:
		CompartmentVector			_nodes;		// associated compartments
		double						_g;			// absolute conductance
		bool						_gIsPreset; // suppress g recalculation

	};

	// ----------------------------------------------------------------
	// CLASS:	ElectricalJunction
	// EXTENDS:	ElectricalCoupling
	// DESC:	Represents electrotonic coupling between
	//			two or more compartments
	// RESP:
	//		1.	Know related compartments
	//		2.	Compute current flows
	//
	// NOTE:	Compartments are assumed to be joined at the
	//			ends while their respective voltages are 
	//			determined in the middle of the compartment.
	//			Ri values in the compartments are used to
	//			extract conductances for the junction.
	//
	//			For this class, g() returns total conductance
	//			for all nodes in the junction.
	// ----------------------------------------------------------------

	class ElectricalJunction : public ElectricalCoupling {

	public:
		// Constructors and destructor
		ElectricalJunction();
		ElectricalJunction(Compartment* comp);

		virtual ~ElectricalJunction();

		// Accessors
		virtual double				g() { return _g; } // total conductance

		// Add/remove a node
		virtual void				add(Compartment* comp);
		virtual void				remove(Compartment* comp);

		// Compute effective conductances for the junction
		virtual void				updateParams();

		// Get the conductance between two peer in the coupling.
		// Except for capacitance, this is the Jacobian entry for
		// the row corresponding to c1 and column corresponding to c2.
		virtual double				gForPair(Compartment* c1, Compartment* c2);

		// Compute outward current with respect to one end of the junction
		virtual double				Iec(Compartment* comp);

		// If exactly two nodes in the junction, return the peer
		// Otherwise, return NULL
		virtual Compartment*		peer(Compartment* comp);

		// Functions that should not be used with junctions
		virtual void				g(double gval);
	};

	// ----------------------------------------------------------------
	// CLASS:	CalciumPool	
	// EXTENDS:	ModelComponent
	// DESC:	Abstract class for a pool model of
	//			calcium concentrations
	//			
	// RESP:
	//		1.	Know related channels providing Ca to pool
	//		2.	Know concentration of Ca ions (via subclass)
	//
	// NOTE:	Units of measure are dependent on the model. 
	//			Some models use dimensionless units, though
	//			agreement between channel and pool models is
	//			required to make such a scheme work.
	// ----------------------------------------------------------------

	class CalciumPool : public ModelComponent {

	public:

		// Constructors and destructor
		CalciumPool();
		virtual ~CalciumPool();

		// Accessors
		inline  Compartment*	container() { return _container; }
		virtual void			container(Compartment* comp);
		inline  Neuron*			neuron() { return container()->neuron(); }
		virtual bool			isCalciumPool() { return true; }

		// Return the domain identifier for calcium pools
		virtual unsigned int	domain() { return 0x0020; }

		// Add to a compartment as an ion channel
		virtual void			addTo(Compartment* comp);

		// Remove from the current compartment
		virtual void			removeFromCompartment();

		// Adjust for changes in compartment parameters (optional)
		virtual void			updateParams() {}

		// Access name for external reporting. Default is "CaPool"
		virtual const char*			componentName();
		virtual void			componentName(char* name);

		// Access accessible concentration of calcium ions in the pool
		virtual Number			CaX() = 0;		// subclass responsibility

		// Maintain relationship with channels supplying Ca++ ions
		virtual void			addSourceChannel(IonChannel* chan);
		virtual void			removeSourceChannel(IonChannel* chan);

		// Sum up the current from all channels supplying Ca++ ions
		virtual Number			ICa();

	protected:
		Compartment*			_container;			// containing compartment
		IonChannelVector		_sourceChannels;	// suppliers of Ca ions
		char*					_componentName;		// non-default component name

		// Add all calcium channels in the compartment
		virtual void			addAllCalciumChannels();
	};

	// ----------------------------------------------------------------
	// CLASS:	SimpleCalciumPool	
	// EXTENDS:	CalciumPool
	// DESC:	Represents a simple model of calcium dynamics
	//			for an active shell just inside the cell membrane.
	//			This shell can (and frequently will) be the
	//			entirety of a compartment.
	// RESP:
	//		1.	Add calcium channels in compartment
	//		2.	Calculate phi value (see below) from compartment size
	//		3.	Calculate beta (extrusion rate) using M-M pump dynamics
	//		4.	Calculate concentration derivative
	//
	// NOTES:	If no source channels are identified using
	//			addSourceChannel prior to starting the simulation,
	//			all calcium channels in the compartment will be added
	//			as sources. Otherwise, only those channels specifically 
	//			indicated as sources are used.
	//
	//			A typical unit here would be [Ca] in nano-molar
	//			(not to be confused with nanoMole). For models
	//			with non-dimensional calcium values, any consistent
	//			unit of measure can be used.
	//
	//			If we assume that a calcium current directly
	//			contributes to the pool concentration we can
	//			compute phi as:
	//
	//			phi = 1/(2*F*v) 
	//
	//			where F is Farady's constant, and v is pool volume.
	//			The factor of 2 is the charge for Ca++. 
	//
	//			The pool volume can be computed based on the assumption of 
	//			a cylindrical compartment with the pool as a shell
	//			inside the exterior of the cylinder. If shell depth
	//			exceeds the compartment radius, the whole compartment's
	//			volume is dedicated to the pool.
	//
	//			Typically we assume a buffering process in which
	//			most of the ions entering the pool are bound in a 
	//			high rate reaction such that there is a stable ratio 
	//			between unbound ions and bound ions. The value 
	//			unboundRatio reflects this and is used to adjust the
	//			computed value of phi in the formula above.
	//
	//			For the simplest model, the ODE for pool dynamics is:
	//
	//			dX/dt = - phi*Ica - beta*(X-Xrest)
	//
	//			where is X is the pool concentration,
	//			Xrest is an equilibrium concentration level, and
	//			Ica is the inward Ca++ current (and hence negative).
	//
	//			If parameters for pumping dynamics are specified,
	//			the terms X and Xrest are replaced respectively by:
	//
	//			X/(1+X/Kd) and Xrest/(1+Xrest/Kd)
	//
	//			where Kd is the concentration for half-maximal extrusion.
	//
	//			To express beta in the form of an inverse time constant,
	//			we can compute beta as follows:
	//
	//			beta=a/v*Vmax/Kd
	//
	//			where Vmax is the maximum pump rate, a is the membrane area,
	//			and v is the pool volume. 
	//
	//			The following arbitary defaults are used during object
	//			construction and can be changed later.
	//
	//			ODE error weight =		1/250 nano-molar^-1
	//			shell depth =			infinite
	//			unbound ratio =			0.001
	//			Vmax =					6e-14 mmol/msec/cm_2
	//			Kd =					1 micro-molar
	//			CaXrest =				50 nano-molar
	//			CaXinit =				not set (defaults to CaXrest)
	//			
	// ----------------------------------------------------------------

	class SimpleCalciumPool : public CalciumPool {

	public:

		// Constructors and destructor
		SimpleCalciumPool();				// init with defaults
		
		SimpleCalciumPool(					// init using specified values
			Number				rest,		// resting concentration
			Number				init,		// initial concentration
			Number				ubr,		// unbound Ca++ ratio
			Number				vmx,		// vmax pump rate
			Number				kd,			// Pump half rate concentration
			Number				shell=numeric_limits<Number>::infinity()); // shell depth
		
		SimpleCalciumPool(					// init with fixed phi and beta
			Number				phi,		// fixed phi value
			Number				beta,		// fixed beta value
			Number				rest=0);	// fixed CaXrest value, if any
											// CaXinit is set to -1 forcing use of CaXrest
		virtual ~SimpleCalciumPool();

		// Accessors
		inline  Number			phi() { return _phi; }
		inline  Number			beta() { return _beta; }

		inline  Number			shellDepth() { return _shellDepth; }
		inline  Number			unboundRatio() { return _unboundRatio; }
		inline  Number			vmax() { return _vmax; }
		inline  Number			Kd() { return _Kd; }

		inline  Number			CaXrest() { return _CaXrest; }
		virtual void			CaXrest(Number x);

		inline  Number			CaXinit() { return _CaXinit; }	// CaXinit<0 if not set
		virtual void			CaXinit(Number x) { _CaXinit = x; }

		inline  Number			weight() { return _weight; }
		virtual void			weight(Number w) { _weight = w; }
	
		// Current concentration of calcium in the pool
		virtual Number			CaX() { return stateValue(0); }

		// Specifying the following values overrides any calculated values.
		virtual void			phi(Number phiValue);
		virtual void			beta(Number betaRateValue);
		
		// Phi can be computed from compartment size by specifying
		// parameter values through the following functions.
		virtual void			shellDepth(Number d);
		virtual void			unboundRatio(Number ubr);

		// Pump dynamics can be specified to compute beta based on
		// current inside concentration and parameters provided.
		virtual void			vmax(Number r);		// max pump rate
		virtual void			Kd(Number kd);		// half rate concentration

		// ODE Interface functions
		virtual int				numStateVar() { return 1; }
		virtual void			simulationStarted();
		virtual void			setInitialState();
		virtual void			setWeightValues();
		virtual void			computeDerivatives();

		// Interfaces for implicit ODE solvers
		virtual bool			canPerformLocalUpdate() { return true; }
		virtual bool			isSecondOrderAccurate() { return true; }
		virtual void			localStateUpdate(SimTime h, CNStepType stepType);

		// Adjust for changes in compartment parameters
		virtual void			updateParams();
		inline  Compartment*	container() { return _container; }
		virtual void			container(Compartment* comp);

		// Reporting accessors
		virtual const char** stateLabels() {
			static const char* sl="cax"; 
			return &sl; }
		virtual Number* unitsOfMeasure() {
			static Number units = UOM::nanoM;
			return &units; }

	protected:

		// Parameter values for this pool
		Number					_phi;				// Mult for current to conc chg
		Number					_beta;				// Mult for conc decay from pool
		Number					_CaXrest;			// Resting Ca++ concentration
		Number					_CaXinit;			// Initial Ca++ concentration (<0 if not set)
		Number					_unboundRatio;		// Fraction of ions not immediately bound
		Number					_shellDepth;		// Shell depth in compartment
		Number					_vmax;				// Pump maximum rate per unit area
		Number					_Kd;				// Half activation concentration
		Number					_weight;			// ODE solver error weight

		// The following flags determine whether to (re)calculate values based
		// on compartment size parameters.
		bool					_calcPhi;
		bool					_calcBeta;

		// Resting effective Ca concentration adjusted for pump dynamics
		Number					_prest;
	};

	// ----------------------------------------------------------------
	// CLASS:	AxonProcess
	// EXTENDS:	none
	// DESC:	Abstract class for representing axonal
	//			connections.
	// RESP:
	//		1.	Know the identifier of associated neuron
	//		2.	Know associated synapses
	//		3.	Know axonal propragation rate
	//		4.	Pass action potentials to synapses for processing.
	//
	// NOTE:	Process in this context refers to the
	//			axonal extension of the cell. It is a
	//			term of biology and not to be confused
	//			with "process" in the computer sense.
	//
	//			This object does not refer directly to its
	//			owning neuron. This is done to allow the
	//			use of proxy objects in distributed support
	//			or replay of recorded spike trains derived
	//			from sources other than Neuron objects.
	// ----------------------------------------------------------------

	class AxonProcess {
	
	public:

		// Constructors and destructor
		AxonProcess(int id = -1);
		virtual ~AxonProcess();

		// Access the identifier of the associated neuron
		inline  int				neuronId() { return _neuronId; }
		virtual void			neuronId(int id) { _neuronId = id; }

		// Get and set AP propagation velocity. Default = 0.5 m/sec.
		inline  Number			propRate() { return _propRate; }
		virtual void			propRate(Number cv) { _propRate = cv; }

		// Maintain synapse relationships
		virtual void			addSynapse(Synapse* syn);
		virtual void			removeSynapse(Synapse* syn);

		// Signal an action potential event.
		// Firing rate is an estimated rate provided by the invoker.
		virtual void			signalSpikeEvent(
			SimTime				t,				// Time of spike origination
			SimTime				isi,			// ISI for this spike
			Number				firingRate);	// Estimated firing rate

	protected:
		int						_neuronId;		// Owning neuron numeric id
		Number					_propRate;		// Propagation rate
		Synapse*				_synapses;		// Header for list of synapses
	};

	// ----------------------------------------------------------------
	// CLASS:	GHKTableEntry	
	// EXTENDS:	none
	// DESC:	Public structure for storing recomputed
	//			values for GHK formulas
	//			
	// RESP:
	//		1.	Store GHK effective potential
	//		2.	Store GHK effective conductance
	//
	// NOTE:	Look-up is by membrane voltage similar to
	//			alpha-beta table in VoltageDepTabChannel.
	//			See Compartment::Vm(v) for setting of indexes
	//			used to access the table.
	// ----------------------------------------------------------------

	class GHKTableEntry {
	public:
		Number			Veff;			// effective potential
		Number			Geff;			// effective conductance
	};

	// ----------------------------------------------------------------
	// CLASS:	IonChannel	
	// EXTENDS:	ModelComponent
	// DESC:	Abstract class for different types of
	//			ion channels.
	//			
	// RESP:
	//		1.	Know containg compartment
	//		2.	Know number of dynamic state variables (subclass)
	//		3.	Know offset into state vector
	//		4.	Know temperature adjustment factors
	//		5.	Calculate current flow (subclass)
	//		6.	Know default temperature for the simulation
	//		7.	Calculate Ca++ current via GHK equations for subclasses
	//
	// NOTES:	Conductance can be specified in either absolute units,
	//			i.e. direct electrical units such as ohms, or in specific 
	//			units, i.e. units scaled by membrane area such as ohms/cm^2.
	//			For conductance to be adjusted following size changes
	//			the specific form must be the last specification provided.
	// ----------------------------------------------------------------

	class IonChannel : public ModelComponent {

	public:

		// Constructors and destructor
		IonChannel(Number gSp=0);
		virtual ~IonChannel();

		// Accessors
		virtual TokenId			componentId();

		inline  Compartment*	container() { return _container; }
		virtual void			container(Compartment* comp);

		inline  CalciumPool*	calciumPool() { return _calciumPool; }
		virtual void			calciumPool(CalciumPool* pool);

		inline  Neuron*			neuron() { return container()->neuron(); }

		// Read only accessor for current channel or gate value
		inline  Number			value() { return stateValue(0); }

		// Return the domain identifier for ion channels
		virtual unsigned int	domain() { return 0x0030; }

		// Add to a compartment as an ion channel
		virtual void			addTo(Compartment* comp);

		// Remove from the current compartment
		virtual void			removeFromCompartment();

		// Set container without notifying containing compartment
		virtual void			setContainer(Compartment* comp) { _container = comp;}

		// Get the absolute conductance of the ion channel within its compartment
		inline Number			g() { return _g; }
		
		// Access the conductance parameter scaled by membrane area
		virtual Number			gSpecific();
		virtual void			gSpecific(Number gsp);
		virtual bool			gIsSpecific() { return _gIsSpecific; }

		// Access the conductance parameter in absolute terms (not scaled by area)
		virtual Number			gAbsolute();
		virtual void			gAbsolute(Number gabs);
		virtual bool			gIsAbsolute() { return !_gIsSpecific; }

		// Access the general conductance multiplier form of neuromodulation
		inline  Number			gModulator() { return _gModulator; }
		virtual void			gModulator(Number x) { _gModulator=x; }

		// Get membrane voltage and derivative
		inline Number			Vm() { return container()->Vm(); }
		inline Number			VmDot() { return container()->VmDot(); }

		// Get calcium concentration of associated calcium pool
		inline Number			CaX() { return calciumPool()->CaX(); }

		// Provide a common unit for measure for reporting currents
		virtual Number			IionUnits() { return UOM::picoA; }

		// Get membrane area (from compartment)
		virtual Number			membraneArea();

		// Update the conductance from current parameter settings.
		// This can be invoked by the containing compartment upon
		// changes in compartment size.
		virtual void			updateParams();

		// Utility functions for temperature adjustments --------------
		
		// Accessors for global defaults
		inline static Number	defaultTempC() { return _DefaultTempC; }
		static void				defaultTempC(Number temp);

		// Constants used in temperature adjustment
		static  Number			FoverRT(Number tempC);	// F/RT arbitrary temp (C)
		virtual Number			FoverRT();				// F/RT at default temp

		// Return a temperature adjustment rate of the form
		// q10^((T2-T1)/10), where
		// T1 is the rated temperature for the channel and
		// T2 is the current simulated temperature.
		// Once computed, the value is cached for faster access.
		virtual Number			Q10Factor();

		// ODE Solver interface functions -----------------------------

		// Indicate that this is an ion channel object
		virtual bool			isIonChannel() { return true; }

		// Indicate that for Crank-Nicolson, this is offset
		// from other objects in terms of evaluation time.
		virtual bool			isOffsetForCN() { return true; }
		
		// Initialize when simulation started. Note that
		// multiple subclasses can invoke this if needed.
		virtual void			simulationStarted();

		// Do any final clean-up operations at the end
		virtual void			simulationEnded();
		
		// Subclass responsibilties/options ---------------------------

		// Get net conductance for channel (pure function)
		virtual Number			conductance() = 0;

		// Apply neuromodulation parameters. Default action handles
		// only gModulator forms of modulation (id=token("gModulator"))
		// setModParam forms are convenience protocols for setModParams.
		virtual void			setModParams(
			TokenId				modId,			// modulation type token id
			int					numValues,		// number of parameter values
			Number*				values);		// parameter values as an array

		virtual void			setModParam(
			TokenId				modId,			// modulation type token id
			Number				value);			// parameter value

		virtual void			setModParam(
			char*				modTypeName,	// modulation type as string
			Number				value);			// parameter value

		// Locate the identified component within this ion channel.
		// If the component is not found, return NULL.
		// Default is to return this object or NULL.
		// Subclasses should override if they contain other components.
		virtual IonChannel*		findIonChannel(TokenId compId);

		// Compute current flow (Ohm's law is default)
		virtual Number			Iion() { return conductance()*(Vm()-Vrev()); }

		// Get both current and conductance at once.
		// Subclasses may override this to improve efficiency.
		virtual void			condAndIion(
								Number&	Gout,			// conductance
								Number& Iout)			// current
								{	
									Gout = conductance();
									Iout = Iion();
								}

		// Return reversal potential (pure function)
		virtual Number			Vrev() = 0;

		// Indicate that this channel admits no calcium ions (default)
		virtual bool			isCalciumChannel() { return false; }

		// Get current flow from calcium ion channels.
		// By default this is the same as Iion() but may be overridden.
		virtual Number			ICa();

		// Return current preparation temperature
		virtual Number			currentTempC() { return defaultTempC(); }

		// Return temperature at which rates were computed.
		// This will typically be overridden by a subclass.
		virtual Number			ratedTempC() { return 22; }	// i.e. room temp

		// Return a Q10 value for temperature adjustments to channel kinetics
		// This will typically be overridden by a subclass.
		virtual Number			Q10() { return 1; }	// i.e. no Q10 value set

		// Provide default names for lookup and external reporting.
		// Subclasses should override to provide meaningful names.
		virtual const char*		componentName() { return "IonChannel"; }
		virtual const char**	stateLabels() { // default labels
			static const char* slbl[] = {"s0","s1","s2","s3","s4","s5","s6","s7"}; 
			return slbl; }

		// Utility functions for calcium channels ---------------------

		// Accessors of values determined by GHK formula
		virtual Number			defaultPeakCaXin() { return _DefaultPeakCaXin; }
		virtual void			defaultPeakCaXin(Number x);

		virtual Number			defaultCaXout() { return _DefaultCaXout; }
		virtual void			defaultCaXout(Number x);

		// GHK functions come in two forms: first where all
		// variables must be supplied and second
		// where defaults are used. The default case can be
		// accomplished much faster through table look-up.

		// Default temperature is the current temperature 
		// as specified for this class. If this or other defaults
		// differ in subclasses, a separate GHK table can be
		// established in the subclass.

		// Attempts at accuracy in evaluating these functions
		// should not be taken as an endorsement of the accuracy
		// of an idealized GHK model to real calcium channels.
		// In practice, GHK equations are highly idealized and
		// are not necessarily good fits for any particular channel.

		// Scale GHK equation suitable for multiplication
		// by a "conductance" to yield current (ala Jaffe, 1994), 
		// which gives this units of voltage.  The formula is:
		// Vca = v*(1-CaIn/CaOut*exp(2vF/RT))/(1-exp(2vF/RT))
		static  Number			ghkCaEffectivePotential(
				Number			v,			// membrane voltage
				Number			CaXin,		// internal [Ca++]
				Number			CaXout,		// external [Ca++]
				Number			tempC);		// temp in centigrade

		// Same thing using default CaXin, CaXout, and tempC.
		// Note that CaXin/CaXOut is typically so low in mammals 
		// that the actual value of CaXin is irrelevant compared
		// to other sources of error in computing currents. V is
		// from the compartment associated with the channel.
		virtual Number			ghkCaEffectivePotential();
		
		// Return the derivative of the Ca effective potential
		// with respect to voltage. This gives a value which when
		// multiplied by a constant "conductance" is a linearized 
		// approximate conductance (e.g. for building a Jacobian).
		static  Number			ghkCaEffectiveCond(
				Number			v,			// membrane voltage
				Number			CaXin,		// internal [Ca++]
				Number			CaXout,		// external [Ca++]
				Number			tempC);		// temp in centigrade

		// Same thing using defaults. V is taken from the
		// compartment associated with the channel.
		virtual Number			ghkCaEffectiveCond();

		// Force (re)loading of GHK tables used in lookups
		virtual void			loadCaGHKTable();

	protected:
		TokenId					_componentId;		// token id from component name
		Compartment*			_container;			// containing compartment
		CalciumPool*			_calciumPool;		// pool for modulating Ca, if any
		Number					_gModulator;		// modulatory conductance multiplier
		Number					_g;					// conductance value for currents
		Number					_gParam;			// conductance param (abs or specific)
		bool					_gIsSpecific;		// indicates if _gParam is per unit area

		// Cached instance values
		Number					_cachedQ10Factor;	// Q10 multiplier at current temp

		// Static global parameters
		static Number			_DefaultTempC;		// default temperature in celcius
		static Number			_FoverRT;			// F/RT using default temp
		static Number			_DefaultCaXout;		// default external Ca++ concentration
		static Number			_DefaultPeakCaXin;	// default peak internal Ca++ conc.

		// Accessor for GHK table (can be overridden by subclasses)
		virtual GHKTableEntry**	pCaGHKTable() { return &_CaGHKTable; }

	private:
		static GHKTableEntry*	_CaGHKTable;		// look-up table for Ca++ GHK formulas
	};

	// ----------------------------------------------------------------
	// CLASS:	DummyIonChannel	
	// EXTENDS:	IonChannel
	// DESC:	Provide an ion channel object that does not
	//			do anything except hold a slot in the array
	//			of ion channels within a compartment.
	//			
	// RESP:
	//		1.	Handle manadatory interface calls.
	//		2.	Know number of state variables to reserve.
	//		3.	Express a dependency on calcium pool, if any.
	//
	// NOTE:	This object provides a handy way to make the
	//			number of state vector entries per compartment 
	//			constant so that externally written state values 
	//			can be more easily interpreted for visualization.
	// ----------------------------------------------------------------

	class DummyIonChannel : public IonChannel {

	public:

		// Constructors and destructor
		DummyIonChannel(int numStateVar=0,CalciumPool* pool=NULL);
		virtual ~DummyIonChannel();

		// Interface functions
		virtual int				numStateVar() { return _numStateVar; }
		virtual void			setInitialState() {}
		virtual void			computeDerivatives() {}
		virtual bool			canPerformLocalUpdate() { return true; }
		virtual void			localStateUpdate(SimTime h, CNStepType stepType) {}
		virtual const char*		componentName() { return "DummyIonChannel"; }
		virtual Number			conductance() {return 0; }
		virtual Number			Iion() {return 0; }
		virtual void			currentAndConductance(Number& Gout, Number& Iout)
								{ Gout = 0; Iout = 0; }

	protected:
		int						_numStateVar;

		// Return reversal potential
		virtual Number			Vrev() {return 0; }
	};

	// ----------------------------------------------------------------
	// CLASS:	MultiGateIonChannel	
	// EXTENDS:	IonChannel
	// DESC:	Abstract class for ion channels made
	//			of multiple gate variables each of which
	//			is a type of ion channel (for inheritance)
	// RESP:
	//		1.	Know gates
	//		2.	Setup relationship between model and gates
	//		3.	Make gates aware of containing compartment
	//		4.	Ensure that exactly the required number of 
	//			gates are provided by the time simulation starts.
	// ----------------------------------------------------------------

	class MultiGateIonChannel : public IonChannel {

	public:
		// Constructors and destructor
		MultiGateIonChannel(Number gSpVal=0);
		virtual ~MultiGateIonChannel();

		// Accessors
		inline  Compartment* container() { return _container; }
		virtual void container(Compartment* comp);

		inline  Model* model() { return _model; }
		virtual void model(Model* m);

		// Typically all state is in the gates themselves
		virtual int numStateVar() { return 0; } // default

		// Polymorphic add function definition
		virtual void add(IonChannel* gate) { addGate(gate); }

		// Make gate variables known
		virtual void addGate(IonChannel* gate);

		// Access a gate from the array of gates starting with 0
		virtual IonChannel* getGate(int n) {return _gates[n]; }

		// Add to components to probe when reporting.
		// This includes the current object and any members.
		virtual void addToComponentsToProbe(ModelComponentVector& comps);

		// Locate the identified component within this ion channel.
		// If the component is not found, return NULL.
		virtual IonChannel*		findIonChannel(TokenId compId);

		// ODE Solver interfaces --------------------------------------

		// Typically these channels will not have state of their own.
		virtual void setInitialState() {}

		// Derivatives are typically computed by individual
		// gate objects by virtue of gates being model components,
		// in which case this definition is a placeholder only.
		virtual void computeDerivatives() {}

		// Verify that the correct number of gates are available
		// at the start of simulation
		virtual void simulationStarted();

		// Subclass responsibilities/options --------------------------

		// Get net conductance for channel
		virtual Number			conductance() = 0;

		// Get number of gates required, which must match those
		// present when simulation started.
		virtual int				requiredNumGates() = 0;

		// Provide a component name. This is taken from the first
		// gate if not overriden.
		virtual const char*		componentName();


	protected:
		IonChannelVector		_gates;

	};
	
	// ----------------------------------------------------------------
	// CLASS:	MnHIonChannel	
	// EXTENDS:	MultiGateIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m^n*h.
	// RESP:
	//		1.	Supply required  number of gates
	//		2.	Compute conductance (via subclass)
	//		3.	Compute currents (via subclass)
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	// ----------------------------------------------------------------

	class MnHIonChannel : public MultiGateIonChannel {

	public:
		// Constructors and destructor
		MnHIonChannel(Number gSpVal=0) : MultiGateIonChannel(gSpVal) {}
		virtual ~MnHIonChannel() {}

		// Indicate number of gates required
		virtual int				requiredNumGates() { return 2; }
	};

	// ----------------------------------------------------------------
	// CLASS:	M1HIonChannel	
	// EXTENDS:	MnHIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	// ----------------------------------------------------------------

	class M1HIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M1HIonChannel(Number gSpVal=0);
		virtual ~M1HIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M2HIonChannel	
	// EXTENDS:	MnHIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	// ----------------------------------------------------------------

	class M2HIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M2HIonChannel(Number gSpVal=0);
		virtual ~M2HIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M3HIonChannel	
	// EXTENDS:	MnHIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	// ----------------------------------------------------------------

	class M3HIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M3HIonChannel(Number gSpVal=0);
		virtual ~M3HIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M4HIonChannel	
	// EXTENDS:	MnHIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*m*m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	// ----------------------------------------------------------------

	class M4HIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M4HIonChannel(Number gSpVal=0);
		virtual ~M4HIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M1HCaIonChannel	
	// EXTENDS:	MultiGateIonChannel
	// DESC:	Abstract class for multigate Ca++ ion channels 
	//			in which the formula for conductance
	//			is of the form m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	//			GHK equations are used to compute currents.
	// ----------------------------------------------------------------

	class M1HCaIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M1HCaIonChannel(Number gSpVal=0);
		virtual ~M1HCaIonChannel();

		// Indicate that this channel is a source of calcium
		virtual bool isCalciumChannel() { return true; }

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();

		// Required function -- unused but typical of Nernst potential
		virtual Number Vrev() { return 140*UOM::mV; }
	};

	// ----------------------------------------------------------------
	// CLASS:	M2HCaIonChannel	
	// EXTENDS:	MultiGateIonChannel
	// DESC:	Abstract class for multigate Ca++ ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	//			GHK equations are used to compute currents.
	// ----------------------------------------------------------------

	class M2HCaIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M2HCaIonChannel(Number gSpVal=0);
		virtual ~M2HCaIonChannel();

		// Indicate that this channel is a source of calcium
		virtual bool isCalciumChannel() { return true; }

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();

		// Required function -- unused but typical of Nernst potential
		virtual Number Vrev() { return 140*UOM::mV; }
	};

	// ----------------------------------------------------------------
	// CLASS:	M3HCaIonChannel	
	// EXTENDS:	MultiGateIonChannel
	// DESC:	Abstract class for multigate Ca++ ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	//
	// NOTE:	m is taken as the first gate added and h is the
	//			second added. There should be exactly two gates.
	//			GHK equations are used to compute currents.
	// ----------------------------------------------------------------

	class M3HCaIonChannel : public MnHIonChannel {

	public:
		// Constructors and destructor
		M3HCaIonChannel(Number gSpVal=0);
		virtual ~M3HCaIonChannel();

		// Indicate that this channel is a source of calcium
		virtual bool isCalciumChannel() { return true; }

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();

		// Required function -- unused but typical of Nernst potential
		virtual Number Vrev() { return 140*UOM::mV; }
	};

	// ----------------------------------------------------------------
	// CLASS:	MnHSIonChannel	
	// EXTENDS:	MultiGateIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m^n*h*s.
	// RESP:
	//		1.	Supply number of required gates
	//		2.	Compute conductance (via subclass)
	//		3.	Compute currents (via subclass)
	//
	// NOTE:	m is taken as the first gate added, h is the
	//			second, and s the third. Typically s is a slow
	//			inactivation gate.
	// ----------------------------------------------------------------

	class MnHSIonChannel : public MultiGateIonChannel {

	public:
		// Constructors and destructor
		MnHSIonChannel(Number gSpVal=0) : MultiGateIonChannel(gSpVal) {}
		virtual ~MnHSIonChannel() {}

		// Indicate number of gates required
		virtual int				requiredNumGates() { return 3; }
	};

	// ----------------------------------------------------------------
	// CLASS:	M1HSIonChannel	
	// EXTENDS:	MnHSIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*h*s.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	// ----------------------------------------------------------------

	class M1HSIonChannel : public MnHSIonChannel {

	public:
		// Constructors and destructor
		M1HSIonChannel(Number gSpVal=0);
		virtual ~M1HSIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M2HSIonChannel	
	// EXTENDS:	MnHSIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*h*s.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	// ----------------------------------------------------------------

	class M2HSIonChannel : public MnHSIonChannel {

	public:
		// Constructors and destructor
		M2HSIonChannel(Number gSpVal=0);
		virtual ~M2HSIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M3HSIonChannel	
	// EXTENDS:	MnHSIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*m*h*s.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	// ----------------------------------------------------------------

	class M3HSIonChannel : public MnHSIonChannel {

	public:
		// Constructors and destructor
		M3HSIonChannel(Number gSpVal=0);
		virtual ~M3HSIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	M4HSIonChannel	
	// EXTENDS:	MnHSIonChannel
	// DESC:	Abstract class for multigate ion channels 
	//			in which the formula for conductance
	//			is of the form m*m*m*m*h.
	// RESP:
	//		1.	Compute conductance
	//		2.	Compute currents
	// ----------------------------------------------------------------

	class M4HSIonChannel : public MnHSIonChannel {

	public:
		// Constructors and destructor
		M4HSIonChannel(Number gSpVal=0);
		virtual ~M4HSIonChannel();

		// Get net conductance and current for channel
		virtual Number			conductance();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	HHIonChannel	
	// EXTENDS:	ModelComponent
	// DESC:	Abstract class for channels with Markov
	//			Hodgkin-Huxley dynamics. A single
	//			state variable is the default.
	//			
	// RESP:
	//		1.	Define alpha and beta functions
	//		2.	Compute initial state
	//		3.	Compute derivatives
	//		4.	Update state using local implicit rule
	// ----------------------------------------------------------------

	class HHIonChannel : public IonChannel {

	public:
		// Constructors and destructor
		HHIonChannel(Number gSpVal=0);
		virtual ~HHIonChannel();

		// Optional alpha/beta functions for Hodgkin-Huxley dynamics 
		virtual Number		alpha() { return 0; }
		virtual Number		beta() { return 1; }

		// Default is one state variable for the gate
		virtual int numStateVar() { return 1; }

		// ODESolver interface functions ------------------------------

		// Set intial state based on t=infinity limit of alpha and beta
		virtual void		setInitialState();

		// Calculate time derivatives of associated state variables
		virtual void		computeDerivatives();

		// Update state by time step h using local implicit rule
		virtual void		localStateUpdate(SimTime h, CNStepType stepType);
		virtual bool		canPerformLocalUpdate() { return true; }
		virtual bool		isSecondOrderAccurate() { return true; }

		// Debug and reporting interface functions --------------------

		// Print the alpha and beta values for valid voltages.
		// This requires using simulation interfaces and should
		// be done for an object in an ongoing simulation.
		virtual void		printAlphaBeta(char* pathName=NULL);

	protected:

		// Utility functions for subclasses ---------------------------

		// Function for computing rates of the form r*v/(1-exp(-v/s))
		// using the value of the limit as v->0.
		double				linoidRate(double r, double v, double s);
	};

	// ----------------------------------------------------------------
	// CLASS:	HHIonGate
	// EXTENDS:	HHIonChannel
	// DESC:	Abstract class for single gate variables
	// RESP:
	//		1.	Allow definition of gates by providing defaults
	//			for pure functions used in channels.
	//
	// NOTE:	The class hierarchy isn't an "is a" pattern.
	// ----------------------------------------------------------------

	class HHIonGate : public HHIonChannel {

	public:
		// Constructors and destructor
		HHIonGate() {}
		virtual ~HHIonGate() {}

		// Accessors -- suppress side-effects of setting container
		inline Compartment*		container() { return _container; }
		virtual void			container(Compartment* comp) { _container = comp; }

		// Dummy functions (unused but required for instantiation)
		virtual Number		Vrev() { return 0; }
		virtual Number		conductance() { return 0; }

		// ODE Solver interface functions -----------------------------

		virtual bool		isIonChannel() { return false; }
		virtual bool		isGateVariable() { return true; }
	};

	// ----------------------------------------------------------------
	// CLASS:	BlendedIonChannel
	// EXTENDS:	HHIonChanell
	// DESC:	Abstract class for single gate channels
	//			with response consisting of a combination
	//			of two other single gate channels.
	// RESP:
	//		1.	Know constituent channels/gates
	//		2.	Connect constituents with containing compartment
	//		3.	Notify constituents on simulation start and end.
	//		4.	Compute blended alpha-beta values.
	//
	// NOTES:	Blending is done by linearly combining alpha
	//			and beta values from the constituent channels.
	//			This approximates a combination of different
	//			channel types with one state variable.
	//
	//			Constituent channel objects are owned by the
	//			blend object and are deleted when the blend
	//			object is deleted.
	//
	//			For gates g1, g2 and ratio r, the blended alpha
	//			and beta values are: (1-r)*alpha1+r*alpha2 and
	//			(1-r)*beta1+r*beta2.
	// ----------------------------------------------------------------

	class BlendedIonChannel : public HHIonChannel {

	public:

		// Constructors and destructor
		BlendedIonChannel(
			Number gSpVal=0,
			HHIonChannel* g1=NULL,
			HHIonChannel* g2=NULL,
			Number ratio=0.5);
		virtual ~BlendedIonChannel();

		// Accessors
		inline	Number			blendRatio() { return _blendRatio; }
		virtual void			blendRatio(Number ratio);

		inline	HHIonChannel*	gate1() { return _gate1; }
		virtual void			gate1(HHIonChannel* g1);

		inline	HHIonChannel*	gate2() { return _gate2; }
		virtual void			gate2(HHIonChannel* g2);

		inline	Compartment*	container() { return _container; }
		virtual void			container(Compartment* comp);

		// Blended alpha-beta values and derivatives
		virtual Number			alpha();
		virtual Number			beta();

		// Blended xinf-tau values
		virtual Number			xinf() { return 1/(1+beta()/alpha()); }
		virtual Number			tau() { return 1/(alpha()+beta()); }

		// ODE solver interface functions -----------------------------
		virtual void			simulationStarted();
		virtual void			simulationEnded();

	protected:
		HHIonChannel*			_gate1;
		HHIonChannel*			_gate2;
		Number					_blendRatio;
	};

	// ----------------------------------------------------------------
	// CLASS:	BlendedIonGate
	// EXTENDS:	HHIonGate
	// DESC:	Abstract class for a single gate object
	//			with aresponse consisting of a combination
	//			of two other gates. This is an specialization
	//			of BlendedIonChannel to gate objects.
	// RESP:
	//		1.	Supress functions unneeded for gates
	//
	// ----------------------------------------------------------------

	class BlendedIonGate : public BlendedIonChannel {

	public:

		// Constructors and destructor
		BlendedIonGate(
			HHIonChannel* g1=NULL,
			HHIonChannel* g2=NULL,
			Number ratio=0.5) 
			: BlendedIonChannel(0,g1,g2,ratio) {}

		virtual ~BlendedIonGate() {}

		// Dummy functions (unused but required for instantiation)
		virtual Number		Vrev() { return 0; }
		virtual Number		conductance() { return 0; }

		// ODE Solver interface functions -----------------------------

		virtual bool		isIonChannel() { return false; }
		virtual bool		isGateVariable() { return true; }
	};

	// ----------------------------------------------------------------
	// CLASS:	Order1BlendedIonChannel
	// EXTENDS:	BlendedIonChannel
	// DESC:	Abstract class for blended ion channels
	//			having a single variable with exponent 1.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// ----------------------------------------------------------------

	class Order1BlendedIonChannel : public BlendedIonChannel {

	public:

		// Constructors and destructor
		Order1BlendedIonChannel(
			Number gSpVal=0,
			HHIonChannel* g1=NULL,
			HHIonChannel* g2=NULL,
			Number ratio=0.5);

		virtual ~Order1BlendedIonChannel();

		// Get conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	AlphaBetaEntry	
	// EXTENDS:	none
	// DESC:	Public structure for storing recomputed
	//			alpha and beta values and their equivalent
	//			xinf and tau values.
	//			
	// RESP:
	//		1.	Know alpha and beta values
	//		2.	Know xinf and tau values
	//
	// NOTE:	This is a structure used in VoltageDepTabChannel
	//			subclasses. It is placed here so that the
	//			declaration is known in those classes.
	// ----------------------------------------------------------------

	class AlphaBetaEntry {
	public:
		Number			alphaValue;		// alpha value at corresponding voltage
		Number			betaValue;		// beta value at corresponding voltage
		Number			xinfValue;		// state value at t=infinity
		Number			tauValue;		// tau value
	};
	
	// ----------------------------------------------------------------
	// CLASS:	VoltageDepTabChannel
	// EXTENDS:	VoltageDepIonChannel
	// DESC:	Abstract class where table lookup is used
	//			for determining state transition rates based
	//			on compartment membrane voltage.
	// RESP:
	//		1.	Calculate table index from voltage (static)
	//		2.	Load table of alpha/beta values
	//		3.	Interconvert alpha/beta and xinf/tau forms
	//		4.	Fastpath implimentations of frequently invoked
	//			inherited functions.
	// ----------------------------------------------------------------

	class VoltageDepTabChannel : public HHIonChannel {

	public:

		// Constructors and destructor
		VoltageDepTabChannel(Number g=0);
		virtual ~VoltageDepTabChannel();

		// Accessors
		virtual int			numStateVar() { return 1; }		// 1 is a default value
		virtual bool		isVoltageDepTabChannel() { return true; }

		// Accessors for alpha-beta values and derivatives from table
		virtual Number		alpha(); 
		virtual Number		beta();

		// Accessors for xinf (state at t=inf value) and time constant from the table
		virtual Number		xinf();
		virtual Number		tau();
		
		// Allocate and load alpha-beta tables as needed.
		virtual void		loadAlphaBetaTable();
		virtual bool		alphaBetaTableLoaded();
		virtual void		deleteAlphaBetaTable();
		virtual bool		alphaBetaLoadInProg();

		// Subclass functions to compute alpha-beta table entries -----
		
		// Depending on the channel definition, either alpha and
		// beta values can be provided directly or else computed
		// from the t->infinity value of the state variable (xinf)
		// and time constant (tau).

		// Subclasses can either override both functions alphaForTable 
		// and betaForTable or else override functions xinfForTable and
		// tauForTable. If neither set of functions is used, 
		// loadAlphaBetaTable must be overridden to load the table.
		
		virtual Number alphaForTable(Number vm) { return -1; }	// -1 = no value
		virtual Number betaForTable(Number vm) { return 0; }

		virtual Number xinfForTable(Number vm) { return -1; }	// -1 = no value
		virtual Number tauForTable(Number vm) { return 0; }

		// Adjustments applied during table loading -------------------

		// Accessor for an exponent that is applied to xinfForTable as it is
		// loaded into the alpha-beta table. This can be used to express the
		// voltage response of multiple independent gates without affecting
		// the net time constants for activation and deactivation. 
		// Non-integer values are supported. xinfExponent must be non-negative.
		virtual Number		xinfExponent() { return 1; }

		// Accessors for a sigmoidal function by which xinfForTable is
		// multiplied before storage in the table. xinf(v) is multiplied by
		// 1/(1+exp(-(v-xinfMultVhalf)/xinfMultK)). xinfMultK()==0 suppresses
		// any multiplication of xinf since 0 is an impossible value. This
		// adjustment to xinf is made after applying any xinfExponent value.
		virtual Number		xinfMultVhalf() { return VMinForIndex; }
		virtual Number		xinfMultK() { return 0; }		// 0 = no value

		// Accessor for an absolute minimum value of xinfForTable. This can be
		// used when modeling an inactivation variable that is bounded away
		// from zero (leaky gate closure). The minimum is only imposed
		// when xinf-tau accessors are supplied by the subclass and is applied
		// after the other adjustments.
		virtual Number		xinfAbsMin() { return 0; }

		// Accessor for an absolute minimum value of tau. This may be
		// meeded when impossibly low tau values force similarly impractical
		// step sizes during ODE integration.
		virtual Number		tauAbsMin() { return 1*UOM::microsec; }

		// ODE Solver Interface functions -----------------------------

		// Load the alpha-beta table upon start
		virtual void		simulationStarted();

		// Delete the alpha-beta table upon end
		virtual void		simulationEnded();

		// Fastpath computation of derivatives
		virtual void		computeDerivatives();

		// Fastpath computation of local state update
		virtual void		localStateUpdate(SimTime h, CNStepType stepType);

		// Debug and reporting interface functions --------------------
		virtual void		printAlphaBetaTable(char* pathName=NULL);

	protected:

		// Flag indicating that load in progress (single thread assumed)
		static bool			_alphaBetaLoadInProg;

		// Keep track of which object is loading the alpha beta table
		// so that caches can be made unique to the load object.
		// This usage is not thread safe.
		static VoltageDepTabChannel*	_alphaBetaLoadObject;

		// Cache values in each instance to speed up processing
		AlphaBetaEntry**	_pABTable;		// cached pAlphaBetaTable

		// Subclass responsibilities ----------------------------------

		// pAlphaBetaTable must be supplied by any subclass using
		// table lookups for alpha-beta computations. Every subclass
		// will need to define its own alphaBetaTable array.
		virtual AlphaBetaEntry**	pAlphaBetaTable() = 0;
	};

	// ----------------------------------------------------------------
	// CLASS:	VoltageDepTabGate
	// EXTENDS:	VoltageDepTabChannel
	// DESC:	Abstract class for single gate variables
	// RESP:
	//		1.	Allow definition of gates by providing defaults
	//			for pure functions used in channels.
	//
	// NOTE:	The class hierarchy isn't an "is a" pattern.
	// ----------------------------------------------------------------

	class VoltageDepTabGate : public VoltageDepTabChannel {

	public:

		// Constructors and destructor
		VoltageDepTabGate() {}
		virtual ~VoltageDepTabGate() {}

		// Accessors -- suppress side-effects of setting container
		inline Compartment*		container() { return _container; }
		virtual void			container(Compartment* comp) { _container = comp; }

		// Dummy functions
		virtual Number		Vrev() { return 0; }
		virtual Number		conductance() { return 0; }

		// ODE Solver interface functions -----------------------------

		virtual bool		isIonChannel() { return false; }
		virtual bool		isGateVariable() { return true; }
	};

	// ----------------------------------------------------------------
	// CLASS:	EnergyBarrierTabChannel
	// EXTENDS:	VoltageDepTabChannel
	// DESC:	Implements an energy barrier style ion channel.
	//			This formulation uses a energy barrier model
	//			as described by Eyring rate theory. This provides a 
	//			consistent way to parameterize channel gates.
	//
	//			The form of the energy barrier equation used here is:
	//
	//			xinf = 1/(1+exp(-(v-vhalf)/k))
	//			tau = tmin+1/(rate*(exp(g*(v-vhalf)/k)+(exp(-(1-g)*(v-vhalf)/k))))
	//
	//			where 
	//				xinf = equilibrium value of state variable
	//				tau = time constant
	//				vhalf = voltage for half response value of xinf
	//				k = param controlling slope of response (1/(4k) = dx/dv)
	//				tmin = minimum value for tau
	//				rate = value set to ensure a specified maximum for tau
	//				g (gamma) reflects the barrier location as proportion
	//
	//			Note that k>0 for activation gates and k<0 for inactivation.
	//
	//			Implicit in these parameters is a F/RT term in the exponents.
	//			Temperature adjustments are handled automatically. However,
	//			an overall adjustment of time constants based on Q10 is also
	//			applied to the tau value derived from a given voltage.
	//
	// RESP:
	//		1.	Get/set channel parameters
	//		2.	Provide xinf and tau for table loading
	//		3.	Compute state variable derivatives
	//
	// NOTE:	See Borg-Graham LJ, 1998. Interpretations of
	//			Data and Mechanisms for Hippocampal Cell Models,
	//			in Cerebral Cortex vol 13. New York: Plenum Press.
	//			Also available via Surf-Hippo web site.
	//			Eyring rate theory is covered in a variety of standard
	//			references including Johnston & Wu and Hille.
	//
	//			The parameters of Eyring rate theory are sometimes
	//			represented differently from the formulas used here. 
	//			Eyring rate theory specifies alpha and beta in the form:
	//
	//			alpha = k0*exp(zeta*gamma*(v-vhalf)*F/(R*T))
	//			beta =  k0*exp(-zeta*(1-gamma)*(v-vhalf)*F/(R*T))
	//
	//			If this form of parameterization is used, subclasses
	//			should override functions rate() and zeta() specifying
	//			parameter values for k0 and zeta.
	//
	//			A subclass in which values are to be instantaneous
	//			should override numStateVar and value. Otherwise a
	//			small value of tau will result in small time step
	//			sizes being needed during ODE evaluation.
	// ----------------------------------------------------------------

	class EnergyBarrierTabChannel : public VoltageDepTabChannel {

	public:

		// Constructors and destructor
		EnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~EnergyBarrierTabChannel();

		// Accessors used by subclasses to provide parameters ---------
		
		// Temperature (deg C) at which these rate parameters were derived
		virtual Number	ratedTempC()=0; // e.g. { return 22; } for room temp.

		// Specify a separate Q10 value for the tauMin value
		// By default, the tauMin value is not temperature sensitive
		// except for the changes implied by Fv/RT terms.
		virtual Number	Q10ForTauMin() { return 1; }

		// Voltage which elicits half the steady state response
		virtual Number	Vhalf() = 0;

		// Response slope (1/(4*slope) = dx/dv at v=vhalf)
		// A value of 0 is impossible and means that no value is set.
		// Do not specify slope if other parameters are voltage dependent.
		virtual Number	slope() { return 0; }

		// Maximum value of tau at any voltage.
		// A value less than 0 means that no value supplied.
		// Do not specify tauMax if other parameters are voltage dependent.
		virtual Number	tauMax() { return -1; }	// -1 = no value set

		// Minimum value of tau. 
		virtual Number	tauMin() { return 0; }

		// Parameter specifying barrier voltage sensor location as proportion.
		virtual Number	gamma() { return 0.5; }

		// Rate derived from tau values given.
		// If Eyring rate formula is used, supply the rate value here.
		// rate=0 indicates that tauMin supplies time constant.
		virtual Number	rate();

		// Zeta value computed from slope and temperature.
		// If Eyring rate formula is used, supply the zeta value here.
		virtual Number	zeta();

		// Some models use voltage dependent parameter values.
		// The following functions should be overridden in that case.
		virtual Number	tauMinAtV(Number v)		{ return tauMin(); }
		virtual Number	gammaAtV(Number v)		{ return gamma(); }
		virtual Number	rateAtV(Number v)		{ return rate(); }
		virtual Number	zetaAtV(Number v)		{ return zeta(); }

		// Compute xinf-tau for table (converted to alpha-beta later)
		virtual	Number	xinfForTable(Number vm);
		virtual Number	tauForTable(Number vm);

		// Functions associated with temperature adjustments ----------

		// Get F/RT at the rated temperature and cache it
		virtual Number	FoverRTAtRatedTemp();

		// Compute the Q10 factor for tauMin and cache it
		virtual Number	Q10FactorForTauMin();
		
		// Other interface functions ----------------------------------
		virtual void	loadAlphaBetaTable();	// wrap superclass loader

	protected:

		// Cached instance values
		Number	_cachedQ10FactorForTauMin;	
	
		// Cached temp values used to speed up table loading
		// These are set during the load process for the current object. 
		static  Number	_cachedFoverRT;
		static  Number	_cachedRate;
		static  Number	_cachedZeta;
	};

	// ----------------------------------------------------------------
	// CLASS:	EnergyBarrierTabGate
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for single gate variables
	//			following the same logic as energy barrier
	//			channels.
	// RESP:
	//		1.	Allow definition of gates by providing defaults
	//			for pure functions used in channels.
	//
	// NOTE:	The class hierarchy isn't an "is a" pattern.
	//			See EnergyBarrierTabChannel for gate parameters.
	// ----------------------------------------------------------------

	class EnergyBarrierTabGate : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		EnergyBarrierTabGate() {}
		virtual ~EnergyBarrierTabGate() {}

		// Accessors -- suppress side-effects of setting container
		inline Compartment*		container() { return _container; }
		virtual void			container(Compartment* comp) { _container = comp; }

		// Dummy functions
		virtual Number		Vrev() { return 0; }
		virtual Number		conductance() { return 0; }

		// ODE Solver interface functions -----------------------------
		virtual bool			isIonChannel() { return false; }
		virtual bool			isGateVariable() { return true; }
	};

	// ----------------------------------------------------------------
	// CLASS:	Order1EnergyBarrierTabChannel
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for energy barrier channels
	//			having a single variable with exponent 1.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// NOTE:	See EnergyBarrierTabChannel for parameters.
	// ----------------------------------------------------------------

	class Order1EnergyBarrierTabChannel : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		Order1EnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~Order1EnergyBarrierTabChannel();

		// Get conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	Order2EnergyBarrierTabChannel
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for energy barrier channels
	//			having a single variable with exponent 2.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// NOTE:	See EnergyBarrierTabChannel for parameters.
	// ----------------------------------------------------------------

	class Order2EnergyBarrierTabChannel : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		Order2EnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~Order2EnergyBarrierTabChannel();

		// Get conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	Order3EnergyBarrierTabChannel
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for energy barrier channels
	//			having a single variable with exponent 3.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// NOTE:	See EnergyBarrierTabChannel for parameters.
	// ----------------------------------------------------------------

	class Order3EnergyBarrierTabChannel : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		Order3EnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~Order3EnergyBarrierTabChannel();

		// Get conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	Order4EnergyBarrierTabChannel
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for energy barrier channels
	//			having a single variable with exponent 4.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// NOTE:	See EnergyBarrierTabChannel for parameters.
	// ----------------------------------------------------------------

	class Order4EnergyBarrierTabChannel : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		Order4EnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~Order4EnergyBarrierTabChannel();

		// Get conductance and current for channel
		virtual Number			conductance();
		virtual Number			Iion();
		virtual void			condAndIion(Number& Gout, Number& Iout);
	};

	// ----------------------------------------------------------------
	// CLASS:	Order1CaEnergyBarrierTabChannel
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for energy barrier channels
	//			having a single variable with exponent 1.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// NOTE:	See EnergyBarrierTabChannel for parameters.
	// ----------------------------------------------------------------

	class Order1CaEnergyBarrierTabChannel : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		Order1CaEnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~Order1CaEnergyBarrierTabChannel();

		// Indicate that this channel is a source of calcium
		virtual bool isCalciumChannel() { return true; }

		// Get conductance for channel (using GHK conductance)
		virtual Number			conductance();

		// Get current for channel (using GHK potential)
		virtual Number			Iion();

		// Dummy function -- use Nernst potential
		virtual Number			Vrev() { return 140*UOM::mV; }
	};

	// ----------------------------------------------------------------
	// CLASS:	Order2CaEnergyBarrierTabChannel
	// EXTENDS:	EnergyBarrierTabChannel
	// DESC:	Abstract class for energy barrier channels
	//			having a single variable with exponent 2.
	// RESP:
	//		1.	Compute conductance.
	//		2.	Compute current.
	//
	// NOTE:	See EnergyBarrierTabChannel for parameters.
	// ----------------------------------------------------------------

	class Order2CaEnergyBarrierTabChannel : public EnergyBarrierTabChannel {

	public:

		// Constructors and destructor
		Order2CaEnergyBarrierTabChannel(Number gSpVal=0);
		virtual ~Order2CaEnergyBarrierTabChannel();

		// Indicate that this channel is a source of calcium
		virtual bool isCalciumChannel() { return true; }

		// Get conductance for channel
		virtual Number			conductance();

		// Get current for channel (performance optimization)
		virtual Number			Iion();

		// Dummy function -- use Nernst potential
		virtual Number			Vrev() { return 140*UOM::mV; }
	};

	// ----------------------------------------------------------------
	// CLASS:	ActionPotentialEvent
	// EXTENDS:	Event
	// DESC:	Event describing an action potential.
	// RESP:
	//		1.	Know event type
	//		2.	Know axon where action potential started
	//		3.	Know recipient synapse
	//		4.	Know quantity
	//		5.	Know estimated firing rate
	//
	// NOTE		The event time is set based on when the
	//			action potential will be received and acted on.
	//			The meaning of quantity is specific to the receptor,
	//			but generally corresponds to mean quanta released.
	//			If set stochastically, then once set it would not
	//			be changed unless different outcomes are possible.
	//			Firing rate and isi are provided by the event originator
	//			though the receiving synapse may use other estimates.
	// ----------------------------------------------------------------

	class ActionPotentialEvent : public Event {
	
	public:

		// Constructors and destructor
		ActionPotentialEvent(
			SimTime t = InfinitePast, 
			AxonProcess* ax = NULL, 
			Synapse* syn = NULL);
		virtual ~ActionPotentialEvent();

		// Accessors
		inline	AxonProcess*		axon() { return _axon; }
		inline  void				axon(AxonProcess* a) { _axon = a; }

		inline  Synapse*			synapse() { return _synapse; }
		inline  void				synapse(Synapse* s) { _synapse = s; }

		inline  Number				quantity() { return _quantity; }
		inline  void				quantity(Number q) { _quantity = q; }

		inline  Number				firingRate() { return _firingRate; }
		inline  void				firingRate(Number q) { _firingRate = q; }

		inline	SimTime				isi() {	return _isi; }
		inline  void				isi(SimTime t) { _isi = t; }

		inline	bool				isFinal() { return _isFinal; }
		inline	void				isFinal(bool flag) { _isFinal = flag; }

		// Identify the class of action potential events
		static unsigned int			eventClassId() { return 0x0001; }

	protected:
		AxonProcess*			_axon;				// axon where AP originated
		Synapse*				_synapse;			// target synapse
		Number					_quantity;			// mean quanta released
		Number					_firingRate;		// estimated rate of firing
		Number					_isi;				// inter-spike interval
		bool					_isFinal;			// one-time adjustments are done
	};

	// ----------------------------------------------------------------
	// CLASS:	ActionPotentialEventQueue
	// EXTENDS:	none
	// DESC:	Provide a specialized queue of time ordered
	//			action potential events
	// RESP:
	//		1.	Maintain queue of events in event time order
	//		2.	Allow inspection of events in order.
	//
	// NOTES:	This event queue is a circular queue optimized 
	//			for use by SynapticResponse objects. Such queues 
	//			are typically not very large and events are added
	//			in a sequence that is almost time ordered.
	//
	//			Heap objects (priority_queue) could perform similar 
	//			functions but do not support examining entries
	//			in sequence before removing them from the heap.
	//
	//			This class assumes that ActionPotentialEvent
	//			destructor is a no-op and does not invoke it
	//			before deleting storage holding such events.
	//
	//			Storing object instances on the queue rather pointers,
	//			as is done in other cases, avoids frequent allocation
	//			and freeing of heap storage.
	// ----------------------------------------------------------------

	class ActionPotentialEventQueue {
	
	public:

		// Constructors and destructor
		ActionPotentialEventQueue(int initialCapacity=0);
		virtual ~ActionPotentialEventQueue();

		// Accessors
		inline  int				size() { return _size; }
		inline  int				capacity() { return _capacity; }

		// Begin, end, next, and previous only return non-NULL pointers
		// to valid entries. This resolves the ambiguity that would otherwise
		// arise when _begin==_end and also allows STD-like iteration loops.
		inline  ActionPotentialEvent*	begin() { return _size>0 ? _begin : NULL; }
		inline  ActionPotentialEvent*	end() { return NULL; }

		// Answer true if the queue is empty
		virtual bool			empty() { return _size==0; }

		// Set the queue to an empty state
		virtual void			clear();

		// Add an action potential event to the queue.
		// A potential side-effect of adding is reallocation of the
		// queue, rendering any extant pointers to its events invalid.
		virtual void			add(ActionPotentialEvent* apEvent);

		// Remove the first (oldest) entry from the queue. The queue cannot
		// be reallocated by removing an entry and pointers to events
		// remain valid as long as they do not point to removed events.
		virtual void			removeFirst();

		// Return a pointer to the next entry beyond the pointer provided.
		// If currentPointer is NULL, return a pointer to the first entry.
		// Return NULL if the end of the queue is reached.
		virtual ActionPotentialEvent* next(ActionPotentialEvent* currentPointer);

		// Return a pointer to the previous entry before the pointer provided.
		// If currentPointer is NULL, return a pointer to the actual last entry.
		// If there are no previous entries, return NULL.
		virtual ActionPotentialEvent* previous(ActionPotentialEvent* currentPointer);

		// Access an event based on an offset from the first entry
		virtual ActionPotentialEvent*	peek(int n);
		
	protected:
		int						_size;		// Number of entries in the queue
		int						_capacity;	// Number of entries allocated less one
		ActionPotentialEvent*	_queue;		// Circular queue of events (first entry)
		ActionPotentialEvent*	_begin;		// First event in the queue
		ActionPotentialEvent*	_end;		// One past last event in the queue

		// Internal routines for incrementing and decrementing pointers
		inline  ActionPotentialEvent* privateNext(ActionPotentialEvent* ptr)
			{ return ptr==&_queue[_capacity-1] ? _queue : ptr+1; }

		inline  ActionPotentialEvent* privatePrev(ActionPotentialEvent* ptr)
			{ return ptr==_queue ? &_queue[_capacity-1] : ptr-1; }
	};

	// ----------------------------------------------------------------
	// CLASS:	SynapticResponse
	// EXTENDS:	IonChannel
	// DESC:	Abstract class representing a pooled
	//			response of one or more synaptic elements
	//			to action potential events.
	// RESP:
	//		1.	Define common protocol for action
	//			potential handling.
	//		2.	Provide a global counter for late arriving
	//			events, that is events that arrive after
	//			the evaluation interval in which they would
	//			have first been included.
	//		3.	Provide synaptic delay time for synapses.
	//			For a discussion of synaptic delay see:
	//			Sabatini BL and Regher WG. 1999.
	//			Timing of synaptic transmission.
	//			Annu. Rev. Physiol. 61:521-542.
	//		4.	Maintain the list of associated synapses.
	//
	// NOTE:	Typically this will be  used to represent
	//			a collection of synaptic receptors of a single
	//			type that share a compartment in common.
	//			Collectively processing the synaptic response
	//			can be much faster than handling each synapse
	//			as an individual entity.
	//
	//			When multiple types of responses are triggered
	//			from a single action potential event (e.g.
	//			AMPA and NMDA receptors in a single synapse),
	//			multiple responses share the list of synapses
	//			owned by a SynapticResponseGroup object.
	//
	//			Processing for the end of the time step is routed
	//			via applyEndOfTimeStep so that, when participating as
	//			a member of a group, a response object can be
	//			assured that purging the AP event queue is done
	//			only after all group members have finished their
	//			end of step processing.
	//
	//			Storage for synapses is allocated with a Synapse
	//			instance immediately followed by an appropriate
	//			subclass of SynapseState. When multiple responses
	//			are triggered, each has its own state instance.
	//			See SynapticResponseGroup for this case.
	// ----------------------------------------------------------------

	class SynapticResponse : public IonChannel {

	public:

		// Constructors and destructor
		SynapticResponse();
		virtual ~SynapticResponse();

		// Respond to basic queries as to type and status
		virtual bool			isSynapticResponse() { return true; }
		virtual bool			isSynapticGroup() { return false; }

		inline  bool			isGroupMember() { return groupOwner() != NULL; }
		inline  bool			isActive() {return synapses()!=NULL; }
		inline  bool			isInactive() {return synapses()==NULL; }

		// Access the owner if this is part of a group (otherwise NULL)
		inline	SynapticGroupResponse* groupOwner() { return _groupOwner; }

		// Accessors to items potentially held by the group owner
		inline  Synapse*		synapses()
		{ return _groupOwner==NULL ? _synapses : groupSynapses(); }

		inline  int				synapseCount()
		{ return _groupOwner==NULL ? _synapseCount : groupSynapseCount(); } 

		inline  ActionPotentialEventQueue& apQueue()
		{ return _groupOwner==NULL ? _apQueue : groupAPQueue(); }

		// Access the offset to state data in a synapse object.
		inline	unsigned int	synapseStateOffset() { return _synapseStateOffset; }
		virtual void			synapseStateOffset(unsigned int nbytes) 
								{ _synapseStateOffset = nbytes; }

		// Response delay for associated synapses (value is default only)
		virtual SimTime			synapticDelay() { return 250*UOM::microsec; }

		// Create a synapse to go with this response. The subclass is asked to
		// create any associated state data to go with the synapse.
		virtual Synapse*		createSynapse(
			AxonProcess*			axon = NULL,				// presynaptic axon process
			Number					wght = 1.0,					// initial weight value
			Number					dist = 100*UOM::micron);	// axonal distance

		// Destroy a synapse and any associated state data.
		// Actual release of storage does not occur until both
		// pre and postsynaptic access has been removed.
		virtual void			destroySynapse(Synapse* syn);

		// Accessors for the synaptic weight within synaptic state data.
		// This assumes that state data is a subclass of SynapticState.
		virtual Number			synapseWeight(Synapse* syn);			// get
		virtual void			synapseWeight(Synapse* syn, Number w);	// set

		// Return the total of all synapses weights for this response type.
		// If requested, also return a count of synapses.
		virtual Number			totalSynapticWeight(int* pCount=NULL);

		// As a convenience for testing, provide an interface for
		// disabling synaptic plasticity
		virtual void			disablePlasticity() {}

		// Other framework interfaces ---------------------------------

		// Set the owner, but only if there are no existing synapses.
		// This is normally only done by framework classes.
		virtual void			groupOwner(SynapticGroupResponse* owner);

		// Handle a new action potential event
		virtual void			signalActionPotential(
								SimTime t,						// Spike time
								ActionPotentialEvent* apEvent);	// Template event

		// At the start of the simulation allow class caches to be allcated.
		virtual void			simulationStarted();

		// For end of time step, notify subclasses and then purge the AP queue.
		// To ensure synchronous operation, this is only done for the queue owner
		// who, in turn, routes to all member responses via applyEndOfTimeStep.
		virtual void			timeStepEnded();
		virtual bool			notifyOnTimeStepEnded() { return !isGroupMember(); } 

		// Provide access to members for reporting. If this is not
		// a group item, return the empty member vector.
		virtual SynapticResponseVector& members() { return _members; } 

		// Subclass responsibilities ----------------------------------

		// Do subclass specific initialization of any class cache.
		virtual void			initializeClassCache() {}

		// Provide the size to allocate for state data.
		// By default, sufficient storage is allowed for a SynapseState instance.
		virtual unsigned int	synapseStateSize();

		// Apply a constructor to synapse state data. By default an instance of
		// synapse state is constructed.
		virtual void			createSynapseState(Synapse* syn, Number wght);

		// Destroy any synapse state data.
		virtual void			destroySynapseState(Synapse* syn) {}

		// Handle changes in synapse list including via a group owner.
		// This allows a subclass to keep its own extra records of synapses.
		virtual void			synapseAdded(Synapse* syn) {}
		virtual void			synapseRemoved(Synapse* syn) {}

		// Handle arrival of an action potential in which the event time
		// precedes the current evaluation time of this response.
		// Default action is to reschedule the event for the current time.
		virtual void			handleLateAPEvent(ActionPotentialEvent* apEvent)
								{ apEvent->eventTime( currentTime() ); }

		// Do any necessary updates before the AP queue is purged.
		virtual void			applyEndOfTimeStep() {}

		// Provide a default component name (should be overridden)
		virtual const char*		componentName() { return "SynapticResponse"; }

		// Static Accessors -------------------------------------------

		static int				lateArrivingEventCount();
		static SimTime			meanLateEventTime();

		// Set late arriving event stats to zero.
		static void				resetLateArrivingEventStatistics();

	protected:
		Synapse*				_synapses;				// list header for synapses
		SynapticGroupResponse*	_groupOwner;			// common owner if part of a group
		int						_synapseCount;			// count of associated synapses
		unsigned int			_synapseStateOffset;	// offset to state data in synapse buffer
		ActionPotentialEventQueue _apQueue;				// time ordered event queue
		SynapticResponseVector	_members;				// group members if any

		// Global counters for late arriving events
		static int				_LateArrivingEventCount;
		static SimTime			_TotalLateEventTime;

		// Access via group owner. This is needed because C++ cannot
		// resolve virtual function calls within inline accessors above.
		virtual Synapse*		groupSynapses();
		virtual int				groupSynapseCount();
		virtual ActionPotentialEventQueue& groupAPQueue();

		// Maintain synapse relationships.
		virtual void			addSynapse(Synapse* syn);
		virtual void			removeSynapse(Synapse* syn);

		// Destroy all current synapses associated with this object
		virtual void			destroyAllSynapses();

		// Subclass responsibilities ----------------------------------

		// Enqueue an AP event after making any needed adjustments
		virtual void			addAPEventToQueue(ActionPotentialEvent* apEvent);
	};

	// ----------------------------------------------------------------
	// CLASS:	SynapticGroupResponse
	// EXTENDS:	SynapticResponse
	// DESC:	Represents a collection of types of
	//			ion channel receptors that occur
	//			together in synapses.
	// RESP:
	//		1.	Maintain a collection of receptor types
	//		2.	Maintain synapse list
	//		3.	Pass synapse updates to members
	//		4.	Respond to Iion etc by summing member values
	//		5.	Notify member of action potential events
	//		6.	Notify members when time step ends
	// ----------------------------------------------------------------
	
	class SynapticGroupResponse : public SynapticResponse {

	public:

		// Constructors and destructor
		SynapticGroupResponse();
		virtual ~SynapticGroupResponse();

		// Polymorphic add/remove for convenience
		inline  void			add(SynapticResponse* resp)		{ addResponse(resp); }
		inline  void			remove(SynapticResponse* resp)	{ removeResponse(resp); }

		// Maintain a collection of group members
		virtual void			addResponse(SynapticResponse* resp);
		virtual void			removeResponse(SynapticResponse* resp);

		// When setting model, pass it along to members
		inline  Model*			model() { return _model; }
		virtual void			model(Model* m);

		// When setting container, pass it along to members
		inline  Compartment*	container() { return SynapticResponse::container(); }
		virtual void			container(Compartment* comp);

		// Provide the size to allocate for state data.
		// This is the sum of the sizes reqested by all group members.
		virtual unsigned int	synapseStateSize();	

		// Access flag that controls whether initial weight assignments
		// are applied to all group members or only the first member.
		virtual bool			applyInitialWeightToAll() { return false; }

		// Apply a constructor to synapse state data.
		// Each member of the group creates its own state in turn.
		// If applyInitialWeightToAll is false, then weight
		// is only provided to the first member response.
		// Otherwise, all members receive the same weight.
		virtual void			createSynapseState(Synapse* syn, Number wght);

		// Apply a destructor to synapse state data.
		// Each member of the group destroys its own state in turn.
		virtual void			destroySynapseState(Synapse* syn);

		// Handle changes in synapse list by notifying members.
		virtual void			synapseAdded(Synapse* syn);
		virtual void			synapseRemoved(Synapse* syn);

		// Inform members of a late arriving action potential event
		virtual void			handleLateAPEvent(ActionPotentialEvent* apEvent);

		// Notify members that the end of the time step was reached
		virtual void			applyEndOfTimeStep();

		// Response delay for associated synapses.
		// By default, value is taken from the first response member.
		virtual SimTime			synapticDelay();

		// Provide an aggregate response from members.
		// When there are no synapses, these return 0 immediately.
		virtual Number			conductance();
		virtual Number			Iion();
		virtual Number			ICa();
		virtual void			condAndIion(
								Number&	Gout,			// conductance
								Number& Iout);			// current

		// Functions required by the ion channel and ODE frameworks
		virtual Number			Vrev() { return 0; }
		virtual int				numStateVar() { return 0; }

		// Indicate that this is a group response
		virtual bool			isSynapticGroup() { return true; }

		// Add to components to probe when reporting.
		// This includes the current object and any members.
		virtual void			addToComponentsToProbe(ModelComponentVector& comps);

		// Locate the identified component within this ion channel.
		// If the component is not found, return NULL.
		virtual IonChannel*		findIonChannel(TokenId compId);

		// Pass any neuromodulation parameters to members
		virtual void			setModParams(
			TokenId				modId,			// modulation type token id
			int					numParams,		// number of parameter values
			Number*				params);		// parameter values as an array

		// As a convenience for testing, provide an interface for
		// disabling synaptic plasticity
		virtual void			disablePlasticity();
	};

	// ----------------------------------------------------------------
	// CLASS:	SynapticConductance
	// EXTENDS:	SynapticResponse
	// DESC:	Represents a collection of individual
	//			synaptic ion channels of a common type.
	//			Conductance is computed for the collection
	//			of synapses since they follow a common
	//			time course in their responses.
	// RESP:
	//		1.	Know conductance for a receptor
	//		2.	Create and locate receptor specific synapse state
	//		3.	Apply synaptic plasticity rules
	//		4.	Get Iion etc as a function of time (in subclasses)
	//		5.	Free storage for state data along with synapse
	//		6.	Adjust weight to compensate for spine neck
	//		7.	Provide subclass utilities for scaling conductance
	//		8.	Provide notifications to plasticity rules
	//
	// NOTES:	In this hierarchy, the conductance multiplier
	//			variable (_g) is intended as the value for 
	//			an individual synapse and does not represent
	//			a value for all synapses together.
	//
	//			Synapse state data is made up of plasticity
	//			state data and other state data. Plasticity
	//			state data is used by the plasticity rule
	//			and is typically a subclass of PlasticityState.
	//			Receptor specific state data is typically a 
	//			subclass of SynapseState and occurs before the 
	//			associated plasticity state data.
	// ----------------------------------------------------------------
	
	class SynapticConductance : public SynapticResponse {

	public:

		// Constructors and destructor
		SynapticConductance();
		virtual ~SynapticConductance();

		// Suppress setting of conductance except via gMax - Fatal error if used.
		virtual void			gSpecific(Number gval);
		virtual void			gAbsolute(Number gval);

		// Accessor for presynaptic plasticity rules
		inline  PlasticityRule* presynapticRule() { return _presynapticRule; }
		virtual void			presynapticRule(PlasticityRule* pr);

		// Accessor for postsynaptic synaptic plasticity rules
		inline  PlasticityRule* postsynapticRule() { return _postsynapticRule; }
		virtual void			postsynapticRule(PlasticityRule* pr);

		// Unhook from either a presynaptic or postsynaptic rule.
		// The invoker is responsible for deleting pr ultimately.
		virtual void			unhookFromRule(PlasticityRule* pr);

		// As a convenience for testing, provide an interface for
		// disabling pre and postsynaptic plasticity. 
		// Current rules objects are deleted.
		virtual void			disablePlasticity();

		// When the time step is over, update state
		virtual void			applyEndOfTimeStep();

		// Handle changes in synapse list including via a group owner.
		virtual void			synapseRemoved(Synapse* syn);

		// When setting model, pass along to any plasticity dependents
		inline  Model*			model() { return _model; }
		virtual void			model(Model* m);

		// Interfaces for managing synaptic state data ----------------

		// Return the size to allocate for state data including any needed
		// to support a plasticity rule.
		virtual unsigned int	synapseStateSize();

		// Apply a constructor to synapse state data.
		// This creates receptor specific state data plus any plasticity state.
		virtual void			createSynapseState(Synapse* syn, Number wght);

		// Apply a destructor to plasticity and receptor specific state data.
		virtual void			destroySynapseState(Synapse* syn);

		// Access the offset to state data in a synapse object.
		// Setting a new value updates plasticity state offsets also.
		inline	unsigned int	synapseStateOffset() { return _synapseStateOffset; }
		virtual void			synapseStateOffset(unsigned int n);

		// Provide offsets for synapse and plasticity state data
		inline  unsigned int	receptorSpecificStateOffset() { return _synapseStateOffset; }

		// Subclass responsibilities ----------------------------------

		// Access the peak conductance achieved in a single response
		virtual Number			gMax() = 0;
		virtual void			gMax(Number gm) = 0;

		// Provide the spine neck resistance to use in adjusting response.
		virtual Number			spineNeckResistance() { return 0; }	// default only

		// Provide a release probability (e.g. for use in presynaptic rules).
		virtual Number			releaseProbability() { return 1; }	// default only

		// Apply any neuromodulation parameters. Default is to defer to
		// superclass and also pass along to plasticity rules.
		virtual void			setModParams(
			TokenId				modId,			// modulation type token id
			int					numParams,		// number of parameter values
			Number*				params);		// parameter values as an array

		// Return the size to allocate for receptor specific state data.
		// By default this is the size of SynapseState.
		virtual unsigned int	receptorSpecificStateSize();

		// Create a receptor specific state object in a synapse buffer.
		// By default this creates an instance of SynapseState.
		virtual void			createReceptorSpecificState(Synapse* syn, Number wght);

		// Destory a receptor specific state object in a synapse buffer.
		// By default, this is a no-op.
		virtual void			destroyReceptorSpecificState(Synapse* syn) {}

	protected:
		PlasticityRule*			_postsynapticRule;			// rule for overall plasticity
		PlasticityRule*			_presynapticRule;			// rule for presynaptic plasticity
				
		// Return a weight adjusted to compensate for spine neck resistance.
		// For purposes of the adjustment, the response is assumed to be
		// proportional to synapse weight and AP quantity.
		virtual Number			adjustedSynapticWeight(ActionPotentialEvent* apEvent);

		// Subclass responsibilities ----------------------------------

		// Update state when the time step ends.
		virtual void			updateStateForTimeStep() = 0;

		// Update the current state to reflect queued AP events.
		// Default is to process each event in order.
		// Interval for events to included is [from to).
		// Subclass will need to invoke as appropriate in its processing loop.
		virtual void			updateForQueuedAPEvents(SimTime from, SimTime to);

		// Update the receptor state to reflect a single AP Event.
		// Subclass must override if this is used, e.g. via updateForQueuedAPEvents.
		// Subclass is responsible for ensuring that isPending is cleared in
		// each event that was processed at least once in a time step.
		virtual void			updateForAPEvent(ActionPotentialEvent* apEvent);

		// Clear state values to zero (default only)
		// Note that this can be invoked during destructor processing
		// for this class and can thus cannot be a pure virtual function. 
		virtual void			clearState() {}

		// Empty any caches (default only)
		virtual void			emptyCaches() {}

		// Utility functions for subclasses ---------------------------

		// Return the peak value for a dual exponent response
		virtual Number			peakDualExpResp(SimTime tau1, SimTime tau2);

		// Return the peak value for a triple exponent formulation:
		// where tau1<tau2, tau1<tau3, and 0<=c<=1.
		// Error tolerance in peak time is rtol*tau1.
		virtual Number			peakTripleExpResp(
			SimTime				tau1, 
			SimTime				tau2, 
			SimTime				tau3,
			double				c,
			double				rtol = 1e-2);
	};

	// ----------------------------------------------------------------
	// CLASS:	SingleExpSynapticCondClassCache 
	// EXTENDS:	ModelEntityClassCache
	// DESC:	Provides a cache of time step dependent values			
	// ----------------------------------------------------------------

	class SingleExpSynapticCondClassCache : public ModelEntityClassCache {

	public:
		// Constructors and destructor
		SingleExpSynapticCondClassCache() {}
		~SingleExpSynapticCondClassCache() {}

		double					Tau1;			// time constant at current temp
		double					Exp1;			// exp(-h/tau1)
	};

	// ----------------------------------------------------------------
	// CLASS:	SingleExpSynapticCond
	// EXTENDS:	SynapticConductance
	// DESC:	Represents a collection of synapses
	//			with an immediate rise time and a single
	//			exponent decay form of synaptic response. 
	//			This can also be extended for more complex
	//			models involving multiple time constants.
	// RESP:
	//		1.	Explicitly solve ODE for conductance
	//		2.	Integrate effects of action potentials
	//		3.	Maintain a cache of values common to instances.
	//
	// NOTES:	Conductance is found by solving
	//
	//			ds1/dt = -s1/tau1 + x(t)
	//
	//			where x(t) is the input (as a weighted impulse).
	//			Each individual synapse has its own weight 
	//			value which is reflected in the value of x(t). 
	//
	//			Resulting conductance is g*s1(t). Since this object
	//			can be used as part of multiple exponent models,
	//			the state variables is termed s1 here where further
	//			state variables might be s2, s3, etc. Note that
	//			s1 is a state variable internal to the instances and
	//			does not appear in the model state vector.
	//
	//			The conductance multiplier of an individual
	//			synapse is carried in the variable _g, though
	//			the synapse object is not obligated to use it
	//			and an adaptive (plastic) synapse could have its
	//			own idea of conductance.
	//
	//			Because inputs are impulses, s1 can be discontiguous.
	//			To avoid ODE solver complications, s1 is solved for
	//			explicitly as s1(t+h)=s1(t)*exp(-h/tau1). Subclasses
	//			may need to use a mean value of s1 over an interval
	//			to achieve higher precision than simply sampling
	//			the most recent value.
	//
	//			A cache of values common to all instances is maintained
	//			in the class cache. This avoids redundant computation
	//			especially when time step sizes are changed. Similarly
	//			the value of tau1 is initially cached after being 
	//			adjusted with the Q10 temperature factor.
	//
	//			Ohm's law is used to compute currents. This permits
	//			faster computation of conductance and current together,
	//			but subclasses may override if Ohm's law does not apply.
	// ----------------------------------------------------------------
	
	class SingleExpSynapticCond : public SynapticConductance {

	public:

		// Constructors and destructor
		SingleExpSynapticCond(Number gMaxValue=0);
		virtual ~SingleExpSynapticCond();

		// Access peak synapse conductance
		inline  Number			gMax() { return _gMax; }
		virtual void			gMax(Number maxCond);

		// Access the s1 value as of the current time
		inline Number			s1() { updateCachedValues(); return _s1Now; }

		// Estimate the mean value of s1 over time interval [t-h t)
		// where t is the current time.
		virtual Number			s1Mean(SimTime h);

		// Model interface --------------------------------------------

		// Define state variables -- none since solution is explicit
		virtual int				numStateVar() { return 0; }

		// Set initial state values to zero
		virtual void			setInitialState();

		// Do nothing for computation of derivatives (if requested)
		virtual void			computeDerivatives() {}

		// Return the total conductance for synapses of this type
		virtual Number			conductance() { return g()*s1(); }

		// Subclass responsibilities ----------------------------------

		// Time constant at channel's rated temperature (before Q10 adjustment).
		// When used as a single time constant, this is a decay rate.
		// When used via dual exponent subclass, this will be the rising rate.
		virtual SimTime			tau1() = 0;

		// Get both current and conductance at once using Ohm's law.
		// Subclass may override as appropriate for alternatives.
		virtual void			condAndIion(
								Number&	Gout,			// conductance
								Number& Iout);			// current

		// Reporting interface ----------------------------------------

		virtual int				numInternalStateVar() { return 1; }
		virtual Number			internalStateValue(int k);

		virtual const char*		componentName() { 
			return "SingleExpSynapticCond"; } // subclass resp

		virtual const char**	stateLabels() {
			static const char* slab[] = {"s1"}; 
			return slab; }

	protected:
		Number					_gMax;			// maximum conductance per synapse

		// Internal state variables
		SimTime					_cachedTime;	// time for state values 
		Number					_s1Start;		// value of s1 at start of time step
		Number					_s1Now;			// cached value of s1 state

		// Update the cached values as of the current time
		virtual void			updateCachedValues();

		// Update state when the time step ends
		virtual void			updateStateForTimeStep();

		// Handle an action potential arriving after the time step
		// in which it would be processed is already completed.
		virtual void			handleLateAPEvent(ActionPotentialEvent* apEvent);

		// Empty cache for time step
		virtual void			emptyCaches();

		// Subclass responsibilities  ---------------------------------
		
		// Address class cache data for a subclass.
		// Typically the subclass will allocate this data as a static.
		// There is no attempt to make accessing the cache thread safe.
		virtual SingleExpSynapticCondClassCache* 
								pSingleExpClassCache()=0;

		// Initialize the class cache (subclass may extend this)
		virtual void			initializeClassCache();

		// Advance state in time by the amount of the time step size
		virtual void			advanceState();

		// Save the current state as the starting state for a time step
		virtual void			saveStartingState();

		// Restore the current internal state to time step start values
		virtual void			restoreToStartingState();

		// Clear the state to initial zero values
		virtual void			clearState();

		// Update cached values to reflect change in time step
		virtual void			updateCacheForStepSize();

		// Update current conductance to reflect an AP event
		virtual void			updateForAPEvent(ActionPotentialEvent* apEvent);
	};

	// ----------------------------------------------------------------
	// CLASS:	DualExpSynapticCondWithVarTau
	// EXTENDS:	SingleExpSynapticCond
	// DESC:	Represents a collection of synapses
	//			with a dual exponent form of synaptic
	//			response in which the falling time constant (tau2)
	//			is a variable dependent on time or other states.
	// RESP:
	//		1.	Use the ODE solver to get conductance.
	//		2.	Set maximum conductance on estimated peak value.
	//
	// NOTE:	Conductance is found by solving
	//
	//			ds1/dt = -s1/tau1 + x(t)
	//			ds2/dt = -s2/tau2 + s1/tau1
	//
	//			where x(t) is the input (as a weighted impulse).
	//
	//			Resulting conductance is g*s2(t).
	//
	//			The value of tau2 is assumed to be a continuous variable
	//			and may depend on other state variables such as membrane
	//			potential or calcium concentration. Q10 adjustments are
	//			applied to both tau1 and tau2.
	//
	//			Note that using s1/tau1 rather than s1 in ds2/dt results
	//			in s1 and s2 both being dimensionless. This form is
	//			useful for scaling s1 and s2 to a consistent range and
	//			gives the same result as the more typical formulation 
	//			in which s1 has units of 1/time.
	//
	//			When tau2 is constant, DualExpSynapticCond performs the
	//			same function with less processing and higher accuracy.
	//			Whenever possible, it should be used instead of this class.
	// ----------------------------------------------------------------

	class DualExpSynapticCondWithVarTau : public SingleExpSynapticCond {

	public:

		// Constructors and destructor
		DualExpSynapticCondWithVarTau(Number gMaxValue=0);
		virtual ~DualExpSynapticCondWithVarTau();

		// Accessors
		inline  Number			s2() { return stateValue(0); }

		// Set synapse conductance multiplier based on maximum value
		// as determined by tau1 and tau2 at simulation start.
		virtual void			gMax(Number maxCond);

		// Return the total conductance for synapses of this type
		virtual Number			conductance() { return g()*s2(); }

		// Model interface --------------------------------------------

		// One state variable is used (other is local to SingleExpSynapticCond)
		virtual int				numStateVar() { return 1; }

		// Compute derivatives by directly applying the ODE
		virtual void			computeDerivatives();

		// Update the local state using an implicit trapezoidal rule
		virtual bool			canPerformLocalUpdate() { return true; }
		virtual void			localStateUpdate(SimTime h, CNStepType stepType);

		// Reporting interface ----------------------------------------

		virtual const char*		componentName() {
			return "DualExpSynapticCondWithVarTau"; } // subclass resp.

		virtual const char**			stateLabels() {
			static const char* slab[] = {"s1","s2"}; 
			return slab; }

		// Subclass responsibilities ----------------------------------

		// Parameter accessors

		// Time constant at channel's rated temperature (before Q10 adjustment).
		// Tau1 is provided via superclass and is rising rate. Tau2 is the decay rate.
		virtual SimTime			tau2() = 0;

		// Return a tau2 value to be used to set the gMax multiplier.
		// This is only accessed when setting the gMax scale factors
		// and is otherwise considered to be constant. By default
		// the initial value of tau2 is used.
		virtual SimTime			nominalTau2() { return tau2(); };

	protected:

		// Clear the state including the state vector value
		virtual void			clearState();
	};

	// ----------------------------------------------------------------
	// CLASS:	DualExpSynapticCondClassCache
	// EXTENDS:	SingleExpSynapticCondClassCache
	// DESC:	Defines public variables for the class cache.
	//			This cache allows individual dual exponent
	//			synapses to save computed values across the
	//			whole class rather than object-by-object.
	// ----------------------------------------------------------------

	class DualExpSynapticCondClassCache 
		: public SingleExpSynapticCondClassCache {

	public:
		double					Tau2;			// cached time constant at current temp
		double					Exp2;			// exp(-h/tau2)
		double					ImpulseResp;	// response impulse after one time step
		bool					useAlpha;		// tau1==tau2 ==> use alpha function

	};

	// ----------------------------------------------------------------
	// CLASS:	DualExpSynapticCond
	// EXTENDS:	SingleExpSynapticCond
	// DESC:	Represents a collection of synapses
	//			with a dual exponent form of synaptic
	//			response.
	// RESP:
	//		1.	Explicitly solve ODE for conductance
	//		2.	Integrate effects of action potentials
	//
	// NOTE:	Conductance is found by solving
	//
	//			ds1/dt = -s1/tau1 + x(t)
	//			ds2/dt = -s2/tau2 + s1/tau1
	//
	//			where x(t) is the input (as a weighted impulse).
	//
	//			Resulting conductance is g*s2(t).
	//
	//			tau1 and tau2 values are assumed to be constant
	//			throughout the simulation. The values returned
	//			during setInitialState processing are cached
	//			for efficiency and used throughout. Q10 factor
	//			adjustments are applied to both tau1 and tau2.
	//
	//			Note that using s1/tau1 rather than s1 in ds2/dt results
	//			in s1 and s2 both being dimensionless. This form is
	//			useful for scaling s1 and s2 to a consistent range and
	//			gives the same result as the more typical formulation 
	//			in which s1 has units of 1/time.
	// ----------------------------------------------------------------

	class DualExpSynapticCond : public SingleExpSynapticCond {

	public:

		// Constructors and destructor
		DualExpSynapticCond(Number gMaxValue=0);
		virtual ~DualExpSynapticCond();

		// Accessors
		inline  Number			s2() { updateCachedValues(); return _s2Now; }

		// Set synapse conductance multiplier based on maximum value
		virtual void			gMax(Number maxCond);

		// Return the group conductance for all synapses
		virtual Number			conductance() { return g()*s2(); }

		// Reporting interface ----------------------------------------

		virtual int				numInternalStateVar() { return 2; }
		virtual Number			internalStateValue(int k);

		virtual const char*		componentName() {
			return "DualExpSynapticCond"; }	// subclass resp.

		virtual const char**			stateLabels() {
			static const char* sl[] = {"s1", "s2"}; 
			return sl; }

		// Parameter accessors (subclass responsibility) --------------

		// Time constant at channel's rated temperature.
		// tau1 is provided via superclass and is rising rate. 
		// tau2 is the decay rate. If tau1==tau2 then an alpha 
		// function response is implied.
		virtual SimTime			tau2() = 0;		// falling time constant

	protected:
		Number					_s2Start;		// starting state value
		Number					_s2Now;			// current state value

		// Advance the current state values for the time step size in _H
		virtual void			advanceState();

		// Clear the state to initial zero values
		virtual void			clearState();

		// Save the current state as the starting state for a time step
		virtual void			saveStartingState();

		// Restore the current internal state to time step start values
		virtual void			restoreToStartingState();

		// Update cached values to reflect change in time step
		virtual void			updateCacheForStepSize();

		// Update current conductance to reflect an AP event
		virtual void			updateForAPEvent(ActionPotentialEvent* apEvent);

		// Address class cache data as used by a superclass
		virtual SingleExpSynapticCondClassCache* pSingleExpClassCache()
								{ return pDualExpClassCache(); }

		// Subclass responsibilities ----------------------------------

		// Address class cache data allocated by a subclass
		virtual DualExpSynapticCondClassCache* pDualExpClassCache()=0;

		// Initialize the class cache (subclass may extend this)
		virtual void			initializeClassCache();
	};

	// ----------------------------------------------------------------
	// CLASS:	DualExpSynapticCurrent
	// EXTENDS:	DualExpSynapticCond
	// DESC:	Represents a collection of synapses
	//			with a dual exponent form of synaptic
	//			response as a pure current source.
	// RESP:
	//		1.	Return ionic current
	//		2.	Return zero conductance
	//
	// NOTES:	Conductance from superclass is translated to the form of a 
	//			conductance using a nominal driving potential of one volt.
	//			Because this is a pure current source, conductance, i.e.,
	//			g in g*(V-Er), is always zero even though current is not.
	//
	//			While not directly a model of a physiological channel, 
	//			current sources are sometimes used in theoretical models
	//			and can represent distant synapses via an implied kernel.
	// ----------------------------------------------------------------

	class DualExpSynapticCurrent : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		DualExpSynapticCurrent(Number Imax=0) : DualExpSynapticCond(Imax) {}
		virtual ~DualExpSynapticCurrent() {}

		// Return current as superclass conductance times one volt
		// UOM conversions cancel out with those of constructor and can be omitted.
		virtual Number			Iion() {return DualExpSynapticCond::conductance(); }

		// Return zero if asked for this channel's conductance
		virtual Number			conductance() {return 0; }

		// Handle a request for conductance and current together
		virtual void			condAndIion(Number& Gout, Number& Iout)
		{ Gout = 0; Iout = DualExpSynapticCond::conductance(); }

		// Supply dummy value for Vrev (required but unused)
		virtual Number			Vrev() {return 0; }
	};

	// ----------------------------------------------------------------
	// CLASS:	TripleExpSynapticCondClassCache
	// EXTENDS:	SingleExpSynapticCondClassCache
	// DESC:	Defines public variables for the class cache.
	// ----------------------------------------------------------------

	class TripleExpSynapticCondClassCache 
		: public SingleExpSynapticCondClassCache {

	public:
		double					C2;				// cached weight value
		double					Tau2;			// cached time constant
		double					Tau3;			// cached time constant
		double					Exp2;			// exp(-h/tau2)
		double					Exp3;			// exp(-h/tau2)
		double					ImpulseResp2;	// response impulse after one time step
		double					ImpulseResp3;	// response impulse after one time step
		bool					useAlpha2;		// tau1 == tau2
		bool					useAlpha3;		// tau1 == tau3
	};

	// ----------------------------------------------------------------
	// CLASS:	TripleExpSynapticCond
	// EXTENDS:	SingleExpSynapticCond
	// DESC:	Represents a collection of synapses
	//			in which there is one rising time
	//			constant and two falling time constants.
	// RESP:
	//		1.	Explicitly solve ODE for conductance
	//		2.	Integrate effects of action potentials
	//
	// NOTE:	Conductance is found by solving
	//
	//			ds1/dt = -s1/tau1 + x(t)
	//			ds2/dt = -s2/tau2 + s1/tau1
	//			ds3/dt = -s3/tau3 + s1/tau1
	//
	//			where x(t) is the input (as a weighted impulse).
	//
	//			Resulting conductance is g*(c2*s2(t)+c3*s3(t))
	//			where c2+c3=1 and c2>0, c3>0.
	//
	//			tau and c values are assumed to be constant
	//			throughout the simulation. The values returned
	//			during setInitialState processing are cached
	//			for efficiency and used throughout. Q10 adjustment
	//			are applied to tau1, tau2, and tau3.
	//
	//			Note that using s1/tau1 results in s1, s2, and s3 being 
	//			dimensionless. This form is useful for scaling to a 
	//			consistent range of variable values.
	// ----------------------------------------------------------------

	class TripleExpSynapticCond : public SingleExpSynapticCond {

	public:

		// Constructors and destructor
		TripleExpSynapticCond(Number gMaxValue=0);
		virtual ~TripleExpSynapticCond();

		// Accessors
		inline  Number			s2() { updateCachedValues(); return _s2Now; }
		inline  Number			s3() { updateCachedValues(); return _s3Now; }

		// Set synapse conductance multiplier based on maximum value
		virtual void			gMax(Number maxCond);

		// Return the total conductance for synapses of this type
		virtual Number			conductance();

		// Reporting interface ----------------------------------------

		virtual int				numInternalStateVar() { return 3; }
		virtual Number			internalStateValue(int k);

		virtual const char*		componentName() {
			return "TripleExpSynapticCond"; }	// subclass resp.

		virtual const char**	stateLabels() {
			static const char* slab[] = {"s1","s2","s3"}; return slab; }

		// Subclass responsibilities ----------------------------------

		// Parameter accessors

		// Time constant at channel's rated temperature (before Q10 adjustment).
		// Tau1 (via superclass) is rising rate. Tau2 and Tau3 are falling rates.
		virtual SimTime			tau2() = 0;		// falling time constant
		virtual SimTime			tau3() = 0;		// falling time constant

		virtual Number			c2() = 0;		// weight constant (c3 is implied)

	protected:
		Number					_s2Start;		// starting state value
		Number					_s2Now;			// current state value
		Number					_s3Start;		// starting state value
		Number					_s3Now;			// current state value

		// Advance the current state values for the time step size in _H
		virtual void			advanceState();

		// Clear the state to initial zero values
		virtual void			clearState();

		// Save the current state as the starting state for a time step
		virtual void			saveStartingState();

		// Restore the current internal state to time step start values
		virtual void			restoreToStartingState();

		// Update cached values to reflect change in time step
		virtual void			updateCacheForStepSize();

		// Update current conductance to reflect an AP event
		virtual void			updateForAPEvent(ActionPotentialEvent* apEvent);

		// Address class cache data as used by a superclass
		virtual SingleExpSynapticCondClassCache* pSingleExpClassCache()
								{ return pTripleExpClassCache(); }

		// Subclass responsibilities ----------------------------------

		// Address class cache data allocated by a subclass.
		virtual TripleExpSynapticCondClassCache* pTripleExpClassCache()=0;

		// Initialize the class cache (subclass may extend this)
		virtual void			initializeClassCache();
	};

	// ----------------------------------------------------------------
	// CLASS:	PlasticityRule
	// EXTENDS:	ModelComponent
	// DESC:	Abstract class for policies to adjust neuron
	//			properties in response to synaptic activity.
	// RESP:
	//		1.	Act on presynaptic spikes (via subclass)
	//		2.	Act on postsynaptic events and state (via subclass)
	//		3.	Know offset to plasticity data within synapse
	//		4.	Create and destroy synapse state data (via subclass)
	//
	// NOTES:	As a model component, subclasses must provide
	//			implementations for the required functions.
	//			However, whether a plasticity rule object is
	//			itself a component of any model is optional.
	//			Because there are requirements for ordering
	//			of processing, rules are typically notified
	//			of events indirectly via the owning synaptic
	//			conductance objects.
	//
	//			Plasticity rules are implemented as a policy
	//			to allow changes of rules without affecting
	//			other aspects of the class hierarchy.
	//
	//			isEnabled and isDisabled are provided so that a 
	//			rule can be temporarily enabled or disabled.
	//			When disabled, the rule does not receive
	//			finalizeAPEvent or applyEndOfStep invocations
	//			unless notifyRule() is overridden. isEnabled
	//			may need to be overridden if internal states
	//			must be updated to reflect rule activity state.
	// ----------------------------------------------------------------

	class PlasticityRule : public ModelComponent {
	
	public:

		// Constructors and destructor
		PlasticityRule();
		virtual ~PlasticityRule();

		// Accessor for associated synaptic conductance
		inline  SynapticConductance*	synapticCond() { return _synapticCond; }
		virtual void					synapticCond(SynapticConductance* sc);

		// Accessor for state offset within synapse
		inline  unsigned int	stateOffset() { return _stateOffset; }
		virtual void			stateOffset(unsigned int n) { _stateOffset=n; }

		// Accessor for current activity state of the rule
		inline  bool			isEnabled() { return _isEnabled; }
		inline  bool			isDisabled() { return !_isEnabled; }

		virtual void			isEnabled(bool x) { _isEnabled=x; }
		virtual void			isDisabled(bool x) { _isEnabled=!x; }

		// Get address for state data for synapse (subclass must cast to use)
		inline  Byte*			plasticityData(Synapse* syn) 
								{ return reinterpret_cast<Byte*>(syn)+stateOffset(); }

		// Accessors for current synapse weight
		inline  Number			synapticWeight(Synapse* syn) 
								{ return synapticCond()->synapseWeight(syn); }

		inline  void			synapticWeight(Synapse* syn, Number w) 
								{ synapticCond()->synapseWeight(syn,w); }

		// Get membrane potential and derivative in containing compartment
		inline  Number			Vm()	{ return synapticCond()->Vm(); }
		inline  Number			VmDot()	{ return synapticCond()->VmDot(); }

		// Get the time step size (at end of step or before model assigned)
		inline	SimTime			timeStepSize() 
		{ return model()==NULL ? 0 : model()->timeStepSize(); }

		// Subclass responsibilities ----------------------------------

		// Return the size to allocate for plasticity state data.
		virtual unsigned int	plasticityStateSize() = 0;

		// Create a plasticity state object in a synapse buffer.
		// This must be overridden if plasticityStateSize > 0.
		virtual void			createPlasticityState(Synapse* syn, Number wght);

		// Destroy a plasticity state object in a synapse buffer.
		// Subclass must override if destructor processing is required.
		virtual void			destroyPlasticityState(Synapse* syn) {};

		// By default, there are no separate state variables for this object.
		virtual int				numStateVar() { return 0; }

		// Default is to not register with model. but subclass can override. 
		// applyEndOfStep is invoked via synaptic conductance in any case.
		virtual bool			joinModelAsComponent() { return false; }

		// Return true if the rule is to receive notification of events.
		// Default is to do so only if the rule is currently enabled.
		virtual bool			notifyRule() { return isEnabled(); }

		// Do any one-time updates of an AP event as it is dequeued for processing.
		// This is typically invoked for presynaptic plasticity rules.
		virtual void			finalizeAPEvent(ActionPotentialEvent* apEvent) {}

		// Take action at the end of the time step before AP purge.
		// Default is to apply the rule one event at a time in order.
		// Subclass may want to get control before or after this process.
		// Also, subclass may want to override if this is a no-op.
		virtual void			applyEndOfStep(ActionPotentialEventQueue& apQueue);

		// Update state one action potential event at a time.
		virtual void			applyRule(ActionPotentialEvent* apEvent) {}

		// Apply any neuromodulation parameters (default is no-op)
		virtual void			setModParams(
			TokenId				modId,			// modulation type token id
			int					numParams,		// number of parameter values
			Number*				params) {}		// parameter values as an array

		// Return a component name for reporting (should be overridden)
		virtual const char*		componentName() {return "PlasticityRule"; }

		// Return an identifier of the type of plasticity type
		virtual TokenId			plasticityTypeId() {return NullTokenId; }

	protected:
		SynapticConductance*	_synapticCond;	// owning synaptic conductance
		unsigned int			_stateOffset;	// offset to plasticity state
		bool					_isEnabled;		// flag indicating currently active
	};

	// ----------------------------------------------------------------
	// CLASS:	RandomReleaseRule
	// EXTENDS:	PlasticityRule
	// DESC:	Implements a presynaptic plasticity rule based
	//			on random release of neurotransmitter.
	// RESP:
	//		1.	Access release probability from synaptic conductance
	//		2.	Set event quantity to 1 or 0 based on a random value.
	// ----------------------------------------------------------------

	class RandomReleaseRule : public PlasticityRule {
	
	public:
		// Constructors and destructor
		RandomReleaseRule(Number probRel=-1);
		virtual ~RandomReleaseRule();

		// Access release probability for this rule.
		// If the release probability specified here is <0, then
		// the owning synaptic conductance value is used instead
		// when determining whether a release should occur or not.
		inline  Number			releaseProbability() { return _releaseProbability; }
		virtual void			releaseProbability(Number p) { _releaseProbability=p; }

		// Set the event quantity based on random selection
		virtual void			finalizeAPEvent(ActionPotentialEvent* apEvent);

		// There is no plasticity state data used here
		virtual unsigned int	plasticityStateSize() { return 0; }

		// No end of step processing is needed
		virtual void			applyEndOfStep(ActionPotentialEventQueue& apQueue) {}

	protected:
		Number					_releaseProbability;	// rule specific Prel
	};

	// ----------------------------------------------------------------
	// CLASS:	Synapse
	// EXTENDS:	none
	// DESC:	Class for representing synaptic
	//			connections, Both pre- and post-synaptic
	//			portions of the synapse are combined into
	//			a one object as a default.
	// RESP:
	//		1.	Know associated postynaptic response
	//		2.	Know presynaptic axon process.
	//		3.	Know delay between afferent firing and 
	//			action potential being effective here.
	//		4.	Pass action potential from pre-synaptic axon
	//			to post-synaptic response object.
	//
	// NOTES:	This class provides the header for a buffer which
	//			holds other synaptic state information. It is not
	//			expected that this class will be subclassed and thus
	//			virtual functions are not used. This reduces storage
	//			requirements by avoiding virtual function pointers.
	//			Such otherwise small optimizations are justified by
	//			the large number of synapses needed when simulating 
	//			a large network with many connections.
	//
	//			State data used by SynapticResponse subclasses is
	//			appended following the area occupied by a Synapse
	//			instance. When multiple responses are used within
	//			a single Synapse instance, their states are placed
	//			one after the other in storage and the responses 
	//			are given different offsets from which to locate 
	//			the data unique to them. Destruction of any state 
	//			data is done by the postsynaptic response object.
	// ----------------------------------------------------------------

	class Synapse {
	
	public:

		// Constructors and destructor
		Synapse(
			AxonProcess* axon,							// presynaptic axon process
			SynapticResponse* resp,						// postsynaptic response process
			Number dist = 100*UOM::micron);				// axonal distance
		
		~Synapse();

		// Override the normal new function to allow placement at a
		// specific location where storage has already been allocated.
		void*						operator new(unsigned int, void *buf) {return buf;}

		// Provide a delete function to be invoked if the above new fails (not likely)
		void						operator delete(void* syn, void* buf ) {}

		// Override the normal delete function to allow deletion of the
		// whole buffer containing synapse plus any state data.
		void						operator delete(void* buf) {
									delete[] reinterpret_cast<Byte*>(buf);}

		// Accessors
		inline  Synapse*			nextPresynaptic() { return _nextPresynaptic; }
		inline  Synapse*			nextPostsynaptic() { return _nextPostsynaptic; }

		inline  AxonProcess*		presynapticProcess() { return _presynapticProcess; }
		inline  SynapticResponse*	postsynapticProcess() { return _postsynapticProcess; }

		inline  Number				delay() { return _delay; } // propagation delay
		
		// Act on an action potential from the presynaptic neuron
		virtual void				signalActionPotential(ActionPotentialEvent* apev);

		// List manipulation functions --------------------------------
		
		// Add the current synapse before another one already in a list
		void						addBeforeInPresynapticList(Synapse* synapseToFollow);
		void						addBeforeInPostsynapticList(Synapse* synapseToFollow);

		// Remove the synapse following this one from a list
		void						removeNextFromPresynapticList();
		void						removeNextFromPostsynapticList();

		// Clear links in preparation for deleting a whole list.
		// Return a link to the next synapse in the list.
		// Delete this object if both pre and post synaptic owners have been cleared.
		Synapse*					clearPresynaptic();
		Synapse*					clearPostsynaptic();

	protected:
		Synapse*					_nextPresynaptic;		// next in presynaptic list
		Synapse*					_nextPostsynaptic;		// next in postsynaptic list
		AxonProcess*				_presynapticProcess;	// owning axon process
		SynapticResponse*			_postsynapticProcess;	// owning synaptic response
		Number						_delay;					// propagation delay time
	};

	// ----------------------------------------------------------------
	// CLASS:	SynapseState
	// EXTENDS:	none
	// DESC:	Class for representing synaptic state data
	//			for an instance of the response object.
	//			Subclasses can extend the basics supplied here.
	// RESP:
	//		1.	Know synaptic weight.
	//
	// NOTE:	This is a minimal version of state. Subclasses
	//			may extend state as required.
	//
	//			It is critical that size be correctly provided during
	//			object creation when this is subclassed. See the
	//			appropriate createSynapse logic. Similarly, if any
	//			destructor processing is needed, either
	//			destroySynapseState or destroyOtherSynapseState must
	//			be overridden as appropriate. See classes
	//			SynapticResponse and SynapticConductance.
	//
	//			This class is intended to provide data only. 
	//			It is one of the few classes where virtual 
	//			functions must not be used. This restriction is
	//			imposed to reduce synapse storage usage.
	// ----------------------------------------------------------------

	class SynapseState {
	
	public:
		Number						weight;	// synaptic weight (unitless)
	};

}; // end of namespace

#endif // #ifndef __BSNF_NMOD_H_
