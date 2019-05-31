// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_sim.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the basic classes used as a
// framework for building biological network simulations.
// As is normal for frameworks, many of these classes are
// abstract and subclassed as needed for the simulation.
//
// This header declares classes used in all simulations.
// These classes implement a generalized implementation
// framework for event driven or ODE-based simulations.
// Some specializations are necessary because of C++
// restrictions -- see ModelComponent below.


// Only include this header once per compilation unit
#ifndef __BSNF_SIM_H_
#define __BSNF_SIM_H_


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
#include <vector>
#include <queue>
#include <map>

// Header for use in manipulating C strings
#include <cstdio>

// Incorporate all the names from the std library by reference
using namespace std;

// Required BNSF headers
#include "bnsf_base.h"
#include "bnsf_math.h"


// ====================================================================
// Primary namespace for the framework
// ====================================================================


namespace BNSF {


	// ================================================================
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ================================================================

	// Simulator classes
	class Model;
	class ModelEntityClassCache;
	class ModelEntity;
		class ModelComponent;

	class Controller;
	class Solver;
		class ODESolver;
			class EulerMethodSolver;
			class RungeKuttaSolver;
			class AdaptiveSolver;
			class ClockSolver;
	class SolverList;
	class SolverListRunOrderFunct;

	class Event;
	class EventOrderFunct;

	// External interface classes
	class Probe;
		class ExternalRecorder;
			class ExternalStateRecorder;

	// ================================================================
	// Common typedefs, enums, and globals
	// ================================================================

	// Define a type for simulation time. Simulation time is basically
	// a real number, but the precision can be changed if necessary
	// to avoid loss of accuracy. Note that common time units such
	// as 1 msec may not have exact binary floating point equivalents.

	typedef double			SimTime;

	// Time values at plus and minus infinity
	static const SimTime	InfiniteFuture = numeric_limits<SimTime>::infinity();
	static const SimTime	InfinitePast = -numeric_limits<SimTime>::infinity();
	
	// Tolerance to comparing two time values as equal (arbitrary small value)
	static const SimTime	EpsilonTime = UOM::picosec;

			
	// ================================================================
	// Vector and iterator typedefs for (most) classes
	// Note that by convention pointers are stored in
	// vectors so this info is not in the typedef name.
	// ================================================================

	typedef vector<Model*>							ModelVector;
	typedef ModelVector::iterator					ModelVectorIt;

	typedef vector<ModelEntity*>					ModelEntityVector;
	typedef ModelEntityVector::iterator				ModelEntityVectorIt;

	typedef vector<ModelComponent*>					ModelComponentVector;
	typedef ModelComponentVector::iterator			ModelComponentVectorIt;

	typedef vector<SolverList*>						SolverListVector;
	typedef SolverListVector::iterator				SolverListVectorIt;

	typedef priority_queue<SolverList*,SolverListVector,SolverListRunOrderFunct>
													SolverListPriorityQueue;

	typedef map<SimTime,SolverList*>				SolverListMap;
	typedef SolverListMap::iterator					SolverListMapIt;
	typedef SolverListMap::value_type				SolverListMapValueType;

	typedef vector<Event*>							EventVector;
	typedef EventVector::iterator					EventVectorIt;

	typedef priority_queue<
		Event*,
		EventVector,
		EventOrderFunct>							EventPriorityQueue;

	typedef vector<Probe*>							ProbeVector;
	typedef ProbeVector::iterator					ProbeVectorIt;

	typedef pair<Probe*,ModelComponent*>			ProbeTargetPair;
	typedef vector<ProbeTargetPair>					ProbeTargetVector;
	typedef ProbeTargetVector::iterator				ProbeTargetVectorIt;

	// Timestep type for localStateUpdate when using Crank-Nicolson solver.
	// Because some objects are offset in time from others in this method,
	// special half steps are needed to get the process started and to
	// reach a common time point at the end.
	enum CNStepType {
		CNFullStep,					// Full time step
		CNStartingHalfStep,			// Half step when starting
		CNEndingHalfStep			// Half step when ending
	} ;
			
	// ================================================================
	// Declare functionals that apparently need to be
	// declared before they are referenced (mere prototypes
	// are not enough).
	// ================================================================

	// Order function used in ordering ODE solver lists for running
	bool solverListRunOrder(SolverList* rhs, SolverList* lhs);

	// Define a functional wrapping the order function
	// permiting its use as a parameter in template,
	class SolverListRunOrderFunct 
		: public binary_function<SolverList*, SolverList*, bool> {
	public:
		bool operator()(SolverList* _X, SolverList* _Y) const
			{return solverListRunOrder(_X, _Y); }
	};

	// Order function used in ordering event queues 
	bool eventOrder(Event* rhs, Event* lhs);

	// Define a functional wrapping the order function
	// permiting its use as a parameter in a template,
	class EventOrderFunct : public binary_function<Event*, Event*, bool> {
	public:
		bool operator()(Event* _X, Event* _Y) const
			{return eventOrder(_X, _Y); }
	};


	// ================================================================
	// Classes that support the infrastructure of the simulation
	// ================================================================


	// ----------------------------------------------------------------
	// CLASS:	Model
	// EXTENDS:	none
	// DESC:	Class containing for top level components
	//			in a simulation. Objects in this class
	//			hierarchy hold information common to all
	//			solution methods so that solution methods
	//			can be changed as needed.
	// RESP:
	//		1.	Know size of state vector
	//		2.	Add components to state vector
	//		3.	Sort components based on dependencies
	//		4.  Interface with ODE solver
	//		5.	Provide access to a stream of random numbers
	//
	// NOTE:	Models do not own their components, that is,
	//			components should remove themselves when
	//			they are destroyed but destroying the
	//			model does not imply destroying the
	//			components.
	// ----------------------------------------------------------------

	class Model {

	public:

		// Constructors and destructor
		Model();
		virtual ~Model();

		// Accessors
		inline int				stateVectorSize() { return _stateVectorSize; }

		inline  ODESolver*		solver() { return _solver; }
		virtual void			solver(ODESolver* newSolver);

		inline  Number*			stateVector() { return _stateVector; }		
		inline  Number*			weightVector() { return _weightVector; }
		inline  Number*			derivVector() { return _derivVector; }

		inline  SimTime			currentTime() { return _currentTime; }
		virtual void			currentTime(SimTime t) { _currentTime = t; }

		inline  SimTime			timeStepStart() { return _timeStepStart; }

		// Access the step size so far in the current step
		inline  SimTime			currentStepSize() { return _currentTime - _timeStepStart; }

		// Access time step size and changed flag at the end of a step
		inline  SimTime			timeStepSize()		{ return _timeStepSize; }
		inline  bool			stepSizeChanged()	{ return _stepSizeChanged; } 

		// Provide direct access to components (use with caution)
		inline  ModelComponentVector& components() { return _components; }

		// Access a generator for uniform random numbers. This allows
		// one or more models to use a common stream of random numbers
		// while other models use different independent streams.
		// If no generator has been supplied, use the default generator.
		inline  UniformRandom*	uniformRandom() 
		{ return _uniformRandom==NULL ? UniformRandom::defaultGen() : _uniformRandom; }

		virtual void			uniformRandom(UniformRandom* unif) { _uniformRandom = unif; }

		// Maintain the collection of probes
		virtual void			addProbe(Probe* pr, ModelComponent* mc);
		virtual void			removeProbe(Probe* pr, ModelComponent* mc);

		// Maintain the collection of associated components
		virtual void			addComponent(ModelComponent*);
		virtual void			removeComponent(ModelComponent*);
		virtual void			clearComponents();

		// ODE Solver interface functions ------------------
		virtual void			simulationStarted();
		virtual void			prepareTimeStepEnded(SimTime timeNow);
		virtual void			timeStepEnded();
		virtual void			simulationEnded() {} // for now a no-op

	protected:
		int						_stateVectorSize;	// size of the state vector to build
		ModelComponentVector	_components;		// current components of the model
		ODESolver*				_solver;			// current solver in use
		ProbeTargetVector*		_probesPtr;			// model-level probes
		UniformRandom*			_uniformRandom;		// random number source

		// State information provided by the ODE solver
		Number*					_stateVector;		// current state vector
		Number*					_derivVector;		// current d/dt vector
		Number*					_weightVector;		// weights used for adaptive methods
		
		SimTime					_timeStepStart;		// beginning of the current step
		SimTime					_currentTime;		// current time in the step

		// Values of time step and flag indicating change this step.
		// These are valid only at the end of the time step.
		SimTime					_timeStepSize;		// size of time step just ended
		bool					_stepSizeChanged;	// size of step was just changed

		virtual ProbeTargetVector&  probes();		// lazy init of probesPtr

		virtual void			allocateVectors();	// allocate state vector etc.
		virtual void			assignSVOffsets();	// set state vector offsets
		virtual void			orderComponents();	// order by dependency
	};


	// ----------------------------------------------------------------
	// CLASS:	ModelEntityClassCache 
	// EXTENDS:	none
	// DESC:	Abstract class providing basic protocol for caching
	//			values associated with all instances of a model entity.
	//			Typical use is to cache values derived from time step
	//			size when this will generally be the same among many
	//			instances of a given class. Since speed is of the
	//			essence here, direct access to cache data is allowed.
	//
	// NOTE:	At present no attempt is made to provide thread safe
	//			access to cache data. This may need to be changed in
	//			the future. One possible approach is to allocate a set
	//			of class caches for each model.
	// ----------------------------------------------------------------

	class ModelEntityClassCache {

	public:
		// Constructors and destructor
		ModelEntityClassCache() : isInitialized(false), H(-1) {}
		~ModelEntityClassCache() {}

		bool					isInitialized;	// initialization done
		SimTime					H;				// size of current time step
	};

	// ----------------------------------------------------------------
	// CLASS:	ModelEntity
	// EXTENDS:	none
	// DESC:	Abstract superclass of modelled objects.
	// RESP:
	//		1.	Declare common protocols
	//		2.	Handle event dispatch as recipient (via subclass)
	//
	// NOTE:	This class is primarily here for future expansion.
	//			Dispatching events involves parallel subclasses here
	//			and also in the Event hierarchy.
	//
	//			A useful enhancement would be to allow generic events
	//			to be queued and dispatched at the appropriate time. 
	// ----------------------------------------------------------------

	class ModelEntity {

	public:
		// Constructors and destructor
		ModelEntity();
		virtual ~ModelEntity();

		// Dispatch a generic received event (see note above).
		// Subclass must override this if it is invoked.
		virtual void			dispatchEvent(Event* ev);
	};

	// ----------------------------------------------------------------
	// CLASS:	ModelComponent
	// EXTENDS:	ModelEntity
	// DESC:	Abstract superclass of dynamical system objects.
	// RESP:
	//		1.	Know model
	//		2.	Know number of state variables
	//		3.	Know offset into overall state vector
	//		4.	Provide interface for ODESolver actions
	//		5.	Support change notifications via subscription
	//		6.	Access random numbers via the model random stream
	//
	// NOTE:	Ideally, this class should be split into a pure simulation
	//			class and a subclass that knows about neurons, channels, etc.
	//			However, this may be hard to do efficiently in C++ because
	//			of the interface with ODE solvers.
	//
	//			Random number access is a pass-through using the
	//			base uniform random number generator held by the
	//			model. Access in this fashion is a utility function
	//			that may be less efficient than building separate
	//			random number objects for each type used.
	// ----------------------------------------------------------------

	class ModelComponent : public ModelEntity {

	public:

		// Typedefs and enums for this class --------------------------

		// Commonly used reasons for signaling changes
		enum {
			terminatedChange=-1,			// Object destroyed
			otherChange=0,					// Other uncategorized change (default)
			parameterChange=1,				// Associated parameter changed
			stateChange=2					// Internal state changed
		} ChangeReason;

		// Constructors and destructor --------------------------------

		ModelComponent(Model* m=NULL);
		virtual ~ModelComponent();

		// Accessors
		inline  Model*			model() { return _model;}
		virtual void			model(Model* m);
		
		virtual ODESolver*		solver() { return model()->solver(); }

		inline  SimTime			currentTime() { return _model->currentTime(); }
		inline  SimTime			timeStepStart() { return _model->timeStepStart(); }

		inline  Number*			stateVector() { return _model->stateVector(); }
		inline  Number*			derivVector() { return _model->derivVector(); }
		inline  Number*			weightVector() { return _model->weightVector(); }

		inline  int				svOffset() { return _svOffset; }
		virtual void			svOffset(int n) {_svOffset=n; _svArray=stateVector()+n; }

		// Accessors for numeric id. subclasses must provide any implementation.
		virtual int				numericIdentifier() { return -1; } // -1 = none
		virtual void			numericIdentifier(int id) {} // placeholder only

		// Access the slots for state vector and derivative
		inline  Number&			stateValue(int n) { return _svArray[n]; }
		inline  Number&			derivValue(int n) { return derivVector()[_svOffset+n];}
		inline  Number&			weightValue(int n) { return weightVector()[_svOffset+n]; }

		// Access a boolean indicating if this object should register with a model
		// to receive events and state vector space or not.
		virtual bool			joinModelAsComponent() { return true; }

		// Functions for ordering ODE function invocations among components.
		// Prerequisites can be returned as a vector or as a single item.
		// Subclasses should override prerequisite() for a single item and
		// prerequisites() if more than one are to be returned.
		// In either case, hasPrerequisites() must be overridden.
		virtual bool			hasPrerequisites() { return false; }
		virtual ModelComponent*	prerequisite() { return NULL; }
		virtual ModelComponentVector prerequisites(); // default is zero or one prereqs

		// Functions for propagating changes --------------------------

		// Return the vector holding subscribers. Subclasses responsibility if used.
		virtual ModelComponentVector*	subscribers() { return NULL; }

		// Add/remove a subscriber to changes from this object
		virtual void			addSubscriber(ModelComponent* mc);
		virtual void			removeSubscriber(ModelComponent* mc);

		// Notify subscribers that this object has changed in some way.
		// The meaning of reason is by agreement between this object and subscribers.
		// See changeReason typedef for commonly used reasons.
		virtual void			changed(int reason=otherChange);

		// Receive notification that a subscribed object has changed.
		virtual void			updateFrom(ModelComponent* mc, int reason) {}
		
		// Solver interface functions ---------------------------------

		// Return number of ODE state vectors to reserve for this component
		virtual int				numStateVar() = 0;		// subclass responsibility

		// Set default initial state values
		virtual void			setInitialState() {};	// subclass responsibility if used

		// Set weight vector for ODE solver -- default is to set all values to 1
		virtual void			setWeightValues();

		// Index of this component in a Jacobian matrix
		//.This is a default only and the matrix may be reordered.
		virtual int				jacobianIndex() { return _svOffset; }
		virtual void			jacobianIndex(int n);	// subclass responsibility if used

		// Identifier of the domain of this component in the
		// sense of domain decomposition ODE solver methods. 
		// These values are used to group related items in the 
		// state vector. Subclasses are responsible for providing 
		// any non-default value. Domains 0 through 4095 are 
		// reserved for BNSF components.
		virtual unsigned int	domain() { return 0; }

		// Update derviatives vector with time derivatives of state.
		// This must be overridden by subclass if invoked by solver.
		virtual void			computeDerivatives();

		// Special inquiries from model, ODE solver, or probes
		virtual bool			raisesEvents() {return false; } // curently unused
		virtual bool			updatesDerivVector() { return numStateVar()>0; }

		// Options to request notification on different ODE solver events
		// so that only components requesting these notifications will get them
		virtual bool			notifyOnStateVectorChanged() { return false; }
		virtual bool			notifyOnTimeStepEnded() { return false; } 

		// Notifications that a solver can make (default is to ignore)
		virtual void			simulationStarted() {}
		virtual void			timeStepEnded() {}
		virtual void			simulationEnded() {}
		virtual void			stateVectorChanged() {}

		// Add this component to a controller using a model and solver
		// that have already been allocated and associated with this object.
		virtual void			addToController(Controller* cont);

		// Interface functions primarily for Crank-Nicolson solvers ---

		// Update state vector by time step h using local information.
		virtual void			localStateUpdate(
			SimTime				h,				// time to advance
			CNStepType			stepType);		// type of step

		// Indicate whether or not the object is capable of a local update.
		virtual bool			canPerformLocalUpdate() { return false; }

		// Indicate whether or not the local update function is at
		// least second order accurate in time.
		virtual bool			isSecondOrderAccurate() { return false; }

		// Specify how this component is treated in a staggered
		// evaluation scheme such as Crank-Nicolson. The meaning of
		// offset  (versus aligned) must be consistently applied but
		// is otherwise by agreement among components of the model.
		virtual bool			isOffsetForCN() { return false; }
		inline  bool			isAlignedForCN() { return !isOffsetForCN(); }

		// Accessors for probing state values -----------------------
		
		// Add/remove a probe for this specific component. 
		// Model must already be defined before probes are added or removed.
		virtual void			addProbe(Probe* pr);
		virtual void			removeProbe(Probe* pr);

		// Add to a vector of the components to be included when
		// probing the contents of this component. By default only 
		// this component is added in the vector, though subclasses
		// may extend this to include related items.
		virtual void			addToComponentsToProbe(ModelComponentVector& comps);

		// Return number of state vector entries included in the object itself
		// For reporting, internal state values conceptually occur after
		// after ODE state values in the sequence of labels and units of measure.
		virtual int				numInternalStateVar() { return 0; }		// default

		// Return number of ODE and internal state variables
		virtual int				numAllStateVar() { return numStateVar()+numInternalStateVar(); }

		// Return the value of an internal state value.
		// n is index of the variable with range 0 through numInternalStateVar-1.
		// Subclass must override this if numInternalStateVar>0.
		virtual Number			internalStateValue(int n);

		// Return a name for the component itself and also its states
		// Actual values are provided by subclasses. 
		// Only defaults are provided here.
		virtual const char*		componentName();
		virtual const char**	stateLabels();

		// Return the component name as a token. Subclasses may want to
		// override this for efficiency.
		virtual TokenId			componentId() { return token( componentName() ); }

		// Return a unit of measure to use in external outputs
		virtual Number*			unitsOfMeasure();

		// Convert a state value for external output
		virtual Number			stateValueConverted(int stateValueOffset);

		// Access random numbers via the model stream -----------------

		// Get the generator from the model
		inline  UniformRandom*	uniformRandom() { return model()->uniformRandom(); }

		// Get a uniformly distributed random number
		inline  double			runif() { return uniformRandom()->next(); }
		inline  double			runif(double min, double max) 
		{ return uniformRandom()->next(min,max); }

		// Get a normally distributed random number given mean and std deviation
		inline  double			rnorm(double mean=0, double sdev=1)
		{ return NormalRandom::value(mean,sdev,uniformRandom()); }

		// Get an exponentially distributed random number given the mean
		inline  double			rexp(double mean=1)
		{ return ExponentialRandom::value(mean,uniformRandom()); }

		// Get a Poisson distributed random number given the mean
		inline  int				rpois(double mean)
		{ return PoissonRandom::ivalue(mean,uniformRandom()); }

		// Get a binomial distributed random number given population size and prob
		inline  int				rbinom(int size, double prob)
		{ return BinomialRandom::ivalue(size,prob,uniformRandom()); }

		// Inquiries specific to neural simulations -----------------

		virtual bool			isNeuron() {return false; }
		virtual bool			isCompartment() { return false; }
		virtual bool			isCalciumPool() { return false; }
		virtual bool			isIonChannel() { return false; }
		virtual bool			isGateVariable() { return false; }
		virtual bool			isVoltageDepTabChannel() { return false; }
		virtual bool			isSynapticResponse() { return false; }
		
	protected:
		Model*					_model;			// owning model
		Number*					_svArray;		// pointer to stateValue(0)
		int						_svOffset;		// offset in full state vector

	private:
		// Set initial default values for use in constructor
		void					initialize();
	};

	// ----------------------------------------------------------------
	// CLASS:	Controller
	// EXTENDS:	none
	// DESC:	Controls scheduling of ODE solvers to ensure
	//			that all solvers are invoked in a consistent
	//			time order with the evaluations effectively
	//			interleaved over time.
	// RESP:
	//		1.	Maintain lists of ODE solvers by time step.
	//		2.	Invoke ODE solution methods in order.
	//
	// NOTE:	There is an assumption that ODE solvers use
	//			a common collection of time steps, that is,
	//			the number of time step values is less than
	//			the number of ODE solvers. This means that
	//			scheduling is for a group of ODE solvers all
	//			of which share the same time step. ODE solvers
	//			can, however, dynamically change their time 
	//			steps in which case they are moved from one
	//			group to another and rescheduled as needed.
	// ----------------------------------------------------------------

	class Controller {

	public:

		// Constructors and destructor
		Controller();
		virtual ~Controller();

		// Accessors
		inline SimTime				beginTime() { return _beginTime; }
		virtual void				beginTime(SimTime t) { _beginTime = _evalTime = t; }

		inline SimTime				endTime() { return _endTime; }
		virtual void				endTime(SimTime t) { _endTime = t; }

		inline SimTime				evalTime() { return _evalTime; }

		// Add an ODESolver to this controller
		virtual void				addSolver(Solver* solver);

		// Add a model based on its solver
		virtual void				addModel(Model* m);

		// Add a model component based on its model's solver
		virtual void				addModelComponent(ModelComponent* mc);

		// Polymorphic add for convenience
		virtual void				add(Solver* s) { addSolver(s); }
		virtual void				add(Model* m) { addModel(m); }
		virtual void				add(ModelComponent* mc) { addModelComponent(mc); }

		// Update the step time of an existing ODESolver
		virtual void				changeSolverTimeStep(Solver* s, SimTime oldStep);

		// Simulation actions
		virtual void				start();		// initialize -- prepare to run
		virtual void				run();			// evaluate to end time
		virtual void				runUpTo(SimTime t);			// run to an end time
		virtual void				runForDuration(SimTime dt);	// run for a time period
		virtual void				finish();		// wrap up -- finalize and disconnect

	protected:
		SimTime						_beginTime;		// starting time
		SimTime						_evalTime	;	// time we have evaluated up to
		SimTime						_endTime;		// ending time for this run
		SolverListMap				_solvers;		// solvers by time step
		SolverListPriorityQueue		_runnableQueue;	// solvers ready to run

	};

	// ----------------------------------------------------------------
	// CLASS:	SolverList
	// EXTENDS:	none
	// DESC:	Provides a header for a list of Solvers.
	// RESP:
	//		1.	Access first and last entries
	//		2.	Add new entry at first or last position
	//		3.	Remove a node of this list
	//		4.	Order entries by time step (via solverListOrder)
	//		5.	Run all solvers in the list up to a specified time
	//
	// NOTE:	Much of this function could be accomplished
	//			using std template lists, but some of the
	//			operations are just easier to understand when
	//			done in the normal way and certainly no slower.
	// ----------------------------------------------------------------

	class SolverList {

	public:
	
		// Constructors and destructor
		SolverList();
		virtual ~SolverList();

		// Accessors
		inline  bool		isEmpty() { return _firstNode == NULL; }
		inline  Solver*		firstNode() { return _firstNode; }
		inline  Solver*		lastNode() { return _lastNode; }

		inline  SimTime		timeStep() { return _timeStep; }
		virtual void		timeStep(SimTime h) { _timeStep = h; }

		inline  SimTime		evalTime() { return _evalTime; }
		virtual void		evalTime(SimTime t) { _evalTime = t; }

		inline	bool		isRunnable() { return _isRunnable; }
		virtual void		isRunnable(bool val) { _isRunnable = val; }

		// Add nodes at beginning or end of list
		virtual void		addFirst(Solver* node);
		virtual void		addLast(Solver* node);

		// Remove a node from this list
		virtual void		removeNode(Solver* node);

		// Run all the solvers in the list
		virtual void		runAllUpTo(SimTime endTime);

	protected:
		Solver*				_firstNode;			// first in the list or NULL
		Solver*				_lastNode;			// last in the list or NULL
		SimTime				_timeStep;			// time step of nodes in the list
		SimTime				_evalTime;			// time of completed evaluations
		bool				_isRunnable;		// is present on the runnable queue
	};


	// ----------------------------------------------------------------
	// CLASS:	Solver
	// EXTENDS:	none
	// DESC:	Abstract class for solvers used to
	//			numerically evaluate dynamical systems.
	// RESP:
	//		1.	Declare common protocol for controller interface.
	//		2.	Add and remove from a list of similar nodes.
	//		3.	Know start and ending simulation times
	//		4.	Know current time step size
	//		5.	Round time steps to a limited set of values
	// ----------------------------------------------------------------

	class Solver {

	public:
			
		// Constructors and destructor
		Solver();
		virtual ~Solver();

		// Accessors
		inline Solver*			nextNode() { return _nextNode; }
		inline Solver*			prevNode() { return _prevNode; }

		inline  Controller*		controller() { return _controller; }
		virtual void			controller(Controller* cont) { _controller = cont; }

		inline  SolverList*		solverList() { return _solverList; }
		virtual void			solverList(SolverList* list) { _solverList = list; }

		inline  SimTime			beginTime() { return _beginTime; }
		virtual void			beginTime(SimTime t) { _beginTime = _currentTime = t; }

		inline  SimTime			endTime() { return _endTime; }
		virtual void			endTime(SimTime t) { _endTime = t; }

		inline  SimTime			timeStep() { return _timeStep; }
		virtual void			timeStep(SimTime t);

		inline  SimTime			currentTime() { return _currentTime; }

		virtual bool			hasController() { return _controller!=NULL; }

		// Set a time step size from a restricted list of possible
		// values and return the value selected.
		virtual SimTime			timeStepLessOrEqual(SimTime h);

		// List manipulation routines
		virtual void			insertBefore(Solver* node);
		virtual void			insertAfter(Solver* node);
		virtual void			removeFromList();

		// Run up to a new end time
		virtual void			runUpTo(SimTime t);

		// Basic simulation actions (subclass responsibility)
		virtual void			start()=0;			// initialize
		virtual void			run()=0;			// run the simulation
		virtual void			finish()=0;			// clean up - finalize

	protected:
		Solver*					_nextNode;			// next in doubly linked list
		Solver*					_prevNode;			// previous in doubly linked list
		Controller*				_controller;		// processing scheduler
		SolverList*				_solverList;		// list this solver is now on
		SimTime					_beginTime;			// initial time value
		SimTime					_endTime;			// final time value
		SimTime					_currentTime;		// time of the current step
		SimTime					_timeStep;			// current time step interval

		// Set the current time
		virtual void			currentTime(SimTime t) { _currentTime = t; }
	};

	// ----------------------------------------------------------------
	// CLASS:	ODESolver
	// EXTENDS:	Solver
	// DESC:	Abstract class for ODE solvers used to
	//			numerically evaluate the differential
	//			equations making up a dynamic model.
	// RESP:
	//		1.	Notify components of simulation state.
	//		2.	Provide utility functions for subclasses.
	// ----------------------------------------------------------------

	class ODESolver : public Solver {

	public:

		// Constructors and destructor
		ODESolver();
		virtual ~ODESolver();

		// Accessors
		inline Model*			model() { return _model; }
		virtual void			model(Model* m);

		inline SimTime			maxTimeStep() { return _maxTimeStep; }
		virtual void			maxTimeStep(SimTime t);

		inline SimTime			minTimeStep() { return _minTimeStep; }
		virtual void			minTimeStep(SimTime t);

		inline Number			errTol() { return _errTol; }
		virtual void			errTol(Number err) { _errTol = err; }

		inline int				stateVectorSize() { return model()->stateVectorSize(); }
		inline int				derivativeEvals() { return _derivativeEvals; }
		inline int				timeStepsDone() { return _timeStepsDone; }

		// Accessor for a safety margin used in estimating ideal step size.
		inline  Number			safetyMargin() { return _safetyMargin; }
		virtual void			safetyMargin(Number x) { _safetyMargin = x; }

		// debugPerformance turns on any performance messages to cerr.
		// Subclass will need to interpret as appropriate for that class.
		inline  bool			debugPerformance() { return _debugPerformance; }
		virtual void			debugPerformance(bool x) { _debugPerformance = x; }

		// Basic simulation actions
		virtual void			start();			// initialize
		virtual void			run();				// run the simulation
		virtual void			finish();			// clean up - finalize

		// Model interface functions
		virtual void			modelComponentsChanged();

	protected:
		Model*					_model;				// model being evaluated
		SimTime					_maxTimeStep;		// maximum allowed step size
		SimTime					_minTimeStep;		// minimum allowed step size
		Number					_errTol;			// error tolerance for solutions
		Number					_safetyMargin;		// fudge factor in estimating step sizes
		bool					_debugPerformance;	// write debug msgs to cerr

		// Performance statistics
		int						_derivativeEvals;	// number of computeDerivatives
		int						_timeStepsDone;		// number of time steps done

		// Subsets of model components to receive notifications
		ModelComponentVector	_stateVectorChangedRecipients;
		ModelComponentVector	_computeDerivativesRecipients;
		ModelComponentVector	_timeStepEndedRecipients;

		// Look through model components and sort by recipient groups
		virtual void			prepareRecipients();

		// Do one time step -- must be supplied if run is not overridden
		virtual void			processStep(SimTime maxDuration);

		// Utility notification functions for subclasses
		virtual void			notifyOnSetInitialState();
		virtual void			notifyOnSimulationStarted();
		virtual void			notifyOnStateVectorChanged();
		virtual void			notifyOnComputeDerivatives();
		virtual void			notifyOnTimeStepEnded();
		virtual void			notifyOnSimulationEnded();

		// Compute derivatives (state vector changed + compute derivatives)
		virtual void			computeDerivatives();

		// Utility vector functions for subclasses
		virtual Number*			allocateNumVector();
		virtual void			deleteNumVector(Number* vect);
		virtual void			zeroNumVector(Number* vect);
		virtual void			axpy(Number a, Number* x, Number* y);
		virtual void			scalex(Number a, Number* x);
		virtual void			cpxy(Number* x, Number* y);

		// Get norm of distance between vectors based on weighted values
		virtual Number			wdistxy(Number* w, Number* x, Number* y); // inf norm
		virtual Number			wdist2xy(Number* w, Number* x, Number* y); // L2 norm
	};

	// ----------------------------------------------------------------
	// CLASS:	ClockSolver
	// EXTENDS:	ODESolver
	// DESC:	A provide minimal solver that is limited
	//			to periodically notifying components
	//			of the passage of time.
	// RESP:
	//		1.	Notify all components of start and end.
	//		2.	Periodically notify of time step end
	//
	// NOTE:	No state vector or derivatives vector is
	//			used. All components are notified of
	//			time step end. notifyOnTimeStepEnded
	//			values are not solicited.
	// ----------------------------------------------------------------

	class ClockSolver : public ODESolver {

	public:

		// Constructors and destructor
		ClockSolver();
		virtual ~ClockSolver();

		// Basic simulation actions
		virtual void			start();			// initialize
		virtual void			finish();			// clean up - finalize

	protected:
		virtual void			processStep(SimTime maxDuration); // take the step
		virtual void			notifyOnTimeStepEnded();	// notify components of step
	};


	// ----------------------------------------------------------------
	// CLASS:	EulerMethodSolver
	// EXTENDS:	ODESolver
	// DESC:	Numerically evaluate the ODE system using
	//			Euler's method with fixed step size
	// RESP:
	//		1.	Invoke derivatives function
	//		2.	Compute new state vector
	// ----------------------------------------------------------------

	class EulerMethodSolver : public ODESolver {

	public:

		// Constructors and destructor
		EulerMethodSolver();
		virtual ~EulerMethodSolver();

	protected:

		// Process one time step
		virtual void		processStep(SimTime maxDuration);
	};

	// ----------------------------------------------------------------
	// CLASS:	RungeKuttaSolver
	// EXTENDS:	ODESolver
	// DESC:	Numerically evaluate the ODE system using
	//			a fourth order Runge-Kutta method with fixed
	//			step size
	// RESP:
	//		1.	Invoke derivatives function
	//		2.	Compute new state vector
	// ----------------------------------------------------------------

	class RungeKuttaSolver : public ODESolver {

	public:

		// Constructors and destructor
		RungeKuttaSolver();
		virtual ~RungeKuttaSolver();

	protected:

		// Process for one time step
		virtual	void			processStep(SimTime maxDuration);
	
	};

	// ----------------------------------------------------------------
	// CLASS:	AdaptiveSolver
	// EXTENDS:	ODESolver
	// DESC:	Numerically evaluate the ODE system using
	//			a method that is adaptive in time step size.
	// RESP:
	//		1.	Get derivatives via components
	//		2.	Compute new state vector
	//
	// NOTE:	This is a multi-order adaptive procedure.
	//			The order is dynamically adjusted based on error
	//			estimates at each time step ranging from order 1
	//			to order 3. Time step is adjusted as needed
	//			to meet error tolerance. The method is related
	//			to classical predictor-corrector schemes but
	//			is a single-step method similar to a multi-value
	//			method without the need to explicitly extract
	//			high order derivatives.
	// ----------------------------------------------------------------

	class AdaptiveSolver : public ODESolver {

	public:

		// Constructors and destructor
		AdaptiveSolver();
		virtual ~AdaptiveSolver();

		// Accessor for the number of stable steps in a row before
		// forcing a reassessment of step order.
		inline  int				stableStepLimit() { return _stableStepLimit; }
		virtual void			stableStepLimit(int n) { _stableStepLimit = n; }

	protected:
		int						_stableSteps;	// count of successive stable steps
		int						_stableStepLimit; // stable steps before forcing order retest 

		// Process for one time step
		virtual	void			processStep(SimTime hMax);
	};

	// ----------------------------------------------------------------
	// CLASS:	Event
	// EXTENDS:	none
	// DESC:	Abstract class for events.
	// RESP:
	//		1.	Know event time
	//		2.	Compare events by time (via eventOrder)
	//		3.	Perform an event specific actions when
	//			the event is acted on (subclass resp)
	//		4.	Track use count
	//		5.	Delete when use is zero
	//
	// NOTE:	In most cases the originator and recipient
	//			are aware of the type of the event shared
	//			between them. Customization of actions and
	//			data would be done by subclassing.
	//			Specialized priority queues can be done
	//			through classes or templates. See
	//			EventPriorityQueue typedef for an example.
	//
	//			When an event is created, the use count is 0.
	//			If the count goes back to 0, the event is deleted.
	//			This allows an event to be handled by multiple
	//			recipients if appropriate.
	//
	// ----------------------------------------------------------------

	class Event {

	public:

		// Constructors and destructor
		Event(SimTime t = InfinitePast);
		virtual ~Event();

		// Accessors
		inline  SimTime			eventTime()		{ return _eventTime; }
		virtual void			eventTime(SimTime t) { _eventTime = t; }

		// Data accessors that can be overridden as needed.
		virtual Number			eventParam()	{ return 0; }
		virtual Number*			eventParams()	{ return NULL; }

		// Do whatever is needed when the event is acted on.
		// If used, the subclass must override this function.
		virtual void			dispatch();

		// Numeric event class identifier to facilitate event polling etc.
		// Event classes 0 through 4095 are reserved for BNSF components
		static unsigned int		eventClassId()	{ return 0; }

		// Maintain use count and delete as needed
		inline  void			addUse() { _useCount++; }
		inline  void			removeUse() { if (--_useCount<=0) delete this; }

	protected:
		int						_useCount;		// count users of this event object
		SimTime					_eventTime;		// effective time of event
	};


	// ================================================================
	// Framework classes for external interface objects
	// ================================================================


	// ----------------------------------------------------------------
	// CLASS:	Probe
	// EXTENDS:	none
	// DESC:	Abstract class for external interfaces to
	//			state changes and events occurring in components.
	// RESP:
	//		1.	Declare functions for recording states and spiking
	//		2.	Monitor time of events and apply filter (via subclass)
	//		3.	Record spiking of neurons etc. (via subclass)
	// ----------------------------------------------------------------

	class Probe {

	public:

		// Constructors and destructor
		Probe();
		virtual ~Probe();

		// Accessors
		inline SimTime				minInterval() { return _minInterval; }
		virtual void				minInterval(SimTime t) { _minInterval = t; }

		inline SimTime				flushInterval() { return _flushInterval; }
		virtual void				flushInterval(SimTime t) { _flushInterval = t; }

		// Interface to trigger probe actions
		virtual void				timeStepEnded(ModelComponent* mc); // for comp
		virtual void				timeStepEnded(Model* m, int id=0); // all comp in model

		// Interface with object being probed as necessary
		// These are here so subclasses can hook into these interfaces
		virtual void				simulationStarted() {}
		virtual void				simulationEnded() {}

		// Allow the probe to be informed of specific events (optional)
		virtual void				signalEvent(unsigned int classId, ModelComponent* source) {}
		virtual void				signalEvent(Event* eventRaised) {}

		// Force output of any buffered probe data
		virtual void				flush() {}						// default no-op

		// Provide an interface for probes to note when a
		// they are added or removed from a component.
		virtual void				addedTo(ModelComponent* mc) {}
		virtual void				removedFrom(ModelComponent* mc) {}

	protected:

		SimTime						_minInterval;		// minimum inter-report interval
		SimTime						_intervalStart;		// start of last reporting interval

		SimTime						_flushInterval;		// time between output flushes
		SimTime						_flushTime;			// time when last flush was done

		// Return true if sufficient time has passed since the last report
		virtual bool				isReportable(SimTime eventTime);

		// Report on a collection of model components
		virtual void				reportOn(
			ModelComponentVector&	comps,		// components to report on
			SimTime					t,			// time of probe
			int						id)=0;		// id of component probed
	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalRecorder
	// EXTENDS:	Probe
	// DESC:	Abstract class for writing state to an
	//			external file. Writing of data is
	//			done by the subclass.
	// RESP:
	//		1.	Know file name and open mode
	//		2.	Open file on simulation started
	//		3.	Close file on simulation ended
	//
	// NOTES:	Actual writing is initiated by Probe via reportOn.
	//			Subclasses must supply the required interface.
	// ----------------------------------------------------------------

	class ExternalRecorder : public Probe {

	public:

		// Constructors and destructor
		ExternalRecorder(char* fn, char* mapfn=NULL);
		virtual ~ExternalRecorder();

		// Accessor for output file name
		inline char*				fileName() { return _fileName; }
		virtual void				fileName(char* fn);

		// Accessor for optional map file defining output columns
		inline  char*				mapFileName() { return _mapFileName; }
		virtual void				mapFileName(char* mfn);

		// Accessor header written flag. Setting this early
		// will bypass writing header data if it is not desired.
		// However, the testing the flag is under subclass control.
		inline  bool				headerWritten() { _headerWritten; }
		virtual void				headerWritten(bool x) { _headerWritten=x; }

		// Handle simulation state changes
		virtual void				simulationStarted();
		virtual void				simulationEnded();

		// Force output of any buffered probe data
		virtual void				flush();

	protected:
		char*						_fileName;
		char*						_mapFileName;
		FILE*						_file;
		FILE*						_mapFile;

		// Flag to track if headers have been written.
		// Initialized to false, but otherwise set by subclasses.
		bool						_headerWritten;

		// Utility functions for subclasses
		virtual void				openMapFile();
		virtual void				closeMapFile();

	};

	// ----------------------------------------------------------------
	// CLASS:	ExternalStateRecorder
	// EXTENDS:	ExternalRecorder
	// DESC:	Records state variables in external file
	//			using comma separated value (CSV) format.
	// RESP:
	//		1.	Write state vector values for all states
	//			in the model of the component attaching 
	//			this recorder.
	// ----------------------------------------------------------------

	class ExternalStateRecorder : public ExternalRecorder {

	public:

		// Constructors and destructor
		ExternalStateRecorder(char* fn, char* mapFileName = NULL);
		virtual ~ExternalStateRecorder();

		// Report on a collection of components
		virtual void				reportOn(ModelComponentVector& comps, SimTime t, int id);

	protected:
		virtual void				writeColumnHeader(ModelComponentVector& comps);
		virtual void				writeState(ModelComponentVector& comps, SimTime t, int id);
	};


}; // end of namespace

#endif // #ifndef __BSNF_H_
