// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_sim.cpp
//
// Release:		0.0.1
// Author:		John Baker
// Updated:		13 Sep 2003
//
// Description:
//
// This file provides the body of BNSF classes.
// The sequence is generally the same as the definitions
// in the associated header file, but a good PDE is
// almost a necessity for any reasonable development effort.
//
// This file implements classes used in all simulations.



#include "bnsf_sim.h"
#include <algorithm>

using namespace std;
using namespace BNSF;


// ====================================================================
// Global functions
// ====================================================================


// Order function used in ordering event queues 
bool BNSF::eventOrder(Event* rhs, Event* lhs) {
	return lhs->eventTime()<rhs->eventTime();
}

// Order function used in ordering Solver list headers for running
bool BNSF::solverListRunOrder(SolverList* lhs, SolverList* rhs)
{
	// Determine what time the current step would end
	const SimTime			lhsEndTime = lhs->timeStep()+lhs->evalTime();
	const SimTime			rhsEndTime = rhs->timeStep()+rhs->evalTime();

	// First, order by the end times
	if (lhsEndTime > rhsEndTime ) return true;
	if (lhsEndTime < rhsEndTime ) return false;

	// Assuming that end times match, order by time step, smallest first
	return lhs->timeStep() >= rhs->timeStep();
}



// ====================================================================
// Model class body
// ====================================================================



Model::Model()
{
	// Size of state vector allocated (will be >0 when allocated)
	_stateVectorSize = 0;

	// Set null pointers as a stand-in
	_stateVector = NULL;
	_derivVector = NULL;
	_weightVector = NULL;
	_solver = NULL;

	// To save space, the probe vector is allocated only when needed
	_probesPtr = NULL;

	// Use the default random number stream to start with.
	_uniformRandom = NULL;
}

Model::~Model() 
{
	ModelComponentVector	temp;
	ModelComponentVectorIt	it;

	// Unhook with any current solver
	solver(NULL);

	// Tell components the model is gone.
	// Use a temp vector to avoid removing
	// components underneath the iterator..
	swap(_components,temp);
	for (it=temp.begin();it!=temp.end();it++) {
		(*it)->model(NULL);
	}

	// Delete any state vectors held by the model
	delete[] _weightVector;
	delete[] _derivVector;
	delete[] _stateVector;

	// Delete any probes remaining
	delete _probesPtr;
}

// Set the ODE solver for this model
void Model::solver(ODESolver* newSolver)
{
	ODESolver*		oldSolver = solver();

	if (oldSolver!=newSolver) {

		// Break any old relationship
		if (oldSolver != NULL && oldSolver->model()==this) {
			_solver = NULL;
			oldSolver->model(NULL);
		}

		// Hook up new relationship
		_solver = newSolver;
		if (_solver != NULL) {
			_solver->model(this);
		}
	}
}

// Add a new component to the model
void Model::addComponent(ModelComponent* comp)
{
	// Add the new component to the list
	_components.push_back(comp);

	// Tell the solver, if any
	if (solver() != NULL) {
		solver()->modelComponentsChanged();
	}
}

// Remove a component from the model.
// Note that state variables would need to ultimately be
// reallocated if any removes are done.
void Model::removeComponent(ModelComponent* comp)
{
	ModelComponentVectorIt		last;

	// Remove from the components
	last=remove(_components.begin(),_components.end(),comp);
	_components.resize(last - _components.begin());

	// Tell the solver, if any
	if (solver() != NULL ) {
		solver()->modelComponentsChanged();
	}
}

// Remove all components. This prevents an N^2 operation
// when deleting components.
void Model::clearComponents()
{
	// Empty the components vector
	_components.clear();
	_components.resize(0);

		// Tell the solver, if any
	if (solver() != NULL ) {
		solver()->modelComponentsChanged();
	}
}

// Allocate a vector for probes using lazy initialization
ProbeTargetVector& Model::probes()
{
	if (_probesPtr==NULL) {
		_probesPtr = new ProbeTargetVector;
	}
	return *_probesPtr;
}

// Associate a probe object with this component
void Model::addProbe(Probe* pr, ModelComponent* mc)
{
	probes().push_back(ProbeTargetPair(pr,mc));
}

// Remove a probe associated with this component
void Model::removeProbe(Probe* pr, ModelComponent* mc)
{
	ProbeTargetPair			ptp(pr,mc);
	ProbeTargetVectorIt		last;
	int						newSize;

	// Remove a matching probe pair from the vector and
	// if there are no more entries, delete the vector itself.

	last=remove(probes().begin(),probes().end(),ptp);
	newSize = last-probes().begin();
	
	if (newSize>0) {
		probes().resize(last - probes().begin());
	}
	else {
		delete _probesPtr;
		_probesPtr = NULL;
	}
}

// Do first time setup before interacting with solver.
void Model::simulationStarted()
{
	// Initialize times
	if (solver()!=NULL) {
		_currentTime = solver()->currentTime();
		_timeStepStart = _currentTime;
	}
	else {
		FatalError("(Model::simulationStarted) Solver not yet allocated");
	}
	_timeStepSize = 0;
	_stepSizeChanged = true;

	// Set up state vector if needed
	if (stateVectorSize() == 0 ) {

		// Put components in order by dependency
		orderComponents();

		// Allocate storage for all state related vectors
		allocateVectors();

		// Assign state vector offsets and tell components
		assignSVOffsets();
	}

	// Do probes notification only if there are probes
	if (_probesPtr!=NULL) {

		ProbeTargetVectorIt pit;

		for (pit=probes().begin();pit!=probes().end();pit++) {
			pit->first->simulationStarted();
		}
	}	
}

// Prepare for ending of the time step by updating
// various internal timestamps.
void Model::prepareTimeStepEnded(SimTime timeNow)
{
	
	// Save the current time
	currentTime(timeNow);

	// Compute the size of the time step just
	// passed and set a flag indicating that
	// the time step changed. Allow a little 
	// tolerance for differences to allow for 
	// rounding errors. A typical use by
	// components is to update a cache of values
	// derived from time step size. When new size
	// is zero, a new step size is indicated to
	// allow components to initialize caches
	// dependent on step size at start-up.

	SimTime newStepSize = currentTime()-timeStepStart();

	if (newStepSize==0 || fabs(newStepSize-timeStepSize())>EpsilonTime) {
		_timeStepSize = newStepSize;
		_stepSizeChanged = true;
	}
}

// Handle end of time step by invoking any probes.
// Also roll over the current time
void Model::timeStepEnded()
{
	// Do notification only if there are probes
	if (_probesPtr!=NULL) {

		ProbeTargetVectorIt pit;

		for (pit=probes().begin();pit!=probes().end();pit++) {
			pit->first->timeStepEnded(pit->second);
		}
	}

	// Roll over the step start time
	_timeStepStart = currentTime();
	_stepSizeChanged = false;
}

// Allocate storage for the state vector, derivative vector,
// and weight vector;
void Model::allocateVectors()
{
	int k;

	// Get the size to be allocated.
	// Note that the first entry is reserved.
	_stateVectorSize = 1;
	for (k=0;k<_components.size();k++) {
		_stateVectorSize += _components[k]->numStateVar();
	}

	// Allocate the storage needed.
	_stateVector = new Number[_stateVectorSize];
	_derivVector= new Number[_stateVectorSize];
	_weightVector= new Number[_stateVectorSize];

	// Clear value to zero initially
	for (k=0;k<_stateVectorSize;k++) {
		_stateVector[k]=0;
		_derivVector[k]=0;
		_weightVector[k]=0;
	}
}

// Compute state vector offsets for current components.
// This should be done before components set their initial
// values but cannot be done when components are added due
// to the limitation on accessing virtual functions during
// construction of an object, which is typically the stage
// when components are associated with a model.
void Model::assignSVOffsets()
{
	typedef pair<int,int>		SortEntType;
	typedef vector<SortEntType>	SortEntVector;

	int							nsv,i,k;
	ModelComponent*				mc;

	SortEntType					sortEnt;
	SortEntVector				sortVect;

	// Sort the components by domain and compartment number.
	// This ordering allows solvers to group items by domain
	// without changing the order in the compartments vector.
	for (k=0;k<_components.size();k++) {
		sortEnt.first = _components[k]->domain();
		sortEnt.second = k;
		sortVect.push_back(sortEnt);
	}
	sort(sortVect.begin(),sortVect.end());

	// Assign state vector numbers in sort order but
	// leaving the first slot empty (simplifies debugging)

	_stateVectorSize = 1;
	for (k=0;k<sortVect.size();k++) {
		i = sortVect[k].second;
		mc = _components[i];
		nsv = mc->numStateVar();
		if (nsv>0) {
			mc->svOffset(_stateVectorSize);
			_stateVectorSize += nsv;
		}
		else {
			mc->svOffset(0); // entry 0 is reserved dummy
		}
	}
}

// Topologically sort components based on dependency relationships.
// The resulting order is reflected in the components array and
// should be used by solvers when invoking component functions.
// See Knuth, Art of Computer Programming Vol I for a description
// of the algorithm used.

void Model::orderComponents()
{

	// Table entry structure for holding component dependencies
	typedef struct ts_s {
		int count;				// Count of prerequisites remaining
		int next;				// Linked list of zero count entries
		vector<int> successors;	// Immediate successors
	} ts_t;

	const int				sz=_components.size();	// Size of components list
	int						n,m,k,s;				// Work variables
	int						nextq = -1;				// Header for list of ready entries
	bool					sortNeeded = false;		// Is a sort needed or not

	ts_t*					ts = new ts_t[sz];		// Topological sort table
	ModelComponentVector	mcv;					// Copy of components vector

	typedef map<ModelComponent*,int>	mcmap_t;
	typedef mcmap_t::iterator			mcmapIt;	// Iterator for the map
	typedef mcmap_t::value_type			mcmapVT;	// Value type inserted into map
	mcmap_t								mcmap;		// Map relating address to index

	// Initialize the topological sort table and map
	for (n=0;n<sz;n++) {
		ts[n].count=0;
		ts[n].next=-1;
		mcmap.insert(mcmapVT(_components[n],n));
	}

	// Populate the sort table to represent immediate successors.
	for (n=0;n<sz;n++) {
		if (_components[n]->hasPrerequisites() ) {

			ModelComponentVector	prereqs(_components[n]->prerequisites());
			ModelComponentVectorIt	pit;
			mcmapIt					mit;

			ModelComponent*			mcdebug = _components[n];

			ts[n].count = prereqs.size();	// save count of prereqs to be worked on

			for (pit=prereqs.begin();pit!=prereqs.end();pit++) {

				// Look up the index of the prereq component in mcmap
				mit=mcmap.find(*pit);
				if (mit==mcmap.end()) {
					FatalError("(Model::orderComponents) Invalid prerequisite returned.");
				}
				m=mit->second;

				// Invert the prereq relationship and save as a successor.
				ts[m].successors.push_back(n);
				if (n<m) {
					sortNeeded = true;
				}
			}
		}
	}				

	// Only do the sort if there are out of order entries
	if (sortNeeded) {

		// Put the current components in mcv leaving _components empty.
		swap(_components,mcv);


		// Build the initial list of entries with a zero count
		// while trying to preserve any original ordering.
		for (n=sz-1;n>=0;n--) {
			if (ts[n].count==0) {
				ts[n].next = nextq;
				nextq = n;
			}
		}

		// Take entries with a zero count off of the list and put them
		// into _components, adjusting counts as we go.
		while (nextq != -1) {
			
			// Take this entry out of the ready queue
			n=nextq;
			nextq=ts[n].next;

			// Put it into _components
			_components.push_back(mcv[n]);
			
			// Decrement the counts of successors and if those counts
			// go to zero put the associated entries in the ready queue.
			// Keep the ready queue sorted to preserve original order.
			for (k=ts[n].successors.size()-1;k>=0;k--) {
				s=ts[n].successors[k];
				if (--ts[s].count==0) {
					int curq = nextq;
					int prevq = -1;
					while (curq!=-1 && curq<s) {
						prevq = curq;
						curq=ts[curq].next;
					}
					if (prevq==-1) {
						ts[s].next = nextq;
						nextq=s;
					}
					else {
						ts[s].next=curq;
						ts[prevq].next=s;
					}
				}
			}
		}
	}

	// Dispose of allocated table
	delete[] ts;

	// Make sure every component was copied.
	// If not, there was a dependency loop that cannot be fixed here.
	if (_components.size() != sz) {
		FatalError("(Model::sortComponents) Loop found in component prerequisites");
	}
}



// ====================================================================
// ModelEntity class body
// ====================================================================



// Construct a new instance
ModelEntity::ModelEntity() {}

// Destroy this instance
ModelEntity::~ModelEntity() {}

// Perform any unique actions when the event is acted on
// If used, this must be overridden in subclasses.
void ModelEntity::dispatchEvent(Event* ev)
{
	FatalError("(ModelEntity::dispatchEvent) No event-specific dispatch function");
}



// ====================================================================
// ModelComponent class body
// ====================================================================


// Constructors and destructor
ModelComponent::ModelComponent(Model* m)
{
	_svOffset = 0;
	_svArray = NULL;
	_model = NULL;
	model(m);
}

ModelComponent::~ModelComponent() 
{
	// Notify subscribers, if any, that this object is terminated.
	if (subscribers()!=NULL) {
		changed(terminatedChange);
	}

	// Unhook from model, if any.
	if (model()!=NULL && joinModelAsComponent() ) {
		model()->removeComponent(this);
	}
}

// Set the model and become a component of it
void ModelComponent::model(Model* m)
{
	// If this object does not actually join
	// the model, just remember the association
	// but do not inform model.
	if (!joinModelAsComponent() ) {
		_model = m;
		return;
	}

	// Otherwise, if there is no change, stop now
	if (_model==m) 
		return;

	// If there is an existing model, unhook from it
	if (_model!=NULL) {
		_model->removeComponent(this);
	}

	// Add into the the new models components
	_model=m;
	if (m!=NULL) {
		m->addComponent(this);
	}
}

// Add a subscriber to this component
void ModelComponent::addSubscriber(ModelComponent* mc)
{
	ModelComponentVector*	subs = subscribers();

	if (subs==NULL) {
		FatalError("(ModelComponent::addSubscriber) Subclass must supply subscriber vector.");
		return;
	}
	subs->push_back(mc);
}

// Remove a subscriber from this component
void ModelComponent::removeSubscriber(ModelComponent* mc)
{
	ModelComponentVector*	subs = subscribers();
	ModelComponentVectorIt	last;

	if (subs==NULL) {
		FatalError("(ModelComponent::removeSubscriber) Subclass must supply subscriber vector.");
		return;
	}

	// Remove the component from the subscriber vector
	last=remove(subs->begin(),subs->end(),mc);
	subs->resize(last - subs->begin());
}

// Signal that a change has occurred.
void ModelComponent::changed(int reason)
{
	ModelComponentVector*	subs = subscribers();
	ModelComponentVectorIt	it;

	// Quietly ignore if no subscribers are defined.
	if (subs!=NULL) {

		// Notify subscribers
		for (it=subs->begin();it!=subs->end();it++) {
			(*it)->updateFrom(this,reason);
		}
	}
}

// Return a vector of prerequisites consisting of zero or one prereqs.
ModelComponentVector ModelComponent::prerequisites()
{
	ModelComponentVector prereqs;

	ModelComponent* prereq = prerequisite();

	if (prereq!=NULL) {
		prereqs.push_back(prereq);
	}
	
	return prereqs;
}

// Add a probe to this component
void ModelComponent::addProbe(Probe* pr)
{
	// Probes added this way are the exception.
	// The model keeps track of all such probes
	// to reduce storage overhead.
	model()->addProbe(pr, this);

	// Tell probe where it was added
	pr->addedTo(this);
}

// Remove a probe from this component
void ModelComponent::removeProbe(Probe* pr)
{
	// Remove probe of this component 
	// the from model's collection.
	model()->removeProbe(pr, this);

	// Let probe know it was removed.
	pr->removedFrom(this);
}

// Associate this component with a controller.
// A model and solver must already be established.
void ModelComponent::addToController(Controller* cont)
{
	cont->addSolver( solver() );
}


// set default values in the weight vector used for tolerances
void ModelComponent::setWeightValues()
{
	// Set each weight value corresponding with a state
	// variable to 1 as a default.
	for (int k=0;k<numStateVar();k++) {
		weightValue(k) = 1;
	}
}

// This must be overridden by subclass if invoked by ODE solver
void ModelComponent::computeDerivatives()
{
	if (numStateVar() > 0) {
		FatalError("(ModelComponent::computeDerivatives) Subclass must override if used");
	}
}

// This must be overridden by subclass if invoked by ODE solver
void ModelComponent::localStateUpdate(SimTime h, CNStepType stepType)
{
	if (numStateVar() > 0) {
		FatalError("(ModelComponent::localStateUpdate) Subclass must override if used");
	}
}

// This must be overridden by subclass if invoked by ODE solver
void ModelComponent::jacobianIndex(int n)
{
	FatalError("(ModelComponent::jacobianIndex) Subclass must override if used");
}

// Provide a default component name
const char* ModelComponent::componentName()
{ 
	return "Other"; 
}

// Provide a default array of state labels
const char** ModelComponent::stateLabels() 
{
	// Provide a few default values. The actual number used is determined by numAllStateVar.
	static const char* labels[] = {"S1", "S2",  "S3",  "S4",  "S5",  "S6",  "S7",  "S8",
							 "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16"};
	return labels;
}


// Provide a default set of units of measure
Number* ModelComponent::unitsOfMeasure()
{
	// Provide a few default values. The actual number used is determined by numAllStateVar.
	static Number units[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	return &(units[0]);
}

// Add to a vector of the components to be included when
// probing the contents of this component and related items.
// By default only this component is added in the vector.
void ModelComponent::addToComponentsToProbe(ModelComponentVector& comps)
{
	comps.push_back(this);
}

// Return an internal state value. Subclass must
// override to provide the value if this is used.
Number ModelComponent::internalStateValue(int n)
{
	FatalError("(ModelComponent::internalStateValue) Subclass must override this function.");
	return 0;	// keep compiler happy
}

// Return a unit of measure converted state value
Number ModelComponent::stateValueConverted(int n) 
{
	if ( n<0 || n>numStateVar()+numInternalStateVar() ) {
		FatalError("(ModelComponent:stateValueConverted) index out of bounds");
	}

	Number	units = unitsOfMeasure()[n];
	Number	value;
	int		nint = numInternalStateVar();

	// Internal state values come before ODE state values.
	// Adjust n accordingly to get the state value.

	value = n<nint ? internalStateValue(n) : stateValue(n-nint);

	return value/units;
}



// ====================================================================
// Controller class body
// ====================================================================



// Initialize a new instance
Controller::Controller()
{
	_beginTime = 0;
	_evalTime = 0;
	_endTime = InfiniteFuture;
}

// Destroy this instance.
Controller::~Controller() 
{
	SolverListMapIt		it;

	// Clean up the ODESolverList instances in case finish not done
	for (it= _solvers.begin(); it != _solvers.end(); it++) {
		delete it->second;
	}
}

// Add a new ODESolver to this Controller
void Controller::addSolver(Solver* s)
{
	SolverListMapIt			it;
	SolverList*				newHeader;
	SimTime					h = s->timeStep();
	SimTime					t;

	// Tell the solver who we are
	s->controller(this);

	// Get a list header based on time step
	it = _solvers.find(h);
	if (it == _solvers.end() ) {

		// Create a new header just for this solver
		// and put it into the solver map and queue.
		newHeader = new SolverList;
		newHeader->addLast(s);

		// Compute a starting eval time as a multiple of
		// the time step from the starting point.
		t = h*floor((s->currentTime()-beginTime())/h) + beginTime();
		newHeader->evalTime( t );

		// Put the new header in the map of all headers and into
		// the runnable queue.
		_solvers.insert(SolverListMapValueType(h,newHeader));
		_runnableQueue.push(newHeader);
		newHeader->isRunnable(true);

		// DEBUG CODE
		// cerr<<"Added list "<<newHeader
		//	<<" step "<<newHeader->timeStep()
		//	<<" end time "<<newHeader->timeStep()+newHeader->evalTime()
		//	<<endl;
	}
	else {
		// Otherwise, add to the header found
		newHeader = it->second;
		if (newHeader->isEmpty() ) {

			// Reset the effective time to now truncated to a multiple
			// of the time step. This improves regularity but is mostly
			// an exercise in neatness so that cases where things get off
			// the track are more apparent.
			t = h*floor((s->currentTime()-beginTime())/h) + beginTime();
			newHeader->evalTime( t );

			// Put on runnable queue now that header will not be empty
			if (!newHeader->isRunnable() ) {
				_runnableQueue.push(newHeader);
				newHeader->isRunnable(true);

				// DEBUG CODE
				// cerr<<"Reinstated list "<<newHeader
				//	<<" step "<<newHeader->timeStep()
				//	<<" end time "<<newHeader->timeStep()+newHeader->evalTime()
				//	<<endl;
			}
		}
		newHeader->addLast(s);
	}
}

// Add a model based on its current solver
void Controller::addModel(Model* m) 
{
	addSolver(m->solver() );
}

// Add a model component based on its model
void Controller::addModelComponent(ModelComponent* mc)
{
	addModel(mc->model() );
}

// Change the header associated with a Solver
// because of a change in time step
void Controller::changeSolverTimeStep(Solver* s, SimTime oldTS)
{
	SolverListMapIt		it;

	// Find the old header, which must already exist
	it = _solvers.find(oldTS);
	if ( it == _solvers.end() ) {
		FatalError("(Controller::changeSolverTimeStep) "
			"Solver does not have matching list header.");
		return;
	}

	// Remove from the old list
	it->second->removeNode(s);

	// Add under a new header
	addSolver(s);
}

// Prepare to run the simulation
void Controller::start()
{
	SolverListMapIt		it;			// iterator which derefs to a pair
	SolverList*			header;
	Solver*				node;

	// Go through all known solvers and tell each of their
	// list members we are starting now.
	for (it=_solvers.begin(); it != _solvers.end(); it++) {
		header = it->second;
		header->evalTime(beginTime() );
		node = header->firstNode();
		while (node!=NULL) {
			node->beginTime( beginTime() );
			node->start();
			node = node->nextNode();
		}
	}
}

// Run the simulation up to the ending time.
void Controller::run()
{
	SolverList*		nextToRun;
	SimTime			solverEndTime;

	for(;;) {

		// Make sure there is something to run
		if (_runnableQueue.empty() ) {

			// this really should not happen - time goes on
			FatalError("(Controller::run) Premature end because nothing is runnable.");
			return;
		}

		// Locate the next to run in sequence
		nextToRun = _runnableQueue.top();
		_runnableQueue.pop();

		// DEBUG CODE
		// cerr<<"Running list "<<nextToRun
		//	<<" step "<<nextToRun->timeStep()
		//	<<" end time "<<nextToRun->timeStep()+nextToRun->evalTime()
		//	<<endl;

		// Skip over any lists that have gone empty
		if (nextToRun->isEmpty() ) {
			nextToRun->isRunnable(false);
			continue;
		}

		// Stop when running would put us over the end time, but
		// allow for case where there has been a slight roundoff of times.
		if (nextToRun->evalTime()+EpsilonTime>endTime() ) {
			_runnableQueue.push(nextToRun);
			break;
		}

		// Set the end time for this solver list but do not
		// go over the current end time for the evaluation.
		solverEndTime = nextToRun->evalTime()+nextToRun->timeStep();
		if (solverEndTime>endTime() ) {
			solverEndTime = endTime();
		}

		// Run all the solvers in this group and update time
		nextToRun->runAllUpTo(solverEndTime);

		// Check to see if we have progressed in time
		if (nextToRun->evalTime() > _evalTime ) {
			_evalTime = nextToRun->evalTime();
		}

		// Queue up for running the next step
		_runnableQueue.push(nextToRun);
	}
}

// Run up to an end time, preserving the old end
void Controller::runUpTo(SimTime t)
{
	// Save the old end time
	SimTime saveEndTime = endTime();

	// Run to the new time point
	endTime(t);
	run();

	// Restore the old end time
	endTime(saveEndTime);
}

// Run for a specified duration from where we last left off
void Controller::runForDuration(SimTime dt)
{
	// Save the old end time
	SimTime saveEndTime = endTime();

	// Run for the duration specified
	endTime(dt+evalTime() );
	run();

	// Restore the old end time
	endTime(saveEndTime);
}

// Wrap-up at the end of the simulation
void Controller::finish()
{
	SolverListMapIt		it;		// iterator which derefs to a pair
	Solver*				node;

	// Go through all known solvers and tell each of their
	// list members we are finishing now.
	for (it=_solvers.begin(); it != _solvers.end(); it++) {
		node = it->second->firstNode();
		while (node!=NULL) {
			node->finish();
			node = node->nextNode();
		}
	}

	// Discard references in the solver map and queue
	SolverListMap				emptyMap;
	SolverListPriorityQueue		emptyQueue;

	swap(_solvers,emptyMap);
	swap(_runnableQueue,emptyQueue);
}



// ====================================================================
// SolverList class body
// ====================================================================



// Create a new list header
SolverList::SolverList()
{
	// create list as empty
	_firstNode = NULL;
	_lastNode = NULL;
	_evalTime = InfiniteFuture;
	_timeStep = 0;
	_isRunnable = false;
}

// Destroy the list header but leave nodes alone
SolverList::~SolverList() {}

// Add a node at the beginning of the list
void SolverList::addFirst(Solver* node)
{
	// Make sure node is not already claimed
	if (node->solverList()!=NULL) {
		FatalError("(SolverList::addFirst) Node must be removed from old list first");
	}

	// Link in the new node
	if (isEmpty() ) {
		_lastNode = node;
		evalTime( node->currentTime() );
		timeStep( node->timeStep() );
	}
	else {
		node->insertBefore(_firstNode);
	}
	_firstNode = node;

	// Tell the node what list it is on
	node->solverList(this);
}

// Add a node at the end of the list
void SolverList::addLast(Solver* node)
{
	// Make sure node is not already claimed
	if (node->solverList()!=NULL) {
		FatalError("(SolverList::addLast) Node must be removed from old list first");
	}

	// Link in the new node
	if (isEmpty() ) {
		_firstNode = node;
		timeStep( node->timeStep() );
	}
	else {
		node->insertAfter(_lastNode);
	}
	_lastNode = node;

	// Tell the node what list it is on
	node->solverList(this);
}

// Remove a node from this list
void SolverList::removeNode(Solver* node)
{
	// Can't remove from an empty list -- always an error
	if (isEmpty() ) {
		FatalError("(SolverList::removeNode) Cannot remove from empty list.");
		return;
	}

	// Unhook if this is the first or last node
	if (_firstNode == node) {
		_firstNode = node->nextNode();
	}
	if (_lastNode == node) {
		_lastNode = node->prevNode();
	}

	// Take the node out of the doubly linked list
	node->removeFromList();
}

// Run all the ODESolvers in the list and adjust evaluation
// when they are all done to reflect the change
void SolverList::runAllUpTo(SimTime endTime)
{
	Solver*			node = firstNode();
	Solver*			next = NULL;

	while (node != NULL) {

		// Get the next node first in case the node
		// is removed because it changes its time step.
		next = node->nextNode();

		// Run the current node
		node->runUpTo(endTime);

		// DEBUG CODE
		// if (node->currentTime()>endTime) {
		//	cerr<<"List "<<this
		//		<<" solver "<<node
		//		<<" end time "<<endTime
		//		<<" current time "<<node->currentTime()
		//		<<endl;
		// }

		// Move on to the next node
		node = next;
	}
	evalTime(endTime);
}



// ====================================================================
// Solver class body
// ====================================================================



// Initialize a new instance
Solver::Solver() 
{
	using namespace UOM;

	_nextNode = NULL;				// no next node
	_prevNode = NULL;				// no previous node
	_controller = NULL;				// no controller
	_solverList = NULL;				// no solver list
	_beginTime = 0;					// by default t0 = 0
	_currentTime = 0;				// current = begin at start
	_endTime = InfinitePast;		// no end time set
	_timeStep = 1*msec;				// default time step (arbitrary)
}

// Destroy this instance
Solver::~Solver() 
{
	// Disconnect from any solver list
	if (_solverList!=NULL) {
		_solverList->removeNode(this);
	}
}


// Insert before the node provided in a list
void Solver::insertBefore(Solver* node)
{
	// Do not link to ourselves
	if (node==this) {
		FatalError("(Solver::insertBefore) Can't add entry before itself");
	}
	_nextNode = node;
	_prevNode = node->_prevNode;
	node->_prevNode = this;
	if (_prevNode != NULL) {
		_prevNode->_nextNode = this;
	}
}

// Insert after the node provided in a list
void Solver::insertAfter(Solver* node)
{
	// Do not link to ourselves
	if (node==this) {
		FatalError("(Solver::insertAfter) Can't add entry after itself");
	}

	// Link in the new node
	_prevNode = node;
	_nextNode = node->_nextNode;
	node->_nextNode = this;
	if (_nextNode != NULL) {
		_nextNode->_prevNode = this;
	}
}

// Remove this node from the list
void Solver::removeFromList()
{
	if (_prevNode != NULL) {
		_prevNode->_nextNode = _nextNode;
	}
	if (_nextNode != NULL) {
		_nextNode->_prevNode = _prevNode;
	}
	_prevNode = NULL;
	_nextNode = NULL;
	_solverList = NULL;
}

// Set the time step and inform controller if necessary
void Solver::timeStep(SimTime h)
{
	SimTime oldTimeStep;

	// If the time step was changed, inform the controller if any
	if (h != timeStep() ) {
		if (controller() != NULL) {
			oldTimeStep = timeStep();
			_timeStep = h;
			controller()->changeSolverTimeStep(this, oldTimeStep);
		}
		else {
			_timeStep = h;
		}
	}
}

// Set a time step less than or equal to the value provided
SimTime Solver::timeStepLessOrEqual(SimTime h)
{
	int			high,low,mid;
	double		hSel;

	// Table of mantissa's used with successive ratio of 2^.125
	const double mTbl[8] = {
		0.50000000000000,		// [0]
		0.54525386633263,		// [1]
		0.59460355750136,		// [2]
		0.64841977732550,		// [3]
		0.70710678118655,		// [4]
		0.77110541270397,		// [5]
		0.84089641525371,		// [6]
		0.91700404320467 };		// [7]

	// Get exponent and mantissa of time step so that mantissa
	// can be rounded to one of a limited set of values.
	int			hExp;
	double		hMant = frexp(h,&hExp);

	// Make sure a zero value was not provided for h (just in case)
	if (hMant==0) {
		FatalError("(Solver::timeStepLessOrEqual) step size must not be zero");
		return 0;
	}

	// Binary search of the mantissa table
	// This makes sense mostly if there are
	// more entries than 8 (just in case).
	low = 0;
	high = 8;
	while (high-low>1) {
		mid = (high+low)/2;
		if (mTbl[mid]>hMant) {
			high=mid;
		}
		else {
			low=mid;
		}
	}

	// Put the mantissa and exponent back together
	hSel = ldexp(mTbl[low],hExp);

	// Set the new time step and return
	timeStep(hSel);
	return hSel;
}

// Run up to a new end time
void Solver::runUpTo(SimTime t)
{
	endTime(t);
	run();
}



// ====================================================================
// ODESolver class body
// ====================================================================



// Initialize a new instance
ODESolver::ODESolver() 
{
	_model = NULL;					// initialize pointer

	_maxTimeStep = InfiniteFuture;	// default max time step (none)
	_minTimeStep = 1*UOM::nanosec;	// default min time step (very very small)
	_errTol = Number( 1e-3 );		// default error tolerance
	_safetyMargin = 0.75;			// default safety margin

	_derivativeEvals = 0;			// zero performance counter
	_timeStepsDone = 0;				// zero performance counter

	_debugPerformance = false;		// reset debug flag
}

// Destroy this instance
ODESolver::~ODESolver() 
{
	// We were once in charge, but now we are dead.
	if (model() != NULL) {
		model()->solver(NULL);
	}
}

// Set up relationship with the model
void ODESolver::model(Model* newModel)
{
	Model*		oldModel = model();

	if (oldModel!=newModel) {

		// Break any old relationship
		if (oldModel != NULL && oldModel->solver()==this) {
			_model = NULL;
			oldModel->solver(NULL);
		}

		// Hook up the new relationship
		_model = newModel;
		if (_model!=NULL) {
			_model->solver(this);
		}
		
	}
}

// Set time step maximum
void ODESolver::maxTimeStep(SimTime h)
{
	_maxTimeStep = h;

	// Adjust current time step down if needed
	if (timeStep() > h) {
		timeStep(h);
	}
}

// Set time step minimum
void ODESolver::minTimeStep(SimTime h)
{
	_minTimeStep = h;

	// Adjust current time step up if needed
	if (timeStep() < h) {
		timeStep(h);
	}
}

// Do preliminary initializations (once only)
void ODESolver::start()
{
	// Set the starting time (once)
	currentTime( beginTime() );

	// Tell model we are starting
	model()->simulationStarted();

	// Build lists of who gets what message
	prepareRecipients();

	// Tell components we are starting
	notifyOnSimulationStarted();

	// Have components set initial conditions
	notifyOnSetInitialState();

	// State vector is now set first time
	notifyOnStateVectorChanged();

	// Compute the initial derivatives
	computeDerivatives();

	// Inform components that initialization is done
	notifyOnTimeStepEnded();

}

// Respond to changes in the model components
void ODESolver::modelComponentsChanged()
{
	// A subclass would need to respond to changes 
	// in the components after things are initialized/
	if (model()->stateVector()!=NULL && hasController() ) {
		FatalError("(ODESolver::modelComponentChanged) "
			"Feature not supported by this solver.");
	}
}

// Run the simulation up until the end of the simulation interval.
void ODESolver::run()
{
	// Process one step at a time until the end of the evaluation
	// interval is reached.
	while( currentTime()+minTimeStep()/2 < endTime() ) { 

		// Do one step at a time up to the current time remaining
		processStep( endTime() - currentTime() );
	}
}

// This must be overridden by subclasses (if invoked) -- i.e. an almost pure function.
void ODESolver::processStep(SimTime maxDuration)
{
	FatalError("(ODESolver::processStep) Subclass must override this function.");
}

// Do final clean up
void ODESolver::finish() 
{
	// Tell components that we are done
	notifyOnSimulationEnded();

	// We are done with this controller
	controller(NULL);

	// Disconnect from model
	model()->solver(NULL);
}


// Build vector of recipients of each type of notification by
// quering each component to see what notifications it should receive.
void ODESolver::prepareRecipients()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = model()->components();

	for (it=comps.begin();it!=comps.end();it++) {
		if ( (*it)->notifyOnStateVectorChanged() ) {
			_stateVectorChangedRecipients.push_back(*it);
		}
		if ( (*it)->updatesDerivVector() ) {
			_computeDerivativesRecipients.push_back(*it);
		}
		if ( (*it)->notifyOnTimeStepEnded() ) {
			_timeStepEndedRecipients.push_back(*it);
		}
	}
}

// Ask components to set initial conditions and weights.
// Note that state vectors must already be in place.
void ODESolver::notifyOnSetInitialState()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = model()->components();

	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->setInitialState();
		(*it)->setWeightValues();
	}
}

// Notify recipients of simulation started condition
void ODESolver::notifyOnSimulationStarted()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = model()->components();

	// Tell model what time it is
	model()->currentTime( currentTime() );
	model()->simulationStarted();

	// Notify components to do initializations
	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->simulationStarted();
	}
}

// Notify all recipients of a (possibly) new state vector with new values.
void ODESolver::notifyOnStateVectorChanged()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = _stateVectorChangedRecipients;

	// Tell the model about the current time
	model()->currentTime(currentTime());

	// Tell any interested components that the state vector was updated
	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->stateVectorChanged();
	}
}

// Notify recipients that they should compute derivatives
void ODESolver::notifyOnComputeDerivatives()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = _computeDerivativesRecipients;

	// tell the model about the current time
	model()->currentTime(currentTime());

	// ask each component to compute its state derivatives
	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->computeDerivatives();
	}
}

// Notify recipients that the time step has ended.
// Any events should be raised at this point.
// Data unique to prior step can be released and any
// necessary setup for the next step done.
void ODESolver::notifyOnTimeStepEnded()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = _timeStepEndedRecipients;

	// Prepare model for ending the step
	model()->prepareTimeStepEnded( currentTime() );

	// Inform each component that the step has ended
	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->timeStepEnded();
	}

	// Inform model that time step has ended
	model()->timeStepEnded();

	// Increment the performance counter
	_timeStepsDone++;
}

// Notify recipients that the simulation has ended
void ODESolver::notifyOnSimulationEnded()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = model()->components();

	// Tell all components that the simulation is over
	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->simulationEnded();
	}

	// Tell model also
	model()->simulationEnded();
}

// Notify that state vector was changed and derivatives are needed
void ODESolver::computeDerivatives()
{
	notifyOnStateVectorChanged();
	notifyOnComputeDerivatives();
	_derivativeEvals++;
}

// allocate a numeric vector initialized to zeros
Number* ODESolver::allocateNumVector()
{
	int				n = stateVectorSize();
	Number*			vect;
	
	vect = new Number[n];
	zeroNumVector(vect);
	return vect;
}

// delete a numeric vector (to go with allocate)
void ODESolver::deleteNumVector(Number* vect)
{
	delete[] vect;
}

// Set all the entries of a numeric vector to zero
void ODESolver::zeroNumVector(Number* vect)
{
	int n = stateVectorSize();
	int k;

	for (k=0;k<n;k++) vect[k]=0;
}

// Multiply each element of x by a and add to y.
// The first entry in x (i.e. x[0]) is not used
// because the first entry is reserved in state vectors
// to allow for cases where no offset was assigned.
// Also the MS debugger seems to have some problems
// showing x[0] correctly when code is optimized.
void ODESolver::axpy(Number a, Number* x, Number* y)
{
	int n = stateVectorSize();
	int k = 1;

	// Loop unrolling -- uncomment if needed
	// for (k=1;k<n-3;k+=4) {
	//	y[k]+=a*x[k];
	//	y[k+1]+=a*x[k+1];
	//	y[k+2]+=a*x[k+2];
	//	y[k+3]+=a*x[k+3];
	// }

	// One at a time loop
	for (;k<n;k++) {
		y[k]+=a*x[k];
	}
}

// Multiply each element of x by a scalar.
// x[0] is ignored (see axpy)
void ODESolver::scalex(Number a, Number* x)
{
	int n = stateVectorSize();
	int k = 1;

	// Loop unrolling -- uncomment if needed
	// for (k=1;k<n-3;k+=4) {
	//	x[k]=a*x[k];
	//	x[k+1]=a*x[k+1];
	//	x[k+2]=a*x[k+2];
	//	x[k+3]=a*x[k+3];
	// }

	// One at a time loop
	for (;k<n;k++) {
		x[k]=a*x[k];
	}

}

// Copy numbers from x to y.
// Entry x[0] is ignored (see axpy).
void ODESolver::cpxy(Number* x, Number* y)
{
	int n = stateVectorSize();
	int k = 1;

	// Loop unrolling -- uncomment if needed
	// for (k=1;k<n-3;k+=4) {
	//	y[k]=x[k];
	//	y[k+1]=x[k+1];
	//	y[k+2]=x[k+2];
	//	y[k+3]=x[k+3];
	// }

	// One at a time loop
	for (;k<n;k++) {
		y[k]=x[k];
	}
}

// Weighted distance between vector x andx vector y.
// Distance function is infinity norm.
// Entry x[0], y[0], and w[0] are ignored (see axpy).
Number ODESolver::wdistxy(Number* w, Number* x, Number* y)
{
	int n = stateVectorSize();
	int k = 1;
	Number maxDist = 0;
	Number d0;

	// Loop unrolling -- uncomment if needed
	// Number d1,d2,d3;
	// for (k=1;k<n-3;k+=4) {
	//
	//	d0=fabs(w[k]*(x[k]-y[k]));
	//	if (d0>maxDist) maxDist = d0;
	//
	//	d1=fabs(w[k+1]*(x[k+1]-y[k+1]));
	//	if (d1>maxDist) maxDist = d1;
	//
	//	d2=fabs(w[k+2]*(x[k+2]-y[k+2]));
	//	if (d2>maxDist) maxDist = d2;
	//
	//	d3=fabs(w[k+3]*(x[k+3]-y[k+3]));
	//	if (d3>maxDist) maxDist = d3;
	// }

	// One at a time loop
	for (;k<n;k++) {
		d0=fabs(w[k]*(x[k]-y[k]));
		if (d0>maxDist) maxDist = d0;
	}
	
	return maxDist;
}

// Weighted distance between vector x andx vector y.
// Distance function is the L2 norm (sum of squares).
// Entry x[0], y[0], and w[0] are ignored (see axpy).
Number ODESolver::wdist2xy(Number* w, Number* x, Number* y)
{
	int n = stateVectorSize();
	int k = 1;
	double squareDist = 0;
	double d0;

	// Loop unrolling -- uncomment if needed
	// double d1,d2,d3;
	// for (k=1;k<n-3;k+=4) {
	//
	//	d0=fabs(w[k]*(x[k]-y[k]));
	//	squareDist += d0*d0;
	//
	//	d1=fabs(w[k+1]*(x[k+1]-y[k+1]));
	//	squareDist += d1*d1;
	//
	//	d2=fabs(w[k+2]*(x[k+2]-y[k+2]));
	//	squareDist += d2*d2;
	//
	//	d3=fabs(w[k+3]*(x[k+3]-y[k+3]));
	//	squareDist += d3*d3;
	// }

	// One at a time loop
	for (;k<n;k++) {
		d0=fabs(w[k]*(x[k]-y[k]));
		squareDist += d0*d0;
	}
	
	return sqrt(squareDist);
}



// ====================================================================
// ClockSolver class body
// ====================================================================



// Initialize a new instance
ClockSolver::ClockSolver() {}

// Destroy this instance
ClockSolver::~ClockSolver() {}

// Notify components of start
void ClockSolver::start()
{
	// Set the starting time (once)
	currentTime( beginTime() );

	// Tell model we are starting
	model()->simulationStarted();

	// Tell components we are starting
	notifyOnSimulationStarted();
}

// Notify component of end
void ClockSolver::finish()
{
	notifyOnSimulationEnded();
}

// Notify all components of a clock tick(s)
// This notifies all components regardless of
// of interest expressed by the component.
void ClockSolver::processStep(SimTime hMax)
{
	SimTime stopTime = currentTime() + hMax;

	while (currentTime()<stopTime) {

		// Bump up current time to account for step taken
		currentTime( currentTime() + timeStep());

		// Inform components that the step is now done
		notifyOnTimeStepEnded();
	}
}

// Notify components that the time step has ended.
// Any events should be raised at this point.
// Data unique to prior step can be released and any
// necessary setup for the next step done.
void ClockSolver::notifyOnTimeStepEnded()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = model()->components();

	// Tell model time step is about to end
	model()->prepareTimeStepEnded( currentTime() );

	// Inform all components that the step has ended
	for (it=comps.begin();it!=comps.end();it++) {
		(*it)->timeStepEnded();
	}

	// Inform model that time step has ended
	model()->timeStepEnded();

	// Increment the performance counter
	_timeStepsDone++;
}




// ====================================================================
// EulerMethodSolver class body
// ====================================================================



// Create a new instance
EulerMethodSolver::EulerMethodSolver()
{
	// set default time step (arbitrary)
	_timeStep = 10*UOM::microsec; 
}

// Destroy the instance and any state vectors
EulerMethodSolver::~EulerMethodSolver() {}

// Take a single simulation step
void EulerMethodSolver::processStep(SimTime hMax)
{
	SimTime hSave = timeStep();
	SimTime h;

	// Choose the current time step
	if (hMax < hSave) {
		timeStep(hMax);
	}
	h = timeStep();

	// Do the Euler step
	axpy(h, model()->derivVector(), model()->stateVector() );

	// Bump up current time to account for step taken
	_currentTime += h;

	// Compute derivatives for the next round
	computeDerivatives();

	// Inform components that the step is now done
	notifyOnTimeStepEnded();

	// Restore the original time step
	timeStep(hSave);
}



// ====================================================================
// RungeKuttaSolver class body
// ====================================================================



// Create and initialize a new instance
RungeKuttaSolver::RungeKuttaSolver()
{
	// Set default time step (arbitrary)
	_timeStep = 25*UOM::microsec; 
}

// Destroy the instance and any of its state vectors
RungeKuttaSolver::~RungeKuttaSolver() {}

// Take a single simulation step
void RungeKuttaSolver::processStep(SimTime hMax)
{
	// Compute using a Runge-Kutta 4th order method.
	//
	// k1 = f(t,y)
	// k2 = f(t+h/2,y+h*k1/2)
	// k3 = f(t+h/2,y+h*k2/2)
	// k4 = f(t+h,y+h*k3)
	// y <- y+h*(k1+2*k2+2*k3+k4)/6

	SimTime			t = currentTime();				// starting time
	SimTime			hSave = timeStep();				// time step at start
	SimTime			h;								// time step to use

	// Work vectors
	Number*			y0 = allocateNumVector();
	Number*			k1 = allocateNumVector();
	Number*			k2 = allocateNumVector();
	Number*			k3 = allocateNumVector();
	Number*			k4 = allocateNumVector();


	// Choose the time step as the smaller of the current step size
	// and the parameter hMax
	if (hMax<hSave ) {
		timeStep(hMax);
	}
	h = timeStep();

	// Save current state values in y0
	cpxy(model()->stateVector(),y0);

	// Save starting deriv as k1
	cpxy(model()->derivVector(),k1);

	// Compute k2=f(t+h/2,y+h*k1/2);
	axpy(h/2,k1,model()->stateVector());
	currentTime(t+h/2);
	computeDerivatives();
	cpxy(model()->derivVector(),k2);

	// Compute k3=f(t+h/2,y+h*k2/2);
	cpxy(y0,model()->stateVector());
	axpy(h/2,k2,model()->stateVector());;
	computeDerivatives();
	cpxy(model()->derivVector(),k3);

	// Compute k4=f(t+h,y+h*k3)
	cpxy(y0,model()->stateVector());
	axpy(h,k3,model()->stateVector());
	currentTime(t+h);
	computeDerivatives();
	cpxy(model()->derivVector(),k4);

	// Compute the next state from the formula
	axpy(2,k2,k1);
	axpy(2,k3,k1);
	axpy(1,k4,k1);
	cpxy(y0,model()->stateVector());
	axpy(h/6,k1,model()->stateVector());

	// Compute derivatives for next time around
	computeDerivatives();

	// Inform components that the step is now done
	notifyOnTimeStepEnded();

	// Restore the original time step if changed
	timeStep(hSave);

	// Delete work vectors
	deleteNumVector(y0);
	deleteNumVector(k1);
	deleteNumVector(k2);
	deleteNumVector(k3);
	deleteNumVector(k4);
}



// ====================================================================
// AdaptiveSolver class body
// ====================================================================



// Create and initialize a new instance
AdaptiveSolver::AdaptiveSolver()
{
	// Set default time step (arbitrary values)
	timeStepLessOrEqual(1*UOM::microsec);	// arbitrary small starting point

	// Set default parameters controlling step size increase
	_stableStepLimit = 16;
	_stableSteps = 0;

}

// Destroy the instance and any of its state vectors
AdaptiveSolver::~AdaptiveSolver() {}

// Take a single simulation step of maximum duration hMax.
// Adjust the size of the step taken and order of method used
// based on error tolerances.
void AdaptiveSolver::processStep(SimTime hMax)
{
	// Misc constants
	const Number		epsilon = Number (1e-12);
	const Number		hAllowance = 1e-3f;
	const int			maxOrder = 3;

	// Constants controlling step size changes
	const Number		maxRelStepChg = 2;
	const Number		maxRelStepOnFailure = 0.75;
	const Number		minRelStepOnFailure = 0.125;

	// Initializations
	const SimTime		t0 = currentTime();
	const SimTime		hSave = timeStep();

	// Working variables
	SimTime				h,hEst,hFom,hLimit;
	SimTime				hEst1 = 0;
	SimTime				hEst2 = 0;
	SimTime				hEst3 = 0;
	double				errEst;

	Number*				y = model()->stateVector();
	Number*				yDot = model()->derivVector();
	Number				dist;

	Number*				y0		= allocateNumVector();
	Number*				y0Dot	= allocateNumVector();
	Number*				y1Dot	= NULL;		// allocate only if needed
	Number*				y2Dot	= NULL;		// allocate only if needed

	bool				stepSizeDecreased = false;
	int					stepOrder;

	// Choose the starting step size constrained by hMax.
	// In some cases hMax is smaller than the time step because
	// we are taking a small step to synchronize the time
	// time step boundaries across multiple solvers.
	// A small extra allowance over hSave is used to allow
	// catching up after round-off errors.
	if (hMax < hSave*(1+hAllowance)) {
		h = hMax;
	}
	else {
		h = hSave;
	}

	// Save starting state in case step size must be adjusted and
	// save initial deriv (from last step) for P/C calculation below.
	cpxy(y, y0);
	cpxy(yDot, y0Dot);

	// Try ever smaller step sizes until error tolerance is satisfied.
	// For each step, try successively higher order methods until one
	// is found that meets the error tolerance. When the stable step
	// limit is reached, try all available orders. This prevents being
	// stuck in a low order method that is working but not that well.
	for(;;) {

		// Compute predictor y1=y0+h*y0Dot as an estimate of the state at
		// time t0+h and then get a new derivative vector at that state
		axpy(h,y0Dot,y);
		currentTime(t0 + h);
		computeDerivatives();

		// See if the first order Euler solution is close enough.
		// The error estimate is h*norm(y2Dot-y1Dot)/2
		// where the error scales according to h^2.
		dist = wdistxy(model()->weightVector(),yDot,y0Dot);
		dist = dist<epsilon ? epsilon : dist;
		errEst = h*dist/2;
		hEst = hEst1 = h*safetyMargin()*sqrt(errTol()/errEst);
		if (_stableSteps < stableStepLimit() ) {
			if (errEst<=errTol() ) {
				// Done with the step using a first order method.
				// Note that errors are second order.
				stepOrder = 1;
				break;
			}
		}
	
		// Compute the corrector y2=y0+h*(y0Dot+y1Dot)/2
		// The actual computation is y2=y1-h/2*y0Dot+h/2*y1Dot;
		if (y1Dot==NULL) {
			y1Dot=allocateNumVector();
		}
		cpxy(yDot,y1Dot);
		axpy(-h/2,y0Dot,y);
		axpy(h/2,y1Dot,y);

		// Get new derivative at state y2 and leave in yDot
		// Note that current time is the same as before.
		computeDerivatives();

		// Compute estimated step size needed.
		// The error estimate is h*norm(y2Dot-y1Dot)/3
		// where the error scales according to h^3.
		dist = wdistxy(model()->weightVector(), yDot, y1Dot);
		dist = dist<epsilon ? epsilon : dist;
		errEst = h*dist/3;
		hEst=hEst2=h*safetyMargin()*pow(errTol()/errEst,1/3.0);
		if (_stableSteps < stableStepLimit() ) {
			if (errEst<=errTol() ) {
				// Done with the step using a second order method
				// Note that errors are third order.
				stepOrder = 2;
				break;
			}
		}

		// Press on to third order method. New value of y is
		// y3 = y0+h*(y0Dot/2+y1Dot/6+y2Dot/3). 
		// Actual computation is y3 = y2+h*(-y1Dot+y2Dot)/3
		if (y2Dot==NULL) {
			y2Dot=allocateNumVector();
		}
		cpxy(yDot,y2Dot);
		axpy(-h/3,y1Dot,y);
		axpy(h/3,y2Dot,y);

		// Get new derivative at state y3 and leave in yDot
		// Note that current time is the same as before.
		computeDerivatives();

		// Compute estimated step size needed.
		// The error estimate is h*norm(y3Dot-y2Dot)/4
		// where the error scales according to h^4.
		dist = wdistxy(model()->weightVector(), yDot, y2Dot);
		dist = dist<epsilon ? epsilon : dist;
		errEst = h*dist/4;
		hEst=hEst3=h*safetyMargin()*pow(errTol()/errEst,1/4.0);
		if (errEst<=errTol() ) {
			// Done with the step using a third order method
			// Note that errors are fourth order.
			stepOrder = 3;
			break;
		}

		// Failed in the attempt to satisfy error tolerance.
		// Try to reduce step size and try again
		// Make sure hEst is actually smaller than h but avoid
		// cases where hEst is too much a change from h in case
		// there is a degenerate condition where hEst ~= 0.
		if (hEst>=h) {
			hEst = maxRelStepOnFailure*h;
		}
		else if (hEst<minRelStepOnFailure*h) {
			hEst = minRelStepOnFailure*h;
		}

		// Make sure step size can be reduced
		if (hEst>=minTimeStep() ) {

			// Restore the initial state
			cpxy(y0,y);

			// Try a smaller step size based on the estimated ideal
			h = hEst;
			stepSizeDecreased = true;
		}
		else {
			// Can't decrease step size, so stop here
			h = minTimeStep();
			break;
		}
	}

	// Inform components that the step is now done
	notifyOnTimeStepEnded();

	// Is increasing the step size a possibility or not
	if (stepSizeDecreased) {

		// Choose a new step size based on what worked the last time.
		timeStepLessOrEqual(h);		

		// This is the end of a run of stable steps
		_stableSteps = 0;
	}
	else {

		if (_stableSteps++ >= stableStepLimit() ) {
			// Start over counting stable steps
			_stableSteps = 0;
		}

		// Avoid trying to make too radical a change in step size
		// limiting maximum new step size relative to current size.
		hLimit = maxRelStepChg*hSave;
		if (hEst1>hLimit) 
			hEst1=hLimit;
		if (hEst2>hLimit) 
			hEst2=hLimit;
		if (hEst3>hLimit)
			hEst3=hLimit;

		// Determine a figure of merit based on net processing speed
		// for each different step order and pick the best step size.
		hFom = hEst1;
		hEst = hEst1;
		if (stepOrder>=2 && hEst2/2 > hFom) {
			hFom = hEst2/2;
			hEst = hEst2;
		}
		if (stepOrder>=3 && hEst3/3 > hFom) {
			hEst = hEst3;
		}

		// Finally, set step size for next time step staying within bounds.
		if (hEst <= minTimeStep() ) {
			timeStep( minTimeStep() );
		}
		else if (hEst >= maxTimeStep() ) {
			timeStep( maxTimeStep() );
		}
		else {
			timeStepLessOrEqual(hEst);
		}
	}

	// Free work vectors
	deleteNumVector(y0);
	deleteNumVector(y0Dot);
	deleteNumVector(y1Dot);
	deleteNumVector(y2Dot);
}



// ====================================================================
// Event class body
// ====================================================================



// Create a new instance with event data
Event::Event(SimTime t)
{
	_useCount = 0;
	_eventTime=t;
}

// Destroy the instance
Event::~Event() {}

// Perform any unique actions when the event is acted on
// If used, this must be overridden in subclasses.
void Event::dispatch()
{
	FatalError("(Event::dispatch) Subclass must override this function");
}



// ====================================================================
// Probe class body
// ====================================================================



Probe::Probe()
{
	using namespace UOM;

	_minInterval = 0*sec;
	_flushInterval = 1*sec;
	_intervalStart = InfinitePast;
	_flushTime = InfinitePast;
}

Probe::~Probe() {}

// Determine if time t is good for reporting
bool Probe::isReportable(SimTime t)
{
	SimTime mi = minInterval();
	SimTime intervalEnd = _intervalStart+mi;

	// Has the current interval expired
	if (intervalEnd > t)
		return false;

	// If t is beyond the end of the current interval but not
	// beyond the end of the next one, just continue the intervals
	// with a period of minInterval. Otherwise, start a new interval
	// at time t.
	if (intervalEnd+mi > t) {
		_intervalStart += mi;
	}
	else {
		_intervalStart = t;
	}
	return true;
}

// Handle end of time step for a single model component
void Probe::timeStepEnded(ModelComponent* mc)
{
	SimTime					t=mc->model()->currentTime();

	if (isReportable(t) ) {

		// Get vector of components to probe
		ModelComponentVector ctp;
		ctp.reserve(128);
		mc->addToComponentsToProbe(ctp);

		// Have subclass handle the reporting
		reportOn(ctp, t, mc->numericIdentifier() );

		// Flush if the time has come
		if (flushInterval() + _flushTime <= t) {
			flush();
			_flushTime = t;
		}		
	}
}

// Handle end of time step for a model. All components
// in the model are typically reported on in this case.
void Probe::timeStepEnded(Model* mod, int id)
{
	// Get time
	SimTime t = mod->currentTime();

	if (isReportable(t) ) {

		// Let subclasses do the reporting
		reportOn( mod->components(), t, id );
	}

	// Flush if the time has come
	if (flushInterval() + _flushTime <= t) {
		flush();
		_flushTime = t;
	}		
}



// ====================================================================
// ExternalRecorder class body
// ====================================================================



// Create a new instance
ExternalRecorder::ExternalRecorder(char* fn, char* mapfn)
{
	// set defaults
	_file = NULL;
	_mapFile = NULL;
	_fileName = NULL;
	_mapFileName = NULL;

	fileName(fn);
	mapFileName(mapfn);

	_headerWritten = false;
}

// Destroy this instance
ExternalRecorder::~ExternalRecorder()
{
	// Close out any files left open
	closeMapFile();
	if (_file != NULL) {
		fclose(_file);
	}
	
	// Delete any held strings
	delete _fileName;
	delete _mapFileName;
}

// Set the file name to write to
void ExternalRecorder::fileName(char* fn)
{
	// Delete any old string
	delete _fileName;

	// Make a local copy of the value passed
	_fileName = new char[strlen(fn)+1];
	strcpy(_fileName,fn);
}

// Save the map file name string
void ExternalRecorder::mapFileName(char* mfn)
{
	// Delete any old string
	delete _mapFileName;

	if (mfn!=NULL) {
		// Make a local copy of the value passed
		_mapFileName = new char[strlen(mfn)+1];
		strcpy(_mapFileName,mfn);
	}
	else {
		_mapFileName = NULL;
	}
}

// Open output file if needed 
void ExternalRecorder::simulationStarted()
{
	if (_file == NULL) {
		// Open the output file
		_file = fopen(_fileName,"w");
		if (_file==NULL) {
			cerr << "Could not open file for output."<<endl;
			cerr << "  file name = " << _fileName <<endl;
			FatalError("(ExternalRecorder::simulationStarted) could not open output file");
		}
	}
}

// Close output file if any
void ExternalRecorder::simulationEnded()
{
	if (_file != NULL) {
		fclose(_file);
	}
	_file = NULL;
	closeMapFile();
}

// Open the map file for output
void ExternalRecorder::openMapFile()
{
	if (_mapFileName != NULL) {
		_mapFile = fopen(_mapFileName,"w");
		if (_mapFile==NULL) {
			cerr << "Could not open file for output."<<endl;
			cerr << "  file name = " << _mapFileName <<endl;
			FatalError("(ExternalRecorder::openMapFile) error on file");
		}
	}
}

// Close the map file
void ExternalRecorder::closeMapFile()
{
	if (_mapFile != NULL) {
		fclose(_mapFile);
	}
	_mapFile = NULL;
}

// Flush the output buffer
void ExternalRecorder::flush()
{
	if (_file != NULL) {
		fflush(_file);
	}
}



// ====================================================================
// ExternalStateRecorder class body
// ====================================================================



// Create a new instance
ExternalStateRecorder::ExternalStateRecorder(char* fn, char* mfn)
: ExternalRecorder(fn,mfn) {}

// Destroy this instance
ExternalStateRecorder::~ExternalStateRecorder() {}

// Write out the state for a collection of components
void ExternalStateRecorder::reportOn(ModelComponentVector& comps, SimTime t, int id)
{
	// Make sure file was opened before reporting
	if (_file != NULL) {

		// Write a header line if needed
		if (!_headerWritten) {
			writeColumnHeader(comps);
			_headerWritten = true;
		}

		// Write the current state values
		writeState(comps,t,id);
	}
	else {
		FatalError("(ExternalStateRecorder::reportOn) output file is not open");
	}
}

// Write out a column header line for documentation
void ExternalStateRecorder::writeColumnHeader(ModelComponentVector& comps)
{
	ModelComponentVectorIt		it;
	int							k,n;
	const char*					name;
	const char**				labels;
	int							colNbr = 0;

	// If a map file was requested, open a file
	if (_mapFileName != NULL) {
		openMapFile();
	}

	// Write standard label for id and simulation time
	fprintf(_file,"id,t");
	if (_mapFile!=NULL) {
		fprintf(_mapFile,"col\tlabel\n1\tid\n2\tt\n");
		colNbr=2;
	}

	// Get state vector labels from components
	for (it=comps.begin();it!=comps.end();it++) {

		// Make name for each state variable ofthe form:
		// "compName_svName". Labels of this form are 
		// not quaranteed to be unique but serve as useful
		// documentation.

		name=(*it)->componentName();
		labels=(*it)->stateLabels();
		n=(*it)->numAllStateVar();

		for (k=0;k<n;k++) {
			fprintf(_file,",%s_%s",name,labels[k]);
			if (_mapFile != NULL) {
				fprintf(_mapFile,"%d %s_%s\n",++colNbr,name,labels[k]);
			}
		}
	}

	// Write end of header line
	fprintf(_file,"\n");

	// Wrap up map file, if any
	closeMapFile();

}

// Write out the current state values preceeded by
// the component identifier and time. Time is converted
// to milliseconds. Other state values are converted
// based on model component units of measure values.

void ExternalStateRecorder::writeState(ModelComponentVector& comps,SimTime t,int id)
{
	ModelComponentVectorIt		it;
	int							k,n;

	// Write out identifier and time
	fprintf(_file, "%d,%.12g",id,t/UOM::msec);
	
	// Get state values from components
	for (it=comps.begin();it!=comps.end();it++) {

		// Write out each state value separated by a comma
		n=(*it)->numAllStateVar();
		for (k=0;k<n;k++) {
			fprintf(_file,",%g",double((*it)->stateValueConverted(k)) );
		}
	}

	// Write the end of line
	fprintf(_file,"\n");
}
