// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_nmod.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file provides the body of BNSF classes.
// The sequence is generally the same as the definitions
// in the associated header file, but a good PDE is
// almost a necessity for any reasonable development effort.
//
// This file contains class bodies associated with neural models.


#include "bnsf_nmod.h"
#include <utility>
#include <stack>


using namespace std;
using namespace BNSF;



// ====================================================================
// NeuronSolver class body
// ====================================================================



// Static data initializations
const int NeuronSolver::_PerPassArraySize = 8;

// Constructor
NeuronSolver::NeuronSolver() 
{
	int			i;

	// Clear step counts (points for neatness)
	for (i=0;i<_PerPassArraySize;i++) {
		_nsteps[i]=0;
	}

	// Initialize current number of passes and related values
	_minSteps = 2;	
	_curPasses = 4;
	_nsteps[0]=9;
	_nsteps[1]=17;
	_nsteps[2]=33;
	_nsteps[3]=65;

	// Other initializations
	_errorEst = 0;
	_compSVOffset = 0;
	_compSVSize = 0;

	// Set other parameter defaults via accessors
	errTol(.001f);					// default error tolerance
	maxPasses(4);					// max passes used in extrapolation	
	maxSteps(1025);					// max steps per pass before time step is changed
	safetyMargin(0.875);			// error relative margin before decreasing order
	timeStep(250*UOM::microsec);	// default macro time step
	assumeFirstOrder(false);		// allow extrapolations to assume second order
	alwaysRebuildJacobian(false);	// save jacobian between steps - avoid rebuild
}

// Destructor
NeuronSolver::~NeuronSolver() {}

// Set maximum passes value
void NeuronSolver::maxPasses(int n)
{
	if (n>_PerPassArraySize) {
		FatalError("(NeuronSolver::maxPasses) Value exceeds maximum allowed.");
	}
	_maxPasses=n;
}

// Set minimum steps value
void NeuronSolver::minSteps(int n)
{
	int k;

	if (n<2) {
		FatalError("(NeuronSolver::minSteps) Must allow at least 2 steps per pass");
	}

	// Save the new minimum
	_minSteps=n;

	// Reset number of steps per pass if needed
	if (_nsteps[0]<_minSteps) {
		_nsteps[0]=n;
		for (k=1;k<_curPasses;k++) {
			_nsteps[k]=2*_nsteps[k-1];
		}
	}
}

// Do preliminary initializations (once only)
void NeuronSolver::start()
{
	// Tell model we are starting
	model()->simulationStarted();

	// Build lists of who gets what message
	prepareRecipients();

	// Set the starting time (once)
	currentTime( beginTime() );

	// Tell components we are starting
	notifyOnSimulationStarted();

	// Set the Jacobian indexes in all compartments
	assignJacobianIndexes();

	// Have components set initial conditions
	notifyOnSetInitialState();

	// State vector is now set first time, so notify affected components.
	notifyOnStateVectorChanged();

	// Compute the initial derivatives
	computeDerivatives();

	// Inform components that initialization is done
	notifyOnTimeStepEnded();

}

// Build vector of recipients of each type of notification by
// querying each component to see what notifications it should receive.
void NeuronSolver::prepareRecipients()
{
	ModelComponentVectorIt		it;
	ModelComponentVector&		comps = model()->components();

	// Build these as we are going along
	_compSVOffset = -1;
	_compSVSize = 0;

	// Customize prepare to handle compartments specially
	for (it=comps.begin();it!=comps.end();it++) {
		if ( (*it)->isCompartment() ) {

			// Save compartments here.
			// By default, all compartments get notification
			// of state changes so they do not need to be
			// included in the vector built below.
			_compartments.push_back(reinterpret_cast<Compartment*>(*it));

			// Save the state vector offset to the first compartment.
			// The other ones should follow this in order in the state vector.
			if (_compSVOffset == -1) {
				_compSVOffset = (*it)->svOffset();
			}

			// Add up the size of the state vector and deriv vector
			// associated with compartments (assumed the same)
			_compSVSize += (*it)->numStateVar();
		}
		else {

			// Compartments all receive notification of state change.
			// Find any other components also requesting this notification.
			if ( (*it)->notifyOnStateVectorChanged() ) {
				_stateVectorChangedRecipients.push_back(*it);
			}

			// For the moment, every component with state must be able 
			// to do local updates. Some general methods might be added later.
			if ( (*it)->numStateVar()>0 &&
				!(*it)->canPerformLocalUpdate() ) {
				FatalError("(NeuronSolver::prepareRecipients) "
					"Found non-compartment component that cannot perform local state update.");
			}

			// Only need to save components that can do local updates
			if ((*it)->canPerformLocalUpdate() ) {

				// See whether this components is evaluated at the same
				// time points as compartments or half a step out of phase.
				if ( (*it)->isAlignedForCN() ) {
					_timeAlignedComponents.push_back(*it);
				}
				else {
					_timeOffsetComponents.push_back(*it);
				}
			}
		}

		// Any component may want notification of time step end
		if ( (*it)->notifyOnTimeStepEnded() ) {
			_timeStepEndedRecipients.push_back(*it);
		}
	}

	// Reallocate vectors to save memory
	CompartmentVector tempCompVect(_compartments);
	swap(tempCompVect,_compartments);

	ModelComponentVector tempSVChg(_stateVectorChangedRecipients);
	swap(tempSVChg,_stateVectorChangedRecipients);

	ModelComponentVector tempAligned(_timeAlignedComponents);
	swap(tempAligned,_timeAlignedComponents);

	ModelComponentVector tempOffset(_timeOffsetComponents);
	swap(tempOffset,_timeOffsetComponents);
}

// Build vector of micro time step sizes
void NeuronSolver::buildHvect(double H)
{
	int			i;

	for (i=0;i<_curPasses;i++) {
		if (_nsteps[i]<2) {
			FatalError("(NeuronSolver::buildHvect) Cannot have fewer than 2 step per pass");
		}
		else {
			_hvect[i]=H/(_nsteps[i]-0.5);
		}
	}

	// Just to be safe, clear the remaining entries to zero
	for (;i<_PerPassArraySize;i++) {
		_hvect[i]=0;
	}
}

// Build extrapolation coefficients given step sizes
void NeuronSolver::buildExtrapolationCoeff()
{
	// A reference on Richardson Extrapolation is:
	// Liem CB,Lu T, and Shih Tm (1995). 
	// The Splitting Extrapolation Method. 
	// Singapore: World Scientific Publishing.

	// Two version of the extrapolation coefficients are built.
	// coeff[0] arrays are for objects that are first order accurate,
	// while coeff[1] arrays are for use by objects that are second
	// order accurate. Both versions are built in the following loop.

	int			i,j,k;
	double		h,hToNth;

	// Matrices and vectors for state and error extrapolation.
	// Sparse matrices are less efficient but keep the implementation
	// self contained without requiring a full matrix support package.

	SparseMatrix	Vst(_curPasses,_curPasses);
	SparseMatrix	Verr(_curPasses-1,_curPasses-1);

	double		stcoeff[_PerPassArraySize];
	double		errcoeff[_PerPassArraySize-1];

	double		a0st[_PerPassArraySize] = {0};
	double		a0err[_PerPassArraySize-1] = {0};

	a0st[0]=1;
	a0err[0]=1;

	// Create coefficients for 1st (k=0) and 2nd (k=1) order
	// extrapolation polynomials.
	for (k=0;k<2;k++) {	

		// Use Vandermonde matrixes for polynomial interpolation.
		// The error coefficients are the difference between those
		// of the number of passes actually used and those of one
		// fewer passes by leaving out the last h value.
		// (see eq. 1.1.54 in Liem). 
		
		// Vandermonde matrixes can be poorly conditioned and hence
		// numerically inaccurate, but only a small number of terms 
		// are used in the extrapolation. Double precision arithmetic
		// should suffice in this case.

		// This process could be made much more efficient by caching
		// generated coefficients. Consider this for a future upgrade.

		// Start by zeroing out all coefficients (just to make sure)
		for (i=0;i<_PerPassArraySize;i++) {
			_errorExtrap[k][i]=0;
			_stateExtrap[k][i]=0;
		}
		for (i=0;i<_curPasses;i++) {

			// Initialize for successive powers of h.
			h=_hvect[i]/_hvect[0];
			hToNth =1 ;

			// For extrapolation of second order accurate state values,
			// drop the h^1 from the extrapolation polynomial since its
			// coefficient is known to be zero (nor nearly so).
			if (k==1) {
				hToNth = h;
			}

			// Build the matrix. To facilitate solving for the h^0
			// coefficients, the matrix is built as the transpose
			// of the usual form of the Vandermonde matrix.
			Vst[0][i]=1;
			for (j=1;j<_curPasses;j++) {
				hToNth *= h;
				Vst[j][i]=hToNth;
			}
		}

		// To estimate error, drop the first (fewest steps) pass and
		// form an estimate of one smaller degree of accuracy than the
		// final estimate. The difference is a estimate of error.
		// See Liem et al. for the justification of this method.
		for (i=0;i<_curPasses-1;i++) {
			for (j=0;j<_curPasses-1;j++) {
				Verr[j][i]=Vst[j][i+1];
			}
		}

		// Solve for the h^0 coefficients (polynomial value at h=0). 
		// These would be the first columns of the inverse of the
		// corresponding (transposed) Vandermonde matrixes.
		Vst.luSolve(stcoeff,a0st);
		Verr.luSolve(errcoeff,a0err);

		// Copy the extrapolation estimates to C arrays.
		// Subtract the all-h case from the leave-out-one-h case 
		// with the difference being the error estimate coefficients.
		for (i=0;i<_curPasses;i++) {
			_stateExtrap[k][i] = stcoeff[i];
			_errorExtrap[k][i] = stcoeff[i] - (i>0 ? errcoeff[i-1] : 0);
		}
	}
}

// Assign Jacobian indexes to minimize LU decomp fill.
// This assigns the most distal nodes the lowest indexes.
// This is roughly the Hines method except that no special
// effort is made to limit bandwidth of the resulting matrix.
void NeuronSolver::assignJacobianIndexes()
{
	typedef pair<Number,Compartment*>	sortPairType;
	typedef vector<sortPairType>		sortVectorType;
	typedef sortVectorType::iterator	sortVectorTypeIt;

	sortPairType			sortEnt;
	sortVectorType			sortVect;
	sortVectorTypeIt		sortIt;

	CompartmentVectorIt		it;
	bool					distSet=false;
	int						nComp=0;
	int						jidx=0;

	// Note: the following uses distance from soma as a means of
	// ordering compartments such that the compartments are at
	// least topologically sorted. A better design would be to
	// do the topological sort here or explicitly compute path
	// distance and sort. The soma could be located as the largest
	// compartment and is, in any case, only one possible root
	// of the tree structure. However, there is no good reason
	// to make these changes at present. Maybe later.

	// Get distance from soma and sort by it (descending)
	for (it=_compartments.begin(); it!=_compartments.end(); it++) {
		sortEnt.first = -( (*it)->distFromSoma() );
		sortEnt.second = *it;
		sortVect.push_back(sortEnt);
		distSet |= (sortEnt.first!=0);
		nComp++;
	}
	if (nComp>1 && !distSet) {
		FatalError("(NeuronSolver::assignJacobianIndexes) "
			"Distance from soma has not bet set for compartments.");
	}
	sort(sortVect.begin(),sortVect.end());

	// Assign the Jacobian indexes
	for (sortIt=sortVect.begin();sortIt!=sortVect.end();sortIt++) {
		sortIt->second->jacobianIndex( jidx++);
	}
}

// Build the Jacobian Matrix
void NeuronSolver::buildJacobian()
{
	int			N = _compartments.size();

	CompartmentVectorIt		it;

	// If there is a current matrix, clear the data without
	// releasing storage unecessarily.
	if (_jacobian.dim1()>0) {
		_jacobian.clear();
	}
	if (_jacobian.dim1()!=N) {
		// Changing size could only occur the first time or
		// if compartments changed between steps.
		_jacobian.resize(N,N,0);
	}

	// Update model's time
	model()->currentTime( currentTime() );

	// Build Jacobian rows compartment by compartment
	for (it=_compartments.begin();it!=_compartments.end();it++) {
		(*it)->buildJacobian(_jacobian);
	}
}

// Update the Jacobian Matrix
void NeuronSolver::updateJacobian()
{
	CompartmentVectorIt		it;

	// Update model's time
	model()->currentTime( currentTime() );

	// Update changed entries in Jacobian
	for (it=_compartments.begin();it!=_compartments.end();it++) {
		(*it)->updateJacobian(_jacobian);
	}
}

// Notify compartments that their state vector values changed.
void NeuronSolver::notifyOnStateVectorChanged()
{
	CompartmentVectorIt		it;

	// Make sure model has current time
	model()->currentTime( currentTime() );

	// Since compartments were not included in recipients, 
	// do them explicitly
	for (it=_compartments.begin();it!=_compartments.end();it++) {
		(*it)->stateVectorChanged();
	}

	ModelComponentVectorIt		mcit;
	ModelComponentVector&		comps = _stateVectorChangedRecipients;

	// Tell any other interested components
	for (mcit=comps.begin();mcit!=comps.end();mcit++) {
		(*mcit)->stateVectorChanged();
	}
}

// Update state via local methods in non-compartment components
// which are (except for the first step) evaluated at a time offset.
void NeuronSolver::updateOffsetStates(SimTime h, CNStepType stepType)
{
	ModelComponentVector&	comp = _timeOffsetComponents;
	ModelComponentVectorIt	compIt;

	// Give model the current time
	model()->currentTime( currentTime() );

	// Relay the update request to affected components
	for (compIt=comp.begin(); compIt!=comp.end(); compIt++) {
		(*compIt)->localStateUpdate(h, stepType);
	}
}

// Update state via local methods in non-compartment components
// which are (except for the first step) evaluated at a time offset.
void NeuronSolver::updateAlignedStates(SimTime h, CNStepType stepType)
{
	ModelComponentVector&	comp = _timeAlignedComponents;
	ModelComponentVectorIt	compIt;

	// Give model the current time
	model()->currentTime( currentTime() );

	// Relay the update request to affected components
	for (compIt=comp.begin(); compIt!=comp.end(); compIt++) {
		(*compIt)->localStateUpdate(h, stepType);
	}
}


// Compute derivatives for compartments
void NeuronSolver::computeDerivatives()
{
	CompartmentVectorIt		it;
	CompartmentVectorIt		compEnd = _compartments.end();

	// Give model the current time
	model()->currentTime( currentTime() );

	// Ask each compartment to compute its derivative
	for (it=_compartments.begin();it!=compEnd;it++) {
		(*it)->computeDerivatives();
	}
}

// Compute the implicit update by solving a linear system.
// This is singled out to allow expansion to cases where
// LU decomp will not work because Jacobian cannot be ordered
// to avoid fill-in or is otherwise not sufficiently accurate.
// See SparseMatrix::pbcgSolve for an alternative to LU decomp.
void NeuronSolver::solveImplicitSystem(SparseMatrix& IhJ, double* compDeriv)
{
	// Solve Ax=b where A is IhJ and b is compDeriv. The output is left in
	// compDeriv to avoid allocating additional storage. IhJ is used as a
	// work area in the process of solving the equation and thus
	// no longer holds it original value at the end of the process.

	// Note that LU decomp is only one possible method. PBCG methods work
	// well but are somewhat slower. However, there is no need to sort
	// compartments when using PBCG, which could be useful in models
	// where there are direct electrical connections that do not form
	// a tree structure. Something to consider later.

	if (!IhJ.luDecompWithoutPivot()) {
		FatalError("(NeuronSolver::solveImplicitSystem) LU decomposition failed");
	}

	// Solve the linear system
	IhJ.luBackSubRight(compDeriv,compDeriv);

}

// Do one time step of the duration indicated.
// A variant of the Crank-Nicolson evaluation scheme is used
// in which evaluation times for compartments are staggered
// by half a time step with those of ion channels and other
// state variables. 
void NeuronSolver::processStep(SimTime maxDuration)
{
	// Direct access to model vectors
	const int		stateVectorSize = model()->stateVectorSize();
	Number*			modelState = model()->stateVector();
	Number*			modelDeriv = model()->derivVector();
	Number*			weight = model()->weightVector();

	// Direct access to the model components vector
	ModelComponentVector& components = model()->components();

	// Copies of the state vectors used in the extrapolation process
	Number*			startingState = NULL;		// all states vector
	Number*			startingDeriv = NULL;		// compartments derivatives

	Number*			endingState[_PerPassArraySize] = {NULL}; // array of arrays

	// Array of indexes in to extrapolation coefficient with one 
	// entry for each state variable. This based on whether or
	// not the associated component is 2nd order accurate.
	short*			orderIndex = NULL;

	// Number of compartments
	const int		numComp = _compartments.size();

	// Derivatives for compartments permuted by Jacobian index order
	double*			compDeriv = NULL;

	// Macro step time. This can be changed if the original
	// interval cannot be crossed within error tolerance.
	SimTime			H = maxDuration;

	// Times for compartments and other components (which can differ)
	SimTime			tStart;		// current time as of step start
	SimTime			tComp;		// current evaluation time for compartments
	SimTime			tOther;		// current time for other components
	SimTime			hComp;		// current time step for compartments
	SimTime			hOther;		// current time step for other components

	// Flag for first time initializations
	bool			initialStep = true;

	// Misc work variables
	int				i,j,k,n;
	int				npass;
	int				nstepsNext;
	int				errorIndex;
	SparseMatrix	IhJ(_compSVSize,_compSVSize,1);
	bool			isStartingHalfStep, isEndingHalfStep;

	// Save the starting state and time
	tStart = currentTime();
	startingState = allocateNumVector();
	cpxy(modelState,startingState);

	// Save the starting derivative values for compartments
	startingDeriv = new Number[_compSVSize];
	for (i=0;i<_compSVSize;i++) {
		startingDeriv[i]=modelDeriv[i+_compSVOffset];
	}

	// Set the indicators for second order accuracy.
	// Any empty slots in the state vector are set to
	// first order to avoid memory referencing exceptions.
	// A possible optimization is to save this vector
	// between steps when it is known that order is unchangeable.

	orderIndex = new short[stateVectorSize];
	for (i=0;i<stateVectorSize;i++) {	// Initialize all entries to 1st order 
		orderIndex[i] = 0;
	}
	if (!assumeFirstOrder() ) {			// Can skip ahead if all are first order
		for (i=0;i<components.size();i++) {	// Set order based on component info

			ModelComponent* pcomp = components[i];
			j = pcomp->svOffset();
			for (k=0;k<components[i]->numStateVar();k++) {
				orderIndex[j++] = pcomp->isSecondOrderAccurate() ? 1 : 0;
			}
		}
	}

	// Allocate a work area for compartment derivative values
	compDeriv = new double[numComp];

	// Allocate vectors for holding ending states
	for (i=0;i<_curPasses;i++) {
		endingState[i] = allocateNumVector();
	}

	// Build vectors containing step sizes and values affected by them
	buildHvect(H);
	buildExtrapolationCoeff();

	// Start with pass 0. Depending on error results, the number
	// of passes may be extended to meet tolerances.
	npass = 0;

	// Loop until error tolerance met
	while (true) {

		// Make one pass for each step size
		while (npass<_curPasses) {

			// Baseline times for each domain
			tComp = tStart;
			tOther = tStart;

			// Except for the first step, where values have not been
			// changed, restore all states to their starting values 
			// and restore compartment derivatives.
			if (!initialStep) {
				cpxy(startingState,modelState);
				for (i=0;i<_compSVSize;i++) {
					modelDeriv[i+_compSVOffset] = startingDeriv[i];
				}
			}

			// Use a first order explicit update to temporarily advance
			// voltages by half of a microstep to facilitate the initial
			// half step for channel states. To be entirely consistent, we
			// should probably also move aligned components half a step
			// ahead temporarily, however, this has not proven necessary
			// in practice and could introduce unnecessary complexity
			// in subordinate component objects.
			hComp=_hvect[npass]/2;
			for (i=0;i<_compSVSize;i++) {
				modelState[i+_compSVOffset] += hComp*startingDeriv[i];
			}
			notifyOnStateVectorChanged();

			// Take micro steps to cross the macro time step interval.
			// Note that compartments and other components are offset
			// by half a time step when crossing the interval except
			// at the beginning and at the end.
			for (n=0;n<_nsteps[npass];n++) {

				// Set time step values for comp and non-comp Test for 
				// half steps at beginning and ending of the pass.
				isStartingHalfStep = n==0;
				isEndingHalfStep = n==_nsteps[npass]-1;

				hComp = isEndingHalfStep ? _hvect[npass]/2 : _hvect[npass];
				hOther = isStartingHalfStep ? _hvect[npass]/2 : _hvect[npass];
			
				// Update states for non-compartmental components
				tOther += hOther;
				currentTime(tOther);
				updateOffsetStates(hOther, 
					isStartingHalfStep ? CNStartingHalfStep : CNFullStep);

				// If this is the starting half step, restore compartment states
				if (isStartingHalfStep) {
					for (i=0;i<_compSVSize;i++) {
						modelState[i+_compSVOffset] = startingState[i+_compSVOffset];
					}
					notifyOnStateVectorChanged();
				}

				// Get compartment derivatives and copy to local
				// vector in Jacobian index order. Allowances are
				// made for compartments to have multiple state
				// variables (just in case -- currently not needed).
				computeDerivatives();
				for (i=0;i<numComp;i++) {
					j=_compartments[i]->jacobianIndex();
					for (k=0;k<_compartments[i]->numStateVar();k++) {
						compDeriv[j+k]=modelDeriv[i+k+_compSVOffset];
					}
				}
		
				// Update or build the Jacobian for compartments
				if (_jacobian.dim1() == 0 ||
					(initialStep && _alwaysRebuildJacobian) ) {
					buildJacobian();
				}
				else {
					updateJacobian();
				}

				// Compute implicit update for compartments.
				// If J is the Jacobian and V the state vector
				// what is being evaluated is:
				// V(t+h)=V(t)+h*inverse(I-h/2*J)*dV/dt
				//
				// When semi-open interval evaluation is done, the
				// final microstep in the passes interpolates the value
				// between two time steps by using:
				// V(t+h)=V(t)+h/2*inverse(I-h/2*J)*dV/dt

				IhJ.setToIdentity();
				_jacobian.addScaledByTo(-_hvect[npass]/2,IhJ);
				solveImplicitSystem(IhJ, compDeriv);

				// Compartment states are updated in two steps.
				// First the state (i.e. voltage) is updated
				// to correspond with the current time (which is
				// offset by half a step) so that currents can be
				// as of the midpoint in the microstep interval.
				// Any aligned components (e.g. calcium pools)
				// are then advanced to the end of the microstep
				// and compartment states are also updated to
				// the end of the microstep. If this is an ending
				// half step, the final update of voltages is not
				// needed and is skipped.

				// Compute updated state (half step) and copy back to model.
				for (i=0;i<numComp;i++) {
					j = _compartments[i]->jacobianIndex();
					for (k=0;k<_compartments[i]->numStateVar();k++) {
						modelState[i+k+_compSVOffset] += _hvect[npass]/2*compDeriv[j+k];
					}
				}
				notifyOnStateVectorChanged();

				// Bring states of other components up to date.
				// Voltages are at the midpoint value as is current time.
				updateAlignedStates(hComp, 
					isEndingHalfStep ? CNEndingHalfStep : CNFullStep);

				// Finish the micro step if not already at the end of the pass.
				if (!isEndingHalfStep) {

					// Do the other half of the update for voltages etc.
					for (i=0;i<numComp;i++) {
						j = _compartments[i]->jacobianIndex();
						for (k=0;k<_compartments[i]->numStateVar();k++) {
							modelState[i+k+_compSVOffset] += _hvect[npass]/2*compDeriv[j+k];
						}
					}

					// Send out notifications of new current time and states.
					tComp += hComp;
					currentTime(tComp);
					notifyOnStateVectorChanged();
				}

				// Phew, one microstep has been done, so count it.
				_derivativeEvals++;
				initialStep = false;
			}

			// Save ending state for this pass and move on to next pass
			cpxy(modelState,endingState[npass]);
			npass++;
		}

		// Compute an error estimate via polynomial extrapolation.
		// Assumed first or second order extrapolation is determined
		// by the orderIndex entry corresponding to each value.
		_errorEst = 0;
		for (i=1;i<stateVectorSize;i++) {
			// Sum error estimate using extrapolation coeff.
			// Since all we want is the worst case error, it is
			// not necessary to build a vector of all errors
			// once we know that error tolerances were exceeded.
			// To find the worst case error for tuning, though,
			// we go ahead and look at all entries.
			double errSum=0;
			for (k=0;k<_curPasses;k++) {
				errSum += _errorExtrap[orderIndex[i]][k]*endingState[k][i];
			}
			errSum *= weight[i];
			if (fabs(errSum)>=_errorEst) {
				errorIndex = i;				// for debug and tuning
				_errorEst = fabs(errSum);
			}
		}

		// If error tolerance is met we can stop the evaluation loop.
		if (_errorEst<errTol() ) {
			break;
		}

		// An idea is to try to detect when the step is actually
		// first order and do the error estimate over again allowing
		// non-zero first order terms. For now this is not yet done.
			
		// Otherwise, could not meet error tolerance constraint.
		// Try to increase steps per pass to reduce error.
		nstepsNext = _nsteps[_curPasses-1]*2-1;

		// Optional debugging display of what caused change
		if (debugPerformance() ) {
			cerr<<"NeuronSolver "
				<<"T="<<tStart
				<<", err="<<_errorEst
				<<", ns["<<_curPasses-1<<"]="<<_nsteps[_curPasses-1]
				<<", svidx="<<errorIndex;
			
			for (int k=1;k<model()->components().size();k++) {
				ModelComponent* comp = model()->components()[k];
				if (errorIndex>=comp->svOffset() && 
					errorIndex<comp->svOffset()+comp->numStateVar() ) {
					cerr<<" sv("<< errorIndex-comp->svOffset()<<")"
						<<" of "<<comp->componentName();
					break;
				}
			}
			cerr<<endl;
		}

		// Make sure the new number of steps is valid
		if (nstepsNext > _maxSteps) {
			
			// Cannot increase steps, so reduce H and try again
			// unless this would take us below the minStep value.
			// Changing H values is a drastic step taken only
			// when absolutely necessary and cannot be efficient.

			if (H/2<minTimeStep()) {
				break;		// Give up and use what we now have
			}

			H /= 2;			// Halve the time step attempted
			npass = 0;		// Start over from the first pass

			if (debugPerformance() ) {
				cerr<<"NeuronSolver "
					<<"T="<<tStart
					<<", err="<<_errorEst
					<<", H="<<H
					<<endl;
			}
		}
		else {
			// Otherwise, are we already at the maximum passes.
			// If so, increase the number of steps for each pass.
			if (_curPasses == _maxPasses) {

				// Shift the nsteps array by 1
				for (k=0;k<_curPasses-1;k++) {
					_nsteps[k] = _nsteps[k+1];
				}
				_nsteps[_curPasses-1] = nstepsNext;

				// We can avoid a lot of work by shifting
				// ending state vectors down to go with
				// the change in nsteps.

				Number* tempState = endingState[0];
				for (k=0;k<_curPasses-1;k++) {
					endingState[k]=endingState[k+1];
				}
				endingState[_curPasses-1] = tempState;

				// Since we have discarded a pass, decrement the count
				// Then allow the loop to continue with the number of
				// steps we just set in _nstep.
				npass--;
			}
			else {
				// Otherwise, increase the order by adding more passes.
				_curPasses++;
				_nsteps[_curPasses-1] = nstepsNext;

				// Allocate a corresponding state vector if needed
				if (endingState[_curPasses-1] == NULL) {
					endingState[_curPasses-1] = allocateNumVector();
				}

				// Allow the loop to continue with an expanded
				// number of passes. Work already done is preserved.
			}
		}

		// Rebuild the vector of time steps and extrap coeff
		buildHvect(H);
		buildExtrapolationCoeff();
	}
		
	// Extrapolate final state values.
	// Assumed first or second order extrapolation is determined
	// by the orderIndex entry corresponding to each value.
	currentTime(tStart+H);
	for (i=1;i<stateVectorSize;i++) {
		double stateValue = 0;
		for (k=0;k<_curPasses;k++) {
			stateValue += _stateExtrap[orderIndex[i]][k]*endingState[k][i];
		}
		modelState[i] = stateValue;
	}
	notifyOnStateVectorChanged();

	// Get derivatives for next step and do final notifications.
	computeDerivatives();
	notifyOnTimeStepEnded();

	// Check to see if the next step can be more efficient.
	// Allow a little cushion to avoid too frequent changes.
	if (_errorEst < safetyMargin()*errTol()) {

		// First, try to reduce the number of steps. Make a
		// rough estimate assuming at least h^2 accuracy to see
		// if this could work before making the change. If it 
		// looks feasible, cut the number of steps per pass in 
		// half, but only if this leaves sufficient steps.

		nstepsNext = (_nsteps[0]+1)/2;

		if (nstepsNext>=minSteps() && 4*_errorEst<errTol()) {

			// Show debug information if requested
			if (debugPerformance() ) {
				cerr<<"NeuronSolver "
					<<"T="<<tStart
					<<", err="<<_errorEst
					<<", ns["<<_curPasses-1<<"]="<<_nsteps[_curPasses-1]
					<<endl;
			}

			// Cut the number of steps in half 
			for (k=_curPasses-1;k>0;k--) {
				_nsteps[k] = _nsteps[k-1];
			}
			_nsteps[0] = nstepsNext;
		}

		// Otherwise, try reducing the order of the process by 
		// dropping the last h value. If this does not work out, 
		// the h value can be restored without much loss of effort.
		// At least two passes must remain for the extrapolation
		// process to work.
		else if (_curPasses>2) {
			_nsteps[_curPasses-1] = 0;	// keep array clean					
			_curPasses--;

			// Show debug information if requested
			if (debugPerformance() ) {
				cerr<<"NeuronSolver "
					<<"T="<<tStart
					<<", err="<<_errorEst
					<<", ns[0]="<<_nsteps[0]
					<<", passes="<<_curPasses
					<<endl;
			}
		}
	}

	// Free work arrays and Jacobian storage (if not saved)
	// It probably doesn't make much difference, but
	// work arrays are deleted in the opposite order from
	// their allocation.
	for (npass=_PerPassArraySize-1;npass>=0;npass--) {
		deleteNumVector(endingState[npass]);
	}
	delete[] compDeriv;
	delete[] orderIndex;
	delete[] startingDeriv;
	delete[] startingState;
	if (_alwaysRebuildJacobian) {
		_jacobian.resize(0,0);
	}
}



// ====================================================================
// ExternalSpikeRecorder class body
// ====================================================================



// Create a new instance
ExternalSpikeRecorder::ExternalSpikeRecorder(char* fn)
: ExternalRecorder(fn) {}

// Destroy this instance
ExternalSpikeRecorder::~ExternalSpikeRecorder() {}

// Write the id and last spike time in milliseconds
void ExternalSpikeRecorder::signalEvent(unsigned int classId, ModelComponent* mc)
{
	if (classId != ActionPotentialEvent::eventClassId()) {
		return;
	}

	SpikingNeuron* neuron = reinterpret_cast<SpikingNeuron*>(mc);

	// Write the output for a spike
	if (!_headerWritten) {
		fprintf(_file,"id,spike_time\n");
		_headerWritten = true;
	}
	fprintf(_file, "%d,%.12g\n",
		neuron->numericIdentifier(),
		neuron->spikeTime()/UOM::msec);
}



// ====================================================================
// ExternalDendriteSpikeRecorder class body
// ====================================================================



// Constructors and destructor
ExternalDendriteSpikeRecorder::ExternalDendriteSpikeRecorder(
	char*				fn,			// file name to write to
	Number				spikeVm)	// voltage threshold for spike
:	ExternalRecorder(fn)
{
	_spikeVm = spikeVm;
	initialize();
}

ExternalDendriteSpikeRecorder::~ExternalDendriteSpikeRecorder() 
{
	// Delete any allocated storage
	delete[] _prevVm;
}

// Set state to an initial condition
void ExternalDendriteSpikeRecorder::initialize()
{
	_neuron = NULL;
	_prevVm = NULL;
	_nDend = 0;
	_prevSomaVm = VMinForIndex;
	_somaSpikeTime = InfinitePast;
}

// Note when removed from a component
void ExternalDendriteSpikeRecorder::removedFrom(ModelComponent* mc)
{
	// If the component is the neuron monitored, reset
	// this probe to an initial state.
	if (mc==_neuron) {
		delete[] _prevVm;
		initialize();
	}
}

// Identify the compartment where spikes initiate
Compartment* ExternalDendriteSpikeRecorder::somaComp()
{
	return _neuron->somaComp();
}

// Report on spikes found in the compartments provided.
void ExternalDendriteSpikeRecorder::reportOn(
	ModelComponentVector&	comps,		// components to report on
	SimTime					t,			// time of probe
	int						id)			// id of component probed
{
	// Make sure file was opened before reporting
	if (_file == NULL) {
		FatalError("(ExternalDendriteSpikeRecorder::reportOn) "
			"output file is not open");
		return;
	}

	// See if the monitored neuron has not yet been found
	if (_neuron==NULL) {

		ModelComponentVectorIt	it;

		// This has to be done the hard way. First we
		// have to look through the components to find
		// the neuron. Then we can report on the dendrites
		// found there. This is necessary since dendrite
		// comparmtents do not know their dendrite number.
		for (it=comps.begin();it!=comps.end();it++) {

			if ((*it)->isNeuron() 
				&& (reinterpret_cast<Neuron*>(*it))->isMorphologicalNeuron()) {

				// We can stop now. Typically the neuron is
				// the first component so this is not so bad.
				_neuron = reinterpret_cast<MorphologicalNeuron*>(*it);
				break;
			}
		}
	}

	// Once the monitored neuron is found, report on it
	if (_neuron!=NULL) {

		// Write a header line if needed
		if (!_headerWritten) {
			writeColumnHeader();
			_headerWritten = true;
		}

		// Allocate and initialize work array if needed
		if (_prevVm==NULL) {
			_nDend = _neuron->numDendriteComp();
			_prevVm = new Number[_nDend+1];

			int i;
			for (i=0;i<=_nDend;i++) {
				_prevVm[i] = VMinForIndex;
			}
		}

		// Check for spikes in the soma
		if (_prevSomaVm<spikeVm() && somaComp()->Vm()>=spikeVm() ) {
			_somaSpikeTime = t;
		}
		_prevSomaVm = somaComp()->Vm();

		// Check for dendritic spikes and write them out
		writeSpikes(t);
	}
}

// Write a column header
void ExternalDendriteSpikeRecorder::writeColumnHeader()
{
	fprintf(_file,"neuron_id,spike_time,dend_nbr,spike_time_diff\n");
}

// Check for dendrite spikes and write them out
void ExternalDendriteSpikeRecorder::writeSpikes(SimTime t)
{
	using namespace UOM;

	double	tInMs = t/msec;
	double	dTInMs;
	int		i;

	// Get the difference between current time and the last cell spike.
	// If there is no such spike, use a large value rather than infinity.
	dTInMs = _somaSpikeTime!=InfinitePast 
		? (t - _somaSpikeTime)/msec
		: 9999e3;

	for (i=1; i<=_nDend; i++) {

		Compartment* dend = _neuron->dendriteComp(i);

		// Check for exceeding the voltage threshold
		if (dend->Vm()>=spikeVm() && _prevVm[i]<spikeVm() ) {

			fprintf(_file,"%d,%.12g,%d,%.8g\n",
				_neuron->numericIdentifier(),	// neuron identifier
				tInMs,							// time at end of step (ms)
				i,								// dendrite number
				dTInMs);						// time difference (ms)
		}

		// Save the voltage for next time
		_prevVm[i] = dend->Vm();
	}
}



// ====================================================================
// ExternalSynapseRecorder class body
// ====================================================================



// Create a new instance using a token identifier
ExternalSynapseRecorder::ExternalSynapseRecorder(char* fn, TokenId synType)
: ExternalRecorder(fn) 
{
	// Save the token identifier of the synapse type to report on
	_synapseType = synType;

	// Set a default minimum reporting interval
	minInterval( 100*UOM::msec );
}

// Create a new instance using a synapse name
ExternalSynapseRecorder::ExternalSynapseRecorder(char* fn, char* synName)
: ExternalRecorder(fn) 
{
	// Convert the synapse name to a token and save it
	if (synName==NULL) {
		_synapseType = NullTokenId;
	}
	else {
		_synapseType = token(synName);
	}

	// Set a default minimum reporting interval
	minInterval( 100*UOM::msec );
}

// Destroy this instance
ExternalSynapseRecorder::~ExternalSynapseRecorder() {}

// Report the current synapse weights
void ExternalSynapseRecorder::reportOn(
	ModelComponentVector&	comps,		// components to report on
	SimTime					t,			// time of probe
	int						id)			// id of component probed
{
	ModelComponentVectorIt		mcit;
	SynapticResponseVectorIt	srit;
	Synapse*					syn;

	// Make sure file was opened before reporting
	if (_file == NULL) {
		FatalError("(ExternalSynapseRecorder::reportOn) output file is not open");
	}

	// Make sure there is a synaptic type to report on
	if (synapseType() == NullTokenId)
		return;

	// Find any synaptic responses and get related weight(s).
	// Since not all synaptic responses are registered with the model
	// as components, get compartments and from these, the appropriate
	// synaptic responses based on token id.

	for (mcit=comps.begin();mcit!=comps.end();mcit++) {
		if ( (*mcit)->isCompartment() ) {

			// Locate the synaptic response within this compartment.
			// Go on to the next compartment if no matching response is found.
			Compartment*		pcomp = reinterpret_cast<Compartment*>(*mcit);
			SynapticResponse*	resp = pcomp->findSynapticResponse( synapseType(), false);
			if (resp== NULL)
				continue;		// move on to next component

			// Get the numeric identifier of the compartment
			int	compId = resp->container()->numericIdentifier();
			
				
			// Write a header line if needed. This must be delayed to the
			// first encounter to get the component names for various weights.
			if (!_headerWritten) {
				fprintf(_file,"post_id,time,pre_id,comp_id");

				// Write the header for a synaptic group
				if (resp->isSynapticGroup() ) {

					SynapticResponseVector& mbrs = resp->members();
					for (srit=mbrs.begin();srit!=mbrs.end();srit++) {
						fprintf(_file,",%s_weight",(*srit)->componentName() );
					}
					fprintf(_file,"\n");
				}
				else { // Otherwise write header for synaptic response
					fprintf(_file,",%s_weight\n",resp->componentName() );
				}
				_headerWritten = true;
			}

			// Go through the synapse list and report weights
			// The user should know how many weight values are reported per synapse
			// based on the type of the synaptic response being reported.
			for(syn=resp->synapses();syn!=NULL;syn=syn->nextPostsynaptic() ) {

				// Write identifiers for this line
				int preId = syn->presynapticProcess()->neuronId();
				fprintf(_file, "%d,%.12g,%d,%d", id, t/UOM::msec, preId, compId);

				// Write all weights for a synaptic group
				if (resp->isSynapticGroup() ) {
					SynapticResponseVector& mbrs = resp->members();
					for (srit=mbrs.begin(); srit!=mbrs.end(); srit++) {
						fprintf(_file, ",%g", double( (*srit)->synapseWeight(syn) ));
					}
					fprintf(_file,"\n");	
				}
				else { // Write weights for a synaptic response
					fprintf(_file, ",%g\n", double( resp->synapseWeight(syn) ));
				}
			}
		}
	}
}



// ====================================================================
// ExternalCompartmentRecorder class body
// ====================================================================



// Create a new instance
ExternalCompartmentRecorder::ExternalCompartmentRecorder(char* fn, char* mapFileName)
: ExternalStateRecorder(fn,mapFileName) {}

// Destroy this instance
ExternalCompartmentRecorder::~ExternalCompartmentRecorder() {}

// Write out a column header line for documentation
void ExternalCompartmentRecorder::writeColumnHeader(ModelComponentVector& comp)
{
	ModelComponentVectorIt		it;
	const char*					name;
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
	for (it=comp.begin();it!=comp.end();it++) {

		// Select just the compartments
		if ( !(*it)->isCompartment() ) {
			continue;
		}

		// Use the component name as a label
		name=(*it)->componentName();
		fprintf(_file,",%s",name);
		if (_mapFile != NULL) {
			fprintf(_mapFile,"%d\t%s\n",++colNbr,name);
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

void ExternalCompartmentRecorder::writeState(ModelComponentVector& comp, SimTime t, int id)
{
	ModelComponentVectorIt		it;

	// Write out identifier and time
	fprintf(_file, "%d,%.12g",id,t/UOM::msec);
	
	// Get state values from components
	for (it=comp.begin();it!=comp.end();it++) {

		// Select only compartments
		if ( !(*it)->isCompartment() ) {
			continue;
		}

		// Write out the compartment value (via subclass)
		writeValue(reinterpret_cast<Compartment*>(*it) );
	}

	// Write the end of line
	fprintf(_file,"\n");
}



// ====================================================================
// ExternalVoltageRecorder class body
// ====================================================================



// Create a new instance
ExternalVoltageRecorder::ExternalVoltageRecorder(char* fn, char* mapFileName)
: ExternalCompartmentRecorder(fn,mapFileName) {}

// Destroy this instance
ExternalVoltageRecorder::~ExternalVoltageRecorder() {}

// Write the compartment value
void ExternalVoltageRecorder::writeValue(Compartment* pcomp)
{
	fprintf(_file,",%8.4f",double( pcomp->stateValueConverted(0) ));
}



// ====================================================================
// ExternalCurrentRecorder class body
// ====================================================================



// Create a new instance
ExternalCurrentRecorder::ExternalCurrentRecorder(char* fn, char* mapFileName)
: ExternalCompartmentRecorder(fn,mapFileName) {}

// Destroy this instance
ExternalCurrentRecorder::~ExternalCurrentRecorder() {}

// Write the compartment value
void ExternalCurrentRecorder::writeValue(Compartment* pcomp)
{
	fprintf(_file,",%8.4f",double( pcomp->Itotal() / pcomp->ItotalUnits() ));
}



// ====================================================================
// ExternalIonCurrentRecorder class body
// ====================================================================



// Create a new instance
ExternalIonCurrentRecorder::ExternalIonCurrentRecorder(char* fn, char* mapFileName)
: ExternalStateRecorder(fn,mapFileName) {}

// Destroy this instance
ExternalIonCurrentRecorder::~ExternalIonCurrentRecorder() {}

// Write out a column header line for documentation
void ExternalIonCurrentRecorder::writeColumnHeader(ModelComponentVector& comp)
{
	ModelComponentVectorIt		it;
	IonChannel*					chan;
	const char*					name;
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
	for (it=comp.begin();it!=comp.end();it++) {

		// Select just the compartments
		if ( !(*it)->isIonChannel() ) {
			continue;
		}

		// Use the component name as a label.
		// Usually this will be preceeded by a compartment name.
		chan = reinterpret_cast<IonChannel*>(*it);
		name = chan->componentName();
		if (chan->container()==NULL) {
			fprintf(_file,",%s",name);
			if (_mapFile != NULL) {
				fprintf(_mapFile,"%d\t%s\n",++colNbr,name);
			}
		}
		else {
			fprintf(_file,",%s_%s",chan->container()->componentName(),name);
			if (_mapFile != NULL) {
				fprintf(_mapFile,"%d\t%s_%s\n",++colNbr,
					chan->container()->componentName(),name);
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

void ExternalIonCurrentRecorder::writeState(
	ModelComponentVector& comp, SimTime t, int id)
{
	ModelComponentVectorIt		it;
	IonChannel*					chan;
	double						current;

	// Write out identifier and time
	fprintf(_file, "%d,%.12g",id,t/UOM::msec);
	
	// Get state values from components
	for (it=comp.begin();it!=comp.end();it++) {

		// Select only compartments
		if ( !(*it)->isIonChannel() ) {
			continue;
		}

		// Write out the current value in appropriate units
		// with a comma delimiter separating the values
		chan = reinterpret_cast<IonChannel*>(*it);
		current = chan->Iion() / chan->IionUnits();
		fprintf(_file,",%g",current);
	}

	// Write the end of line
	fprintf(_file,"\n");
}



// ====================================================================
// ExternalCalciumRecorder class body
// ====================================================================



// Create a new instance
ExternalCalciumRecorder::ExternalCalciumRecorder(char* fn)
: ExternalStateRecorder(fn) {}

// Destroy this instance
ExternalCalciumRecorder::~ExternalCalciumRecorder() {}

// Write out a column header line for documentation
void ExternalCalciumRecorder::writeColumnHeader(ModelComponentVector& comp)
{
	ModelComponentVectorIt		it;

	const char*					poolName;
	const char*					compName;

	// Write standard label for id and simulation time
	fprintf(_file,"id,t");

	// Get state vector labels from components
	for (it=comp.begin();it!=comp.end();it++) {

		// Select just the calcium pools
		if ( (*it)->isCalciumPool() ) {

			// Make a label from the compartment and
			// calcium pool names.
			poolName=(*it)->componentName();
			compName=(reinterpret_cast<CalciumPool*>(*it))->container()->componentName();
			fprintf(_file,",%s_%s",compName,poolName);
		}
	}
	// write end of header line
	fprintf(_file,"\n");

}

// Write out the current state values preceeded by
// the component identifier and time. Time is converted
// to milliseconds. Other state values are converted
// based on model component units of measure values.

void ExternalCalciumRecorder::writeState(ModelComponentVector& comp, SimTime t, int id)
{
	ModelComponentVectorIt		it;

	// Write out identifier and time
	fprintf(_file, "%d,%.12g",id,t/UOM::msec);
	
	// Get state values from components
	for (it=comp.begin();it!=comp.end();it++) {

		// Select only calcium pools and write the
		// calcium concentration (assumed to be state 0)
		if ( (*it)->isCalciumPool() ) {
			// Write out the first state value (membrane voltage)
			// with a comma delimiter separating the values
			fprintf(_file,",%g",double((*it)->stateValueConverted(0)));
		}
	}

	// Write the end of line
	fprintf(_file,"\n");
}




// ====================================================================
// Neuron class body
// ====================================================================




// Create a new instance and hook to the model provided, if any.
// The ODESolver is found via the model and is allocated later
// if needed.
Neuron::Neuron(Model* m)
{
	// Connect with a model, creating one if needed.
	if (m==NULL) {
		_usesOwnModel = true;
		model(new Model);
	}
	else {
		_usesOwnModel = false;
		model(m);
	}

	// Set initial values 
	_allocatedSolver = NULL;
	_numericIdentifier = -1;;
}

// Destroy the instance and its associated model and solver,
// but only if the model was allocated by this instance.
Neuron::~Neuron() 
{
	// If this neuron created its own model,
	// quickly break the relationship with
	// all components.
	if (_usesOwnModel) {
		model()->clearComponents();
	}

	// Delete all owned components
	deleteAllCompartments();

	// If solver was allocated here, delete it
	if (_allocatedSolver != NULL) {
		delete _allocatedSolver;
		_allocatedSolver=NULL;
	}

	// Unhook relationship with the model, if any
	model(NULL);

}

// Associate with the model and pass along to compartments
void Neuron::model(Model* m)
{
	CompartmentVectorIt		it;

	// if there is a prior model allocated by this
	// instance, clear it out first along with the
	// associated solver.
	if (_model != NULL && _usesOwnModel) {
		_usesOwnModel = false;

		if (_allocatedSolver != NULL && _model->solver() == _allocatedSolver) {
			delete _allocatedSolver;
			_allocatedSolver = NULL;
		}
		delete _model;
		_model = NULL;
	}

	// Hook up with the new model, if any
	ModelComponent::model(m);

	// Pass new model along to subordinates
	for (it=_compartments.begin();it!=_compartments.end(); it++) {
		(*it)->model(m);
	}
}

// Accessor for ODE solver. If there is no current solver, a default
// solver is allocated as a lazy initialization.
ODESolver* Neuron::solver() 
{ 
	if (model()->solver() == NULL) {
		if (_allocatedSolver==NULL) {
			_allocatedSolver = defaultSolver();
		}
		_allocatedSolver->model( model() );
	}
	return model()->solver(); 
}

// Set ODE solver after disposing of any prior solvers
void Neuron::solver(ODESolver* newSolver) 
{ 
	// Dispose of any old allocated solver, but only if different
	if (newSolver != _allocatedSolver && _allocatedSolver!=NULL) {

		// Watch out for sequence error in setting controller and
		// new solver. Must clear old solver if already set up with
		// a controller. This catches an easy error of assuming
		// that controllers go with neuron rather than solvers.
		if (_allocatedSolver->hasController()) {
			FatalError("(Neuron::solver) Controller already assigned to current solver.");
		}

		delete _allocatedSolver;
		_allocatedSolver=NULL;
	}
		
	// Hook up model with new solver
	newSolver->model( model() );
}

// Create a default ODE solver for this neuron.
// This can be overridden if a different solver is needed.
ODESolver* Neuron::defaultSolver()
{
	return new AdaptiveSolver;
}

// Associate a probe object with this neuron
void Neuron::addProbe(Probe* pr)
{
	_probes.push_back(pr);
}

// Remove the probe association
void Neuron::removeProbe(Probe* pr)
{
	ProbeVectorIt		last;

	last=remove(_probes.begin(),_probes.end(),pr);
	_probes.resize(last - _probes.begin());
}

// Add a new compartment
void Neuron::addCompartment(Compartment* comp)
{
	// Add as a compartment held by the neuron
	_compartments.push_back(comp);
	comp->neuron(this);
	
	// Also add as a component of the model
	comp->model( model() );
}

// Remove a compartment
void Neuron::removeCompartment(Compartment* comp)
{
	CompartmentVectorIt		last;

	// Remove from the compartments of the neuron
	last=remove(_compartments.begin(),_compartments.end(),comp);
	_compartments.resize( last - _compartments.begin() );
	comp->neuron(NULL);

	// Also remove from the model
	comp->model(NULL);
}

// Delete all associated compartments
void Neuron::deleteAllCompartments()
{
	CompartmentVectorIt		it;
	CompartmentVector		temp;

	// Use a temp collection to avoid problems.
	// This also avoids problems of removing entries
	// underneath the iterator.
	swap(_compartments,temp);

	// Delete held compartments
	for (it=temp.begin(); it!=temp.end(); it++) {
		delete *it;
	}
}

// Sum the membrane area of all compartments
Number Neuron::membraneArea()
{
	CompartmentVectorIt		it;
	Number					ma = 0;

	for (it=_compartments.begin();it!=_compartments.end();it++) {
		ma += (*it)->membraneArea();
	}
	return ma;
}

// Do initialization on simulation start
void Neuron::simulationStarted()
{
	ProbeVectorIt			pit;

	// Notify probes
	for (pit=_probes.begin();pit!=_probes.end();pit++) {
		(*pit)->simulationStarted();
	}
}

// Do clean up on simulation ended
void Neuron::simulationEnded()
{
	ProbeVectorIt			pit;

	// Notify probes
	for (pit=_probes.begin();pit!=_probes.end();pit++) {
		(*pit)->simulationEnded();
	}
}

// When time step is finalized, notify affected parties
void Neuron::timeStepEnded()
{
	ProbeVectorIt			pit;

	// Notify probes
	for (pit=_probes.begin();pit!=_probes.end();pit++) {
		(*pit)->timeStepEnded( model(), numericIdentifier() );
	}
}

// Set each compartment's distance from a root compartment
void Neuron::setDistFromSoma(Compartment* root)
{
	CompartmentVector		peers;	// connected compartments
	CompartmentVectorIt		compIt;	// compartment iterator
	stack<Compartment*>		wip;	// work in progress
	Compartment*			next;	// next commpartment to work on
	Number					dist;	// distance from soma

	// A NULL value of root implies use of the default root
	if (root==NULL) {
		root = _compartments[0];
	}

	// Clear any current distances -- this gives a way to know
	// what nodes have been visited since only the root has a
	// distance of 0.
	for (compIt=_compartments.begin(); compIt!=_compartments.end(); compIt++) {
		(*compIt)->distFromSoma(0);
	}

	// Set the root distance to 0 and start expanding from there
	root->distFromSoma(0);
	wip.push(root);

	// Expand from the contents of the stack until all done
	while (!wip.empty() ) {
		// Get the next compartment to work on and its peers
		next=wip.top();
		wip.pop();
		peers=next->connectedCompartments();

		// Assign a distance to each peer and put it on the stack
		// The distance is the sum of half the lengths of the
		// two connected compartments. However, don't do any
		// compartment twice.
		for (compIt=peers.begin(); compIt!=peers.end(); compIt++) {
			if (*compIt!= root && (*compIt)->distFromSoma()==0) {
				dist = next->distFromSoma();
				dist += next->length()/2;
				dist += (*compIt)->length()/2;
				(*compIt)->distFromSoma(dist);
				wip.push(*compIt);
			}
		}
	}
}



// ====================================================================
// SpikingNeuron class body
// ====================================================================



// Create an instance and use the model provided, if any
SpikingNeuron::SpikingNeuron(Model* m)
:	Neuron(m)
{
	// Set initial values
	_spikeTime = InfinitePast;
	_axonProcess = new AxonProcess;
}

// Destroy this instance
SpikingNeuron::~SpikingNeuron()
{
	delete _axonProcess;
}

// Pass any changes in numeric identifier along to axon.
void SpikingNeuron::numericIdentifier(int id)
{
	_numericIdentifier = id;
	_axonProcess->neuronId(id);
}

// Signal that a spike occurred.
void SpikingNeuron::signalSpikeEvent(SimTime t)
{
	ProbeVectorIt			probeIt;
	SimTime					isi;

	// Compute inter spike interval and then 
	// set the time of the current spike
	isi = t - _spikeTime;
	_spikeTime = t;

	// Have the axon propagate the firing event
	axonProcess()->signalSpikeEvent(_spikeTime, isi, firingRate() );

	// Notify probes that a spike was fired
	for (probeIt=_probes.begin();probeIt!=_probes.end();probeIt++) {
		(*probeIt)->signalEvent(ActionPotentialEvent::eventClassId(),this);
	}
}



// ====================================================================
// VoltageTriggeredSpikingNeuron class body
// ====================================================================



// Create an instance and use the model provided, if any
VoltageTriggeredSpikingNeuron::VoltageTriggeredSpikingNeuron(Model* m)
:	SpikingNeuron(m)
{
	initialize();
}

// Destroy this instance
VoltageTriggeredSpikingNeuron::~VoltageTriggeredSpikingNeuron() {}


// Initialize the instance
void VoltageTriggeredSpikingNeuron::initialize()
{
	using namespace UOM;

	// Set initial state
	_previousFiringVm = VMinForIndex;
	_previousFiringVmDot = 0;
	_firingRateS1 = 0;
	_firingRateS2 = 0;
	_ExpHFRTau = 0;

	// Set default values for params
	_firingRateTau = 125*msec;
}

// Set the firing rate estimator time constant
void VoltageTriggeredSpikingNeuron::firingRateTau(SimTime tau)
{
	// Save the parameter and force recomputations 
	_firingRateTau = tau;
	_ExpHFRTau = 0;
}


// Return an estimate of the firing rate
Number VoltageTriggeredSpikingNeuron::firingRate()
{
	return _firingRateS2/firingRateTau();
}

// Update the state variables used to estimate firing rate.
void VoltageTriggeredSpikingNeuron::updateFiringRate(bool spikeNow)
{
	// The estimate of firing rate is done by using an
	// alpha function kernel scaled so that the integral
	// of the estimate over time is unity.

	// The equations being solved are:
	// ds1/dt = -s1/tau + spike(t)
	// ds2/dt = -s2/tau +s1/tau;
	//
	// where tau is the time constant and spike(t) is a unit
	// impulse when there is a spike at time t and zero otherwise.
	//
	// Since integral(t/tau*exp(-t/tau)*dt) from 0 to infinity is tau,
	// a final scaling by 1/tau is needed (see firingRateEst). This
	// could be avoided, but the present scheme keeps s1 and s2 as 
	// unitless quantities.

	const SimTime h = model()->timeStepSize();
	const SimTime frtau = firingRateTau();

	// See if the time step size has changed
	// or derived value must be computed.
	if (model()->stepSizeChanged() || _ExpHFRTau == 0) {
		_ExpHFRTau = exp(-h/frtau);
	}

	// Update the state by solving the ODE explicitly
	_firingRateS2 = _ExpHFRTau * (_firingRateS2 + h/frtau*_firingRateS1);
	_firingRateS1 = _ExpHFRTau * _firingRateS1;

	// Update S1 if there is a spike. For simplicity, the exact spike time
	// is not used -- the spike is treated as the end of the timestep.
	if (spikeNow) {
		_firingRateS1 += 1;
	}
}

// Test for firing event at time step end. Default test
// provided here is voltage exceeding a threshold with
// a constraint on minimum interspike interval.
void VoltageTriggeredSpikingNeuron::timeStepEnded()
{
	const Number	vth = firingThreshold();
	const Number	vm = firingVoltage();
	const Number	vmDot = firingDeriv();
	const Number	prevVm = _previousFiringVm;
	const Number	prevVmDot = _previousFiringVmDot;

	const SimTime	stepStart = timeStepStart();
	const SimTime	stepEnd = currentTime();
	const SimTime	h = stepEnd-stepStart;

	Number			peakVm = maxval(vm,prevVm);

	// Save the current voltage and deriv for next time.
	_previousFiringVm = vm;
	_previousFiringVmDot = vmDot;

	// Estimate a peak voltage in the interval.
	// This is done by fitting linear estimates
	// at the beginning and end of the step to
	// a common value in the middle. There might
	// not always be such a common value, in which
	// case, the maximum at the end points is used.
	if (prevVmDot>0 && vmDot<=0) {

		SimTime dt = (vm-h*vmDot-prevVm)/(prevVmDot-vmDot);
		if (dt>0 && dt<h) {
			peakVm = prevVm+dt*prevVmDot;
		}
	}
	
	// If we do not meet the criteria for firing, stop now
	if (peakVm<vth) {

		// Update firing rate state without spike
		updateFiringRate(false);

		// Let super class notify probes before we exit here
		Neuron::timeStepEnded();
		return;
	}

	// Estimate a firing time by examining different cases.
	// Note that we might already be above the threshold 
	// at the start of the step because of minISI restrictions.

	SimTime			tspike;

	if (prevVm>=vth) {
		// Already passed threshold at the start of the step.
		tspike = stepStart;
	}
	else if (vm>=vth) {
		// Passed vth during the step. Use a linear interpolation
		// to estimate the crossing time ignoring derivatives.
		tspike = stepStart+h*(vth-prevVm)/(vm-prevVm);
	}
	else {
		// The interval had a peak voltage greater than
		// voltages at the start or end of the step.
		// Use a linear interpolation from the beginning
		// of the step to find the crossing point.
		tspike = stepStart+(vth-prevVm)/prevVmDot;
	}

	// Check that the minimum spike interval is satisfied
	if (tspike<spikeTime()+minISI()) {

		// There is an ISI rule conflict in the spike time.
		// See if the allowed spike time falls in the time step.
		// If so, delay the spike until the allowed time. 
		// If not, there is no spike in this time step.
		if (spikeTime()+minISI()>=stepStart &&
			spikeTime()+minISI()<=stepEnd ) {
			tspike = spikeTime()+minISI();
		}
		else {
			// For this step, the spike is lost to the ISI constraint.
			// Update firing rate state and then let the superclass
			// notify probes before we exit here.
			updateFiringRate(false);
			Neuron::timeStepEnded();
			return;
		}
	}

	// Update firing rate including the current spike
	updateFiringRate(true);

	// Signal that a spike occurred, notifying both other neurons and probes.
	signalSpikeEvent(tspike);

	// Let subclasses know that a firing occurred and potentially make
	// adjustments to the final state (e.g. by changing Vm).
	postFiringActions(tspike);

	// Notify probes (via superclass) of state at the end of the step
	Neuron::timeStepEnded();
}



// ====================================================================
// MorphologicalNeuron class body
// ====================================================================



// Constructors and destructor
MorphologicalNeuron::MorphologicalNeuron(Model* m)
: VoltageTriggeredSpikingNeuron(m) 
{
	initialize();
}

MorphologicalNeuron::~MorphologicalNeuron() {}

// Minimally initialize the instance during construction
void MorphologicalNeuron::initialize()
{
	// Set default values
	_morphology = NULL;
	_morphologySize = 0;
	_numSomaComp = 0;
	_dendriteCompOffset = 0;
	_numDendriteComp = 0;
	_axonCompOffset = 0;
	_numAxonComp = 0;
	_numISComp = 0;

	// Set a typical orientation
	_orientationX = 0;
	_orientationY = 1;
	_orientationZ = 0;
}

// Access a compartment in the soma starting with number 1.
// First allocated compartment is assumed to be the soma root.
Compartment* MorphologicalNeuron::somaComp(int n)
{
	if (n<1 || n>_numSomaComp) {
		FatalError("(MorphologicalNeuron::axonComp) index out of range.");
	}

	return _compartments[n-1];
}

// Access a dendrite compartment starting with number 1.
Compartment* MorphologicalNeuron::dendriteComp(int n)
{
	if (n<1 || n>_numDendriteComp) {
		FatalError("(MorphologicalNeuron::dendriteComp) index out of range.");
	}

	return _compartments[n-1+_dendriteCompOffset];
}

// Return the associated morphology entry for a dendrite
MorphologyEntry* MorphologicalNeuron::dendriteMorphologyEntry(int n)
{
	// Check that n has a valid value
	if (n<1 || n>_numDendriteComp) {
		FatalError("(MorphologicalNeuron::dendriteMorhologyEntry) index out of range.");
	}
	return &(_morphology[n-1+dendriteMorphologyOffset()]);
}

// Access a dendrite compartment by branch id. The dendrite
// with the soma distance closest to the value provided is returned.
// If the branch is not found, NULL is returned.
Compartment* MorphologicalNeuron::dendriteCompByBranch(int branchId, Number dist)
{
	int					k;
	Number				foundDist;

	Compartment*		foundComp = NULL;
	Compartment*		pcomp;
	MorphologyEntry*	pmorph;	

	// Must do a search of the whole neuron to match branchId
	for (k=1;k<=numDendriteComp();k++) {
		pcomp = dendriteComp(k);
		pmorph = dendriteMorphologyEntry(k);

		// Look for a compartment on this branch with a center
		// nearer to dist than already found.
		if (pmorph->branch == branchId) {
			if (foundComp==NULL || fabs(pcomp->distFromSoma()-dist)<foundDist) {
				foundComp = pcomp;
				foundDist = fabs(pcomp->distFromSoma()-dist);
			}
		}
	}
	return foundComp;
}

// Return whether the dendrite indicated is a basal dendrite or not
bool MorphologicalNeuron::isBasalDendrite(int n)
{
	return dendriteMorphologyEntry(n)->type==MorphologyEntry::basalDendriteType;
}

// Access an axon compartment starting with number 1.
Compartment* MorphologicalNeuron::axonComp(int n)
{
	if (n<1 || n>_numAxonComp) {
		FatalError("(MorphologicalNeuron::axonComp) index out of range.");
	}

	return _compartments[n-1+_axonCompOffset];
}

// Access an initial segment compartment starting with number 1.
Compartment* MorphologicalNeuron::ISComp(int n)
{
	if (n<1 || n>_numISComp) {
		FatalError("(MorphologicalNeuron::ISComp) index out of range.");
	}

	return _compartments[n-1+_ISCompOffset];
}

// Return offset in morphology table of first dendrite entry
int MorphologicalNeuron::dendriteMorphologyOffset()
{
	int k;

	for (k=0;k<_morphologySize;k++) {
		switch (_morphology[k].type) {

		case MorphologyEntry::basalDendriteType:
		case MorphologyEntry::apicalDendriteType:
			return k;
		}
	}

	FatalError("(MorphologicalNeuron::dendriteMorphologyOffset) dendrite entry not found");
	return 0;	// keep compiler happy
}


// Set the morphology table
void MorphologicalNeuron::morphology(MorphologyEntry* mt)
{
	const int		maxSize = 0x10000; // 64k entries assumed for max
	int				n;

	// Save the morphology table starting address.
	// Note that a neuron shares the morphology table with other
	// instances and thus does not delete the array in a destructor.
	_morphology = mt;

	// Find morphology table size by looking for an end marker
	for (n=0;mt[n].type!=MorphologyEntry::endType; n++) {
		if (n>maxSize) {
			FatalError("(MorphologicalNeuron::morphology) table end marker not found.");
		}
	}
	_morphologySize = n;
}

// Allocate the soma as a unipotential ball with membrane area the
// same as the first morphology entry compartment.
void MorphologicalNeuron::createSoma(Number cm_sp, Number rm_sp)
{
	SphericalCompartment*	soma;
	Number					radius, area;

	// The soma is treated as an equipotential ball with a
	// radius chosen to preserve membrane area.
	area=2*Pi*_morphology[0].r*_morphology[0].len*UOM::micron_2;
	radius=sqrt(area/(4*Pi));

	soma = new SphericalCompartment(radius,cm_sp,rm_sp);
	soma->componentName("Soma");

	// Make the soma part of the cell model
	add(soma);
	_numSomaComp = 1;
}

// Electrically connect dendrite compartments using common default of
// soma as first compartment and morphology entry. Dendrites are added
// following existing compartments.
void MorphologicalNeuron::createDendrites(
			Number cm_sp,		// Cm_specific
			Number ri_sp,		// Ri_specific
			Number rm_sp,		// Rm_specific
			Number areaAdjust)	// Membrane area adjustment
{
	const int			dmo = dendriteMorphologyOffset();
	int					k;
	char				nameBuf[16];
	Compartment*		pcomp;
	
	// Set starting offset and initialize number for counting
	_dendriteCompOffset = _compartments.size();
	_numDendriteComp = 0;

	// Add dendrite compartments stopping if a non-dendrite entry is found
	for (	k=dmo;	
			k<_morphologySize &&
				(_morphology[k].type == MorphologyEntry::basalDendriteType ||
				 _morphology[k].type == MorphologyEntry::apicalDendriteType);
			k++) {
		
		// Create a new dendrite compartment
		pcomp = new Compartment(
			_morphology[k].r * UOM::micron,			// radius
			_morphology[k].len * UOM::micron,		// length
			cm_sp,ri_sp,rm_sp,areaAdjust);			// electrical parameters etc.

		// Assign a name to the compartment for reporting
		sprintf(nameBuf,"D%04d",k);
		pcomp->componentName(nameBuf);

		// Add the compartment to the neuron's collection
		add(pcomp);

		// Count the compartment added
		_numDendriteComp++;
	}

	// Make the connections for dendrites
	connectCompartments(
		dmo,						// starting entry in morphology table		
		k-1,						// ending entry
		_dendriteCompOffset,		// offset in compartments vector
		somaComp() );					// root compartment in tree

	// Assign a path distance from soma to all compartments
	setDistFromSoma( somaComp() );
}

// Connect compartments using default of a single compartment soma
void MorphologicalNeuron::connectDendrites()
{
	const int	dmo = dendriteMorphologyOffset();
	int			k;

	// Find out how many dendrite entries there are in the morphology table.
	// All dendrite entries must occur continuously for this to work.
	k=dmo;
	while ( k<_morphologySize &&
			(_morphology[k].type == MorphologyEntry::basalDendriteType ||
			 _morphology[k].type == MorphologyEntry::apicalDendriteType) ) {
		k++;
	}

	// Make the connections for dendrites
	connectCompartments(
		dmo,						// starting entry in morphology table		
		k-1,						// ending entry
		_dendriteCompOffset,		// offset in compartments vector
		somaComp() );				// root compartment in tree
}

// Electrically connect compartments based on a morphology table.
void MorphologicalNeuron::connectCompartments(
	int					fromIdx,	// starting morphology index
	int					toIdx,		// ending index
	int					offset,		// starting offset into compartments vector
	Compartment*		root)		// root compartment to connect to
{
	// Working variables
	int						k,n,p;				// misc indexes
	Compartment*			pcomp;				// pointer to current neurite
	CompartmentVector*		children;			// children of a neurite
	ElectricalJunction*		junct;				// electrical junction to be added

	// If toIdx is -1, then use the remainder of the table
	if (toIdx == -1) {
		toIdx = _morphologySize-1;
	}

	// Sanity check that we are not going to go outside the compartments vector
	if (toIdx-fromIdx+offset>=_compartments.size() ) {
		FatalError("(MorphologicalNeuron::connectCompartments) "
			"Compartments vector size exceeded");
	}

	// Allocate work arrays cross referencing parents and children
	children = new CompartmentVector[toIdx-fromIdx+1];

	// Scan morphology table for parent-child relationships
	for (k=fromIdx; k<=toIdx; k++) {

		// Keep track of parent-child relationships in tree structure
		// Parent=0 indicates connection with the tree root.
		// Parentage out of the current from-to range is ignored.
		if (_morphology[k].parent!=0) {
			p = _morphology[k].parent;
			if (fromIdx<=p && p<=toIdx) {
				children[p-fromIdx].push_back(_compartments[k-1+offset]);
			}
		}		
	}

	// Make electrical connections between the compartments to
	// form the dendritic tree.
	for (k=fromIdx; k<=toIdx; k++) {

		pcomp=_compartments[k-fromIdx+offset];
		p = _morphology[k].parent;

		// If this compartment has no parent, connect it with the tree root
		// An assumption here is that the root is an equipotential ball
		// and thus multiple couplings have the same effect as a junction.
		if (p==0) {
			new ElectricalCoupling(root,pcomp);
		}

		// If there is only one child in the tree, use a simple coupling
		if (children[k-fromIdx].size()==1) {
			new ElectricalCoupling(pcomp,children[k-fromIdx][0]);
		}

		// If there are multiple children, create an n-way junction
		else if (children[k-fromIdx].size()>1) {
			junct = new ElectricalJunction(pcomp);
			for (n=0;n<children[k-fromIdx].size();n++) {
				junct->add(children[k-fromIdx][n]);
			}
		}

		// Otherwise, this compartment is a terminal node -- do nothing
	}

	// Delete work array(s)
	delete[] children;
}

// Set ion channel conductance multiplier in soma and all dendrites
// for the identified ion channel to be the value provided.
// mustMatch esures that either the soma or at least one compartment
// constains the indicated channel.
void MorphologicalNeuron::setGModulator(
	TokenId			chanId,	
	Number			value,
	bool			mustMatch)
{
	IonChannel*		chan;
	int				k;
	bool			matched = false;

	static TokenId	gModulatorId = token("gModulator");

	// Apply to soma
	chan = somaComp()->findIonChannel(chanId, false);
	if (chan!=NULL) {
		chan->setModParam(gModulatorId, value );
		matched = true;
	}

	// Apply to initial segment
	for (k=1;k<=numISComp();k++) {

		// Look for channel and set modulation if found
		chan = ISComp(k)->findIonChannel(chanId, false);
		if (chan!=NULL) {
			chan->setModParam(gModulatorId, value );
			matched = true;
		}
	}

	// Apply to axon
	for (k=1;k<=numAxonComp();k++) {

		// Look for channel and set modulation if found
		chan = axonComp(k)->findIonChannel(chanId, false);
		if (chan!=NULL) {
			chan->setModParam(gModulatorId, value );
			matched = true;
		}
	}

	// Apply to dendrites
	for (k=1;k<=numDendriteComp();k++) {

		// Look for channel and set modulation if found
		chan = dendriteComp(k)->findIonChannel(chanId, false);
		if (chan!=NULL) {
			chan->setModParam(gModulatorId, value );
			matched = true;
		}
	}

	// Check to see if the channel was found somewhere
	// and thus is at least a plausibly correct id.
	if (mustMatch && !matched) {
		FatalError("(MorphologicalNeuron::setGModulator) "
			"ChanId was not matched in any compartment.");
	}
}

// Set neuromodulation for channels using a Michaelis-Mentor formula:
// g_as_mod = g_normal * (1+a*X/(X+Kd)) 
void MorphologicalNeuron::setMichaelisMentenMod(
	TokenId			chanId,		// token id of channel component
	Number			X,			// Concentration
	Number			a,			// a in MM formula
	Number			Kd,			// Kd in above formula
	bool			mustMatch)	// must be matched somewhere
{
	setGModulator(chanId, 1+a*X/(X+Kd), mustMatch );
}



// ====================================================================
// Compartment class body
// ====================================================================



// Create a compartment
Compartment::Compartment(bool doInit)
{
	if (doInit) {
		setDefaults();
		initialize();
	}
}
 
Compartment::Compartment(
	Number r,				// radius
	Number len,				// length
	Number cm_sp,			// Cm_specific
	Number ri_sp,			// Ri_specific
	Number rm_sp,			// Rm_specific
	Number areaAdjustment)	// Membrane area adjustment factor
{
	// Set parameter values
	_radius = r;
	_length = len;
	_areaAdjustment = areaAdjustment;
	_Cm_specific = cm_sp;
	_Ri_specific = ri_sp;
	_Rm_specific = rm_sp;

	// Initialize this instance
	initialize();
}

// Destroy a compartment
Compartment::~Compartment() 
{
	// Delete constituents
	deleteAllCalciumPools();
	deleteAllIonChannels();
	deleteAllCouplings();

	// Unhook from the neuron
	if (neuron()!=NULL) {
		neuron()->removeCompartment(this);
	}

	// Delete the name string, if any
	if (_name != NULL) {
		delete[] _name;
	}
}

// Get the name string for the component
const char* Compartment::componentName()
{
	return _name == NULL ? "Compartment" : _name;
}

// Set the name string for the component
void Compartment::componentName(char* nameStr)
{
	// Copy the string 
	_name = new char[strlen(nameStr)+1];
	strcpy(_name,nameStr);
}


// Set default value for parameters
void Compartment::setDefaults()
{
	using namespace UOM;

	// A somewhat arbitrary list of defaults.
	// Your mileage may vary, but it's a start.
	_Cm_specific=1*microF/cm_2;
	_Rm_specific=50000*ohm*cm_2;
	_Ri_specific=100*ohm*cm;

	// Obviously arbitrary, but in some
	// models, only relative sizes matter and
	// we need some values as placeholders.
	_radius=10*micron;
	_length=10*micron;
	_areaAdjustment = 1;
}

// Initialize the object following parameter setting.
void Compartment::initialize()
{
	// Set NULL neuron at the outset
	_neuron = NULL;
	_numericIdentifier = -1;

	// Clear the default name
	_name = NULL;

	// Set default leak potential and initial potential
	_Vleak = -65*UOM::mV;
	_Vinit = -65*UOM::mV;
	_Vclamp = 0;

	// Set injected current to zero by default
	// and indicate that no voltage clamp is in effect
	_Iinjected = 0;
	_isVoltageClamped = false;
	_forceJacobianRebuild = false;

	// Other defaults
	_distFromSoma = 0;
	_jacobianIndex = 0;
	_gChanTotal = 0;
	_Itotal = 0;
	_Vm=0;
	_VmRem=0;

	// Adjust for the compartment size
	updateParams();
}

// Set up the compartment at the start of simulation
void Compartment::simulationStarted()
{
	// Set Vm based on Vinitial. This is done here
	// to allow channels to access Vm during
	// their setInitialState processing.
	Vm( _isVoltageClamped ? _Vclamp : Vinit() );

	// Initially set jacobian index to the svOffset
	_jacobianIndex = svOffset();
}

// Set initial states. This is somewhat redundant since
// the same thing was done in simulationStarted.
void Compartment::setInitialState()
{
	Vm( _isVoltageClamped ? _Vclamp : Vinit() );
}

// Set model and pass along to ion channels
void Compartment::model(Model* m)
{
	IonChannelVectorIt		itchan;
	CalciumPoolVectorIt		itpool;

	ModelComponent::model(m);

	// Pass model to ion channels
	for (itchan=_ionChannels.begin();itchan!=_ionChannels.end(); itchan++) {
		(*itchan)->model(m);
	}

	// Pass model to calcium pools
	for (itpool=_calciumPools.begin();itpool!=_calciumPools.end(); itpool++) {
		(*itpool)->model(m);
	}
}

// Set weight values for membrane potential
void Compartment::setWeightValues()
{
	weightValue(0) = 1.0/(10*UOM::mV);
}

// Set membrane voltage
void Compartment::Vm(Number v)
{
	_Vm=v;
	_VmIndex=VTableIndex(v);
	_VmRem=VTableRem(v, _VmIndex);
	stateValue(0)=v;
}

// Set radius
void Compartment::radius(Number r)
{
	_radius = r; 
	updateParams(); 
}

// Set length
void Compartment::length(Number r)
{
	_length = r; 
	updateParams(); 
}

// Set spine factor
void Compartment::areaAdjustment(Number x) 
{
	_areaAdjustment = x;
	updateParams();
}

// Set specific capacitance
void Compartment::Cm_specific(Number c) 
{
	_Cm_specific = c; 
	updateParams(); 
}

// Set specific membrane resistivity
void Compartment::Rm_specific(Number r) 
{
	_Rm_specific = r; 
	updateParams(); 
}

// Set specific internal resistivity
void Compartment::Ri_specific(Number r) 
{
	_Ri_specific = r; 
	updateParams(); 
}

// Set a voltage clamp to a given voltage
void Compartment::setVoltageClamp(Number vm)
{
	// Indicate the voltage clamp is in effect
	_isVoltageClamped = true;
	_Vclamp = vm;

	// Because voltage clamps affect the Jacobian,
	// entries related to this compartment are
	// built anew when next needed by the solver.
	_forceJacobianRebuild = true;

	// If simulation is started, set the new state
	if (model()!=NULL && model()->stateVector()!=NULL) {
		Vm( vm );
	}
}

// Clear a voltage clamp
void Compartment::clearVoltageClamp()
{
	_isVoltageClamped = false;
	_forceJacobianRebuild = true;
}

// Calculate the cross section
Number Compartment::crossSectionArea()
{
	Number r = radius();

	// calculate cross sectional area of a cylinder
	return Pi*r*r;
}

// Calculate the compartmental volume
Number Compartment::volume()
{
	return crossSectionArea()*length();
}

// Calculate the volume of a subshell of a given depth
Number Compartment::subshellVolume(Number depth)
{
	Number	r = radius();
	Number	d = depth<r ? depth : r;

	return Pi*(r*r - (r-d)*(r-d))*length(); 
}

// Calculate the area of a shell surrounding the compartment
// with an additional radius over that of the compartment.
Number Compartment::extraShellArea(Number dist)
{
	return 2*Pi*length()*(radius() + dist);
}

// Calculate absolute RC values from
// specific ones assuming a cylinder.
void Compartment::updateParams() 
{
	// Calculate membrane area assuming a cylinder
	// with an adjustment to allow for spines etc, which increase
	// membrane area but neither compartment volume or cross section.
	_membraneArea = 2*Pi*radius()*length()*areaAdjustment();

	// Calculate absolute values from
	// corresponding specific values.
	_Cm=Cm_specific()*membraneArea();
	_Rm=Rm_specific()/membraneArea();
	_Ri=Ri_specific()*length()/crossSectionArea();

	// Pass along updates to associated dependents
	propagateUpdates();
}

// Pass update event along to dependents
void Compartment::propagateUpdates()
{
	IonChannelVectorIt ionIt;
	ElectricalCouplingVectorIt ecIt;
	CalciumPoolVectorIt capIt;

	// Propagate changes to any ion channels
	for (ionIt=_ionChannels.begin(); ionIt!=_ionChannels.end(); ionIt++) {
		(*ionIt)->updateParams();
	}

	// Propagate changes to any couplings
	for (ecIt=_couplings.begin();ecIt!=_couplings.end();ecIt++) {
		(*ecIt)->updateParams();
	}

	// Propagate changes to any calcium pools
	for (capIt=_calciumPools.begin(); capIt!=_calciumPools.end(); capIt++) {
		(*capIt)->updateParams();
	}
}

// Add an ion channel to the compartment
void Compartment::addIonChannel(IonChannel* ic)
{
	// add the ion channel here
	_ionChannels.push_back(ic);
	ic->container(this);

	// also add the ion channel to the model
	ic->model( model() );
}

// Remove an ion channel from this compartment
void Compartment::removeIonChannel(IonChannel *ic)
{
	IonChannelVectorIt		last;

	// Remove the ion channel here (if present)
	last=std::remove(_ionChannels.begin(),_ionChannels.end(),ic);
	_ionChannels.resize( last - _ionChannels.begin() );
	ic->container(NULL);

	// Clear the model in the ion channel
	ic->model(NULL);
}

// Add an electrical coupling to the compartment
void Compartment::addCoupling(ElectricalCoupling* ec)
{
	_couplings.push_back(ec);
}

// Remove an electrical coupling from this compartment
void Compartment::removeCoupling(ElectricalCoupling *ec)
{
	ElectricalCouplingVectorIt		last;

	// remove the coupoing here
	last=std::remove(_couplings.begin(),_couplings.end(),ec);
	_couplings.resize( last - _couplings.begin() );
}

// Add an calcium pool to the compartment
void Compartment::addCalciumPool(CalciumPool* pool)
{
	// add the ion channel here
	_calciumPools.push_back(pool);
	pool->container(this);

	// also add the calcium pool to the model
	pool->model( model() );
}

// Return all compartments coupled with this one
CompartmentVector Compartment::connectedCompartments()
{
	CompartmentVector			connected;
	CompartmentVector			ecNodes;
	CompartmentVectorIt			compIt;
	ElectricalCouplingVectorIt	ecIt;

	// Go through all couplings -- each peer should be in only one
	for (ecIt=_couplings.begin(); ecIt!=_couplings.end(); ecIt++) {
		ecNodes = (*ecIt)->nodes();
		for (compIt=ecNodes.begin(); compIt!=ecNodes.end(); compIt++) {
			// Add every node but this one to the connected list
			if ( (*compIt)!=this ) {
				connected.push_back(*compIt);
			}
		}
	}

	return connected;
}

// Remove an calcium pool from this compartment
void Compartment::removeCalciumPool(CalciumPool *pool)
{
	CalciumPoolVectorIt		last;

	// Remove the ion channel here
	last=std::remove(_calciumPools.begin(),_calciumPools.end(),pool);
	_calciumPools.resize( last - _calciumPools.begin() );
	pool->container(NULL);

	// Clear the model in the ion channel
	pool->model(NULL);
}

// Delete held objects -- ion channels
void Compartment::deleteAllIonChannels()
{
	IonChannelVectorIt		it;
	IonChannelVector		temp;

	// Use a temporary vector to avoid
	// removing underneath the iterator
	// and also avoid extra remove overhead.
	swap(_ionChannels, temp);
	
	for (it=temp.begin(); it!=temp.end(); it++) {
		(*it)->container(NULL);
		delete *it;
	}
}

// Delete held object -- couplings
// Note that any peer can initiate the delete
void Compartment::deleteAllCouplings()
{
	ElectricalCouplingVectorIt	it;
	ElectricalCouplingVector	temp;

	// Use a temporary vector to avoid
	// removing underneath the iterator
	// and also avoid extra remove overhead.
	swap(_couplings, temp);

	// Now delete the coupling objects themselves
	for (it=temp.begin();it!=temp.end();it++) {
		delete *it;
	}
}

// Delete held objects -- calcium channels
void Compartment::deleteAllCalciumPools()
{
	CalciumPoolVectorIt		it;
	CalciumPoolVector		temp;

	// Use a temporary vector to avoid
	// removing underneath the iterator
	// and also avoid extra remove overhead.
	swap(_calciumPools, temp);
	
	for (it=temp.begin(); it!=temp.end(); it++) {
		(*it)->container(NULL);
		delete *it;
	}
}

// Get all calcium channels that are calcium sources and
// add them to the vector of channels provided.
IonChannelVector Compartment::getCalciumChannels()
{
	IonChannelVector		chan;
	IonChannelVectorIt		it;

	for (it=_ionChannels.begin();it!=_ionChannels.end();it++) {
		if ((*it)->isCalciumChannel() ) {
			chan.push_back(*it);
		}
	}
	return chan;
}

// Get an ion channel by matching on the component name.
// Only the first matching channel is returned.
IonChannel* Compartment::findIonChannel(char* channelName, bool mustMatch)
{
	return findIonChannel(token(channelName),mustMatch);
}

// Get an ion channel by matching on the component name.
// Only the first matching channel is returned.
IonChannel* Compartment::findIonChannel(TokenId channelId, bool mustMatch)
{
	IonChannelVectorIt	it;
	IonChannel*			found;

	// Look for a matching channel
	for (it=_ionChannels.begin(); it!=_ionChannels.end(); it++) {
		found = (*it)->findIonChannel(channelId);
		if (found!=NULL) {
			return found;
		}
	}

	// Handle the no match case
	if (mustMatch) {
		FatalError("(Compartment::findIonChannel) no matching channel found.");
	}
	return NULL;

}

// Find a synaptic response by matching on component name
SynapticResponse* Compartment::findSynapticResponse(char* channelName, bool mustMatch)
{
	return findSynapticResponse(token(channelName),mustMatch);
}


// Find a synaptic response by matching on component token id
SynapticResponse* Compartment::findSynapticResponse(TokenId channelId, bool mustMatch)
{
	IonChannel*		found = findIonChannel(channelId, false);

	// End now if nothing was found and no error occurred
	if (found==NULL) {
		if (mustMatch) {
			FatalError("(Compartment::findSynapticResponse) "
				"Synaptic response was not found.");
		}		
		return NULL;
	}

	// Make sure this is a synaptic response
	if (!found->isSynapticResponse()) {
		if (mustMatch) {
			FatalError("(Compartment::findSynapticResponse) "
				"Channel found is not a synaptic response.");
		}
		else {
			return NULL;
		}
	}
	return reinterpret_cast<SynapticResponse*>(found);
}

// Create a synapse give its component name
Synapse* Compartment::createSynapse(
	char*					name,			// synapse component name
	AxonProcess*			axon,			// presynaptic axon process
	Number					wght,			// initial weight value
	Number					dist)			// axonal distance
{
	return createSynapse(token(name),axon,wght,dist);
}

// Create a synapse give its component token id
Synapse* Compartment::createSynapse(
	TokenId					id,				// synapse component token id
	AxonProcess*			axon,			// presynaptic axon process
	Number					wght,			// initial weight value
	Number					dist)			// axonal distance
{
	SynapticResponse*		resp;
	
	resp = findSynapticResponse(id);

	return resp->createSynapse(axon,wght,dist);
}
								

// Cache voltage and its index for faster access
void Compartment::stateVectorChanged()
{
	_Vm=stateValue(0);
	_VmIndex=VTableIndex(_Vm);
	_VmRem=VTableRem(_Vm,_VmIndex);
}

// Compute the total current from all channels
// Leave the total conductance in _gIonTotal
double Compartment::Ichan()
{
	IonChannelVectorIt				it;
	double							isum = 0;

	Number							Ichan;
	Number							Gchan;
	
	// Clear the total conductance
	_gChanTotal = 0;

	// Get the current flow from each channel
	for (it=_ionChannels.begin(); it !=_ionChannels.end();it++) {
		(*it)->condAndIion(Gchan,Ichan);
		isum += Ichan;
		_gChanTotal += Gchan;
	}
	
	return isum;
}

// Compute the current from all electrical couples
double Compartment::Icouple()
{
	// Handle the next most common case via explicit code
	if (_couplings.size()==2) {
		return _couplings[0]->Iec(this)+_couplings[1]->Iec(this);
	}

	// Otherwise, use a loop to get all coupling currents

	ElectricalCouplingVectorIt		icoup;
	double							isum = 0;

	// Get current flow via each adjacent compartment
	for (icoup = _couplings.begin(); icoup!=_couplings.end(); icoup++) {
		isum += (*icoup)->Iec(this);
	}

	return isum;
}

// Compute derivatives after the fact. This is just a hook
// in case there is some special action needed to avoid
// redoing something expensive or double counting.
void Compartment::recomputeDerivatives()
{
	computeDerivatives();
}

// Compute derivatives of local state variables.
// This is a part of the performance path and some
// levels of indirection have been removed.
void Compartment::computeDerivatives()
{
	// Start with the leak current
	_Itotal = Ileak();

	// Add the current flow from each channel
	_Itotal += Ichan();

	// Add the current from adjacent compartments
	_Itotal += Icouple();

	// Adjust for injected currents (inward negative sign convention)
	_Itotal -= Iinjected();


	// Check for a voltage clamp, in which case
	// the derivative is always zero. The above logic
	// is allowed to run anyway to permit side-effects.
	if (isVoltageClamped() ) {
		derivValue(0) = 0;
	}
	else {
		// Actually compute the derivative
		derivValue(0)= -_Itotal / _Cm;
	}
}

// Build the Jacobian row associated with this compartment
// using the Jacobian matrix provided in J.
void Compartment::buildJacobian(SparseMatrix& J)
{
	int			myIndex = jacobianIndex();
	int			nodeIndex;

	double		gCouple,gCoupleTotal;

	ElectricalCouplingVectorIt	coupleIt;

	CompartmentVector			nodes;
	CompartmentVectorIt			nodeIt;

	// Get each coupling/junction conductance
	// and build the associated Jacobian entry.
	// Since membrance conductances are assumed to be constant,
	// we can save the computation for this node with itself
	// for later use in the update.
	gCoupleTotal = 0;

	for (coupleIt=_couplings.begin();coupleIt!=_couplings.end();coupleIt++) {
		
		nodes = (*coupleIt)->nodes();
		for (nodeIt=nodes.begin();nodeIt!=nodes.end();nodeIt++) {
		
			// Make the jacobian entry for each connected
			// node as it is found (there should be no duplicates).
			// Total the conductance found and include it the
			// diagonal term to offset, e.g. cm*dv/dt = g12*(V2-V1) 
			// gives g12/cm for J[1][2] and -g12/cm for J[1][1].
			if (*nodeIt != this) {

				// Make the Jacobian entry. Note change of signs
				// corresponding to sign of derivative value.
				
				gCouple = (*coupleIt)->gForPair(this,*nodeIt);
				gCoupleTotal += gCouple;

				nodeIndex = (*nodeIt)->jacobianIndex();

				// If a voltage clamp is in effect, the result
				// is still a zero entry. Otherwise, set to
				// the partial derivative value.
				if (isVoltageClamped() ) {
					J[myIndex][nodeIndex] = 0;
				}
				else {
					J[myIndex][nodeIndex] = gCouple/Cm();
				}
			}
		}
	}

	// Can now make the diagonal entry for this compartment
	// and then we are done with this row. If there is a voltage
	// clamp in effect, this term is always zero.
	if (isVoltageClamped() ) {
		J[myIndex][myIndex] = 0;
	}
	else {
		J[myIndex][myIndex] = -(gLeak()+_gChanTotal+gCoupleTotal)/Cm();
	}

	// Clear any pending request to rebuild the Jacobian entry
	_forceJacobianRebuild = false;
}

// Update the Jacobian diagonal entry for this compartment.
// Connections are assumed constant and do not need to be updated.
// Compute derivatives must be done before this to update
// conductance totals.
void Compartment::updateJacobian(SparseMatrix& J)
{
	// If a Jacobian rebuild is required, do that rather than
	// just update the diagonal term.
	if (_forceJacobianRebuild) {
		buildJacobian(J);
		return;
	}

	int			myIndex = jacobianIndex();
	double		gCoupleTotal = 0;

	ElectricalCouplingVectorIt		coupleIt;

	for (coupleIt=_couplings.begin();coupleIt!=_couplings.end();coupleIt++) {
		gCoupleTotal += (*coupleIt)->gForPair(this,this);
	}

	// Can now make the diagonal entry for this compartment
	// and then we are done with this row. If there is a voltage
	// clamp in effect, this term is always zero.
	if (isVoltageClamped() ) {
		J[myIndex][myIndex] = 0;
	}
	else {
		J[myIndex][myIndex] = -(gLeak()+_gChanTotal+gCoupleTotal)/Cm();
	}
}

// Add to components to probe when reporting on
// this compartment. By default this includes
// ion channels and calcium pools as well as
// the compartment itself.
void Compartment::addToComponentsToProbe(ModelComponentVector& comps)
{
	int	i;
	
	// Put the compartment itself into the vector
	comps.push_back(this);

	// Put in ion channels and calcium pools.
	for (i=0;i<_ionChannels.size();i++) {
		_ionChannels[i]->addToComponentsToProbe(comps);
	}
	for (i=0;i<_calciumPools.size();i++) {
		_calciumPools[i]->addToComponentsToProbe(comps);
	}
}



// ====================================================================
// SphericalCompartment class body
// ====================================================================



// Create a compartment
SphericalCompartment::SphericalCompartment(bool doInit)
{
	if (doInit) {
		setDefaults();
		initialize();
	}
}
 
SphericalCompartment::SphericalCompartment(
	Number r,				// radius
	Number cm_sp,			// Cm_specific
	Number rm_sp,			// Rm_specific
	Number areaAdjustment)	// Membrane area adjustment factor
{
	// Set parameter values
	_radius = r;
	_areaAdjustment = areaAdjustment;
	_Cm_specific = cm_sp;
	_Rm_specific = rm_sp;

	// Initialize this instance
	initialize();
}

// Destroy a compartment
SphericalCompartment::~SphericalCompartment() {}

// Calculate the compartmental volume
Number SphericalCompartment::volume()
{
	Number	r = radius();

	return 4.0/3.0*Pi*r*r*r;
}

// Calculate the volume of a subshell of a given depth
Number SphericalCompartment::subshellVolume(Number depth)
{
	Number	r = radius();
	Number	d = depth<r ? depth : r;

	return 4.0/3.0*Pi*(r*r*r - (r-d)*(r-d)*(r-d)); 
}

// Calculate the area of a shell surrounding the compartment.
Number SphericalCompartment::extraShellArea(Number dist)
{
	Number r = radius() + dist;

	return 4.0*Pi*r*r;
}

// Calculate absolute RC values from
// specific ones assuming a cylinder.
void SphericalCompartment::updateParams() 
{
	// Set values that are not used by spheres
	_Ri = 0;
	_length = 0;

	// Calculate membrane area assuming a sphere
	// with an adjustment to allow for spines etc, which increase
	// membrane area but not compartment volume.
	_membraneArea = 4*Pi*radius()*radius()*areaAdjustment();

	// Calculate absolute values from
	// corresponding specific values.
	_Cm=Cm_specific()*membraneArea();
	_Rm=Rm_specific()/membraneArea();

	// Pass along updates to associated dependents
	propagateUpdates();
}

// Indicate error if unsupported function is used
void SphericalCompartment::Ri_specific(Number x)
{
	FatalError("(SphericalCompartment::Ri_specific) Function not supported for this class");
}

// Indicate error if unsupported function is used
void SphericalCompartment::length(Number x)
{
	FatalError("(SphericalCompartment::length) Function not supported for this class");
}

// Indicate error if unsupported function is used
Number SphericalCompartment::crossSectionArea()
{
	FatalError("(SphericalCompartment::crossSectionArea) Function not supported for this class");
	return 0;
}



// ====================================================================
// AxonProcess class body
// ====================================================================



// Construct a new instance
AxonProcess::AxonProcess(int id) 
{
	// Set the neuron owning this axon
	_neuronId = id;

	// Initialize variables to default values
	_synapses = NULL;
	_propRate = 0.5*UOM::meter/UOM::sec;
}

// Destroy this instance
AxonProcess::~AxonProcess() 
{
	Synapse*	syn;

	// Clear the list of synapses associated with this axon
	syn=_synapses;
	while (syn!=NULL) {
		syn = syn->clearPresynaptic();
	}
}

// Add an associated synapse
void AxonProcess::addSynapse(Synapse* syn)
{
	// Add to the beginning of the list
	syn->addBeforeInPresynapticList(_synapses);
	_synapses = syn;
}

// Remove an associated synapse
void AxonProcess::removeSynapse(Synapse* syn)
{
	Synapse*		next;

	// See if the synapse to be removed is the first one
	if (_synapses==syn) {
		_synapses = syn->nextPresynaptic();
	}
	else {
		// Scan the list looking for the synapse before the one to be removed
		next = _synapses;
		while (next!=NULL && next->nextPresynaptic()!=syn) {
			next=next->nextPresynaptic();
		}
		if (next!=NULL) {
			next->removeNextFromPresynapticList();
		}
		else {
			FatalError("(AxonProcess::removeSynapse) Synapse to be removed not found on list");
		}
	}

	// Clear the presynaptic link in the synapse removed
	syn->clearPresynaptic();
}

// Signal a new action potential event
void AxonProcess::signalSpikeEvent(SimTime t, SimTime isi, Number firingRate)
{
	// Create a template event object. This template is passed to
	// each synapse where a local copy made so that the synapse
	// can make local adjustments such as spike time received
	// and quantity as determined by presynaptic plasticity.
	// This strategy results in fewer dynamic memory allocations
	// at the expense of more bytes copied and more queue storage.
	ActionPotentialEvent apev(t,this);
	apev.isi(isi);
	apev.firingRate(firingRate);

	// Traverse the list, giving the template event to each synapse.
	Synapse* syn=_synapses;
	while (syn!=NULL) {
		apev.synapse(syn);
		syn->signalActionPotential(&apev);
		syn = syn->nextPresynaptic();
	}
}



// ====================================================================
// ElectricalCoupling class body
// ====================================================================



// Create a default (empty) coupling object
ElectricalCoupling::ElectricalCoupling()
{
	_gIsPreset = false;			// no value of g set yet
	_g=0;						// initialize g for consistency
}

// Create a coupling between two compartments and compute
// the coupling conductance from their internal resistances.
ElectricalCoupling::ElectricalCoupling(Compartment* c1, Compartment* c2)
{

	_nodes.push_back(c1);
	_nodes.push_back(c2);

	_gIsPreset = false;
	updateParams();

	c1->addCoupling(this);
	c2->addCoupling(this);
}

// Create a coupling object between two compartments with a specified
// conductance. Note that gVal must be scaled for compartment sizes so
// that it represents an absolute conductance (e.g. mS units) not a
// specific conductance (e.g. mS/cm_2). This is required since the specific
// conductance from compartments A to B is not the same as from
// B to A if they are different sizes.
ElectricalCoupling::ElectricalCoupling(Compartment* c1, Compartment* c2, Number gVal)
{

	_nodes.push_back(c1);
	_nodes.push_back(c2);

	_g = gVal;
	_gIsPreset = true;

	c1->addCoupling(this);
	c2->addCoupling(this);
}

// Destroy a coupling after unhooking it from associated compartments
ElectricalCoupling::~ElectricalCoupling()
{
	CompartmentVectorIt		it;

	// unhook from compartments
	for (it=_nodes.begin();it!=_nodes.end();it++) {
		if (*it!=NULL) { // should always be true, but be careful
			(*it)->removeCoupling(this);
		}
	}
}

// Get the first compartment of the coupling.
// If none, return NULL
Compartment* ElectricalCoupling::comp1()
{
	return _nodes.size()<1 ? NULL : _nodes[0];
}

// Return the number of nodes in the coupling
int ElectricalCoupling::nodeCount()
{
	return _nodes.size();
}

// Set comp1 and try to calculate g
void ElectricalCoupling::comp1(Compartment* c1) 
{
	if (_nodes.size()<1) {
		_nodes.resize(1);
	}
	_nodes[0]=c1; 
	updateParams();

	c1->addCoupling(this);

}

// Get the second compartment of the coupling.
// If none, return NULL
Compartment* ElectricalCoupling::comp2()
{
	return _nodes.size()<2 ? NULL : _nodes[1];
}

// Set comp2 and try to calculate g
void ElectricalCoupling::comp2(Compartment* c2) 
{
	if (_nodes.size()<2) {
		_nodes.resize(2);
	}
	_nodes[1]=c2; 
	updateParams();

	c2->addCoupling(this);
}

// Get the peer (i.e alternate) node to the one provided.
// If none, return NULL
Compartment* ElectricalCoupling::peer(Compartment* comp)
{
	if (_nodes.size()<2) {
		return NULL;
	}

	return _nodes[1]==comp ? _nodes[0] : _nodes[1];
}

// Set g absolutely and suppress recalculation
void ElectricalCoupling::g(double gval)
{
	_g=gval;
	_gIsPreset=true;
}

// (re)compute g if it was not specified in a constructor
void ElectricalCoupling::updateParams()
{
	if (_gIsPreset) return;

	// Reset g to a default value pending sufficient data
	// for calculation from compartment parameters.

	_g = 0;

	if (comp1() == NULL) return;
	if (comp2() == NULL) return;

	// Compute a conductance by adding resistance values
	// based on assumption that each component contributes
	// half its internal resistance to the couple.
	
	Number r1 = comp1()->Ri();
	Number r2 = comp2()->Ri();
	_g=2/(r1+r2);
}

// Get the conductance between two peers in the coupling.
// Except for capacitance, this is the Jacobian entry for
// the row corresponding to c1 and column corresponding to c2.
double ElectricalCoupling::gForPair(Compartment* c1, Compartment* c2)
{
	// Since the conductance is symmetric, return it
	return _g;
}

// Compute the current flowing with respect to one end
// of the connection. Since this class only models
// a two-way junction, the computation is straightforward
// and symmetric. Some levels of indirection are removed
// for efficiency.
double ElectricalCoupling::Iec(Compartment* comp1)
{
	if (_nodes[0]==comp1) {
		return _g*(_nodes[0]->Vm() - _nodes[1]->Vm());
	}
	else {
		return _g*(_nodes[1]->Vm() - _nodes[0]->Vm());
	}
}



// ====================================================================
// ElectricalJunction class body
// ====================================================================


// Constructors and destructor
ElectricalJunction::ElectricalJunction() {}

ElectricalJunction::ElectricalJunction(Compartment* comp)
{
	add(comp);
}

ElectricalJunction::~ElectricalJunction() {}

// Add a new compartment to the junction.
void ElectricalJunction::add(Compartment* comp)
{
	// Put the node in the list
	_nodes.push_back(comp);

	// Get new total conductance
	updateParams();

	// Hook up with the compartment added
	comp->addCoupling(this);
}

// Remove a compartment from the junction.
void ElectricalJunction::remove(Compartment* comp)
{
	CompartmentVectorIt		last;

	// Remove the compartment
	last=std::remove(_nodes.begin(),_nodes.end(),comp);
	_nodes.resize( last - _nodes.begin() );

	// Update total conductance
	updateParams();

	// Remove junction from the compartment
	comp->removeCoupling(this);
}

// Compute the total conductance for the junction.
void ElectricalJunction::updateParams()
{
	CompartmentVectorIt		it;
	
	_g = 0;
	for (it=_nodes.begin();it!=_nodes.end();it++) {
		_g += 2 / (*it)->Ri();
	}
}

// Get the conductance between two peers in the coupling.
// Except for capacitance, this is the Jacobian entry for
// the row corresponding to c1 and column corresponding to c2.
// See Compartment::buildJacobian for use of this function.
double ElectricalJunction::gForPair(Compartment* c1, Compartment* c2)
{
	// See function Iec for the equations from which this is derived.

	double		gself = 2/(c1->Ri());
	double		gratio = gself/_g;
	double		gnode;

	if (c1==c2) {
		return gself*(1-gratio);
	}
	else {
		gnode = 2/(c2->Ri());
		return gratio*gnode;
	}
}

// Compute outward current flow for a compartment
double ElectricalJunction::Iec(Compartment* comp)
{
	CompartmentVectorIt		it;

	// Because offsetting currents are computed, some extra
	// numberical precision is needed to make this accurate.
	// Hence the use of double.

	double					gself = 2/comp->Ri();
	double					gratio = gself/_g;
	double					Isum;
	double					gnode,vmnode;

	// Compute current flow by summing up individual
	// currents derived from superposition of the
	// effects of each individual compartment.

	// Note that a good bit of this could be precomputed
	// at start of simulation if this implementation turns 
	// out to be a performance problem.

	// Start with the compartment itself
	Isum = gself * (1-gratio) * comp->Vm();

	// Now look at all the other nodes.
	for(it=_nodes.begin();it!=_nodes.end();it++) {
		if (*it!=comp) {
			gnode = 2/(*it)->Ri();
			vmnode = (*it)->Vm();
			Isum -= gratio * gnode * vmnode;
		}
	}

	return Isum;
}

// If exactly two nodes in the junction, return the peer
// Otherwise, return NULL
Compartment* ElectricalJunction::peer(Compartment* comp)
{
	if (nodeCount()==2) {
		return ElectricalCoupling::peer(comp);
	}
	else {
		return NULL;
	}
}

// Inherited functions that do not apply to junctions
void ElectricalJunction::g(double gval)
{
	FatalError("(ElectricalJunction::g) Function does not apply to this class.");
}




// ====================================================================
// IonChannel class body
// ====================================================================



// Global values
Number IonChannel::_DefaultTempC = 37;					// Default body temperature
Number IonChannel::_DefaultCaXout = 2*UOM::mM;			// Ca++ external concentration
Number IonChannel::_DefaultPeakCaXin = 1*UOM::microM;	// Ca++ internal concentration

GHKTableEntry* IonChannel::_CaGHKTable = NULL;			// Ca++ GHK look-up table

Number IonChannel::_FoverRT = IonChannel::FoverRT(_DefaultTempC);

// Constructor
IonChannel::IonChannel(Number gSpecific)
{
	// Clear any pointers
	_container=NULL;
	_calciumPool=NULL;

	// Clear any cached values and set defaults
	_cachedQ10Factor = -1;			// No value set yet
	_g=0;							// Cannot get g yet (no container)
	_gModulator = 1;				// Default is no modulation
	_gParam=gSpecific;				// Use parameter value provided
	_gIsSpecific = true;			// Indicate specific vs absolute

	// Set the component id to null and allow
	// lazy initialization to catch up later
	// after constructor has run its course.
	_componentId = NullTokenId;
}

// Destructor
IonChannel::~IonChannel() 
{
	// Remove this channel from the compartment
	if (container()!=NULL) {
		container()->removeIonChannel(this);
	}
}

// Get the component token id
TokenId IonChannel::componentId()
{
	// Use lazy initialization to allow
	// component name to be retrieved via
	// subclass after base class construction.
	if (_componentId == NullTokenId) {
		_componentId = token( componentName() );
	}

	return _componentId;
}

// Save containing compartment
void IonChannel::container(Compartment* comp)
{
	// Save new containing compartment
	_container = comp;

	// Compute conductance values based on compartment
	if (container()!=NULL) {
		updateParams();
	}
}

// Set the pool that supplies Ca++ concentrations
// for channels that are modulated by Ca++.
void IonChannel::calciumPool(CalciumPool* pool)
{
	_calciumPool = pool;
}

// Add this channel to a compartment
void IonChannel::addTo(Compartment* comp)
{
	comp->addIonChannel(this);
}

// Remove this channel from its compartment
void IonChannel::removeFromCompartment()
{
	if (container()!=NULL) {
		container()->removeIonChannel(this);
	}
}

// Get membrane area to use in computing distributions
Number IonChannel::membraneArea()
{
	// Default assumption is that channel is distributed
	// uniformly over the whole membrane area of the compartment.
	if (container()!=NULL) {
		return container()->membraneArea();
	} 
	else {
		FatalError("(IonChannel::membraneArea) This channel has no container");
		return 0;
	}
}

// Get the specific conductance for this channel
Number IonChannel::gSpecific()
{
	// See which way g is specified and convert accordingly
	return gIsSpecific() ? _gParam : _gParam / membraneArea();
}

// Get the absolute conductance for this channel
Number IonChannel::gAbsolute()
{
	// See which way g is specified and convert accordingly
	return gIsAbsolute() ? _gParam : _gParam * membraneArea();
}

// Set the specific conductance for this channel
void IonChannel::gSpecific(Number gsp) 
{ 
	_gParam=gsp;
	_gIsSpecific = true;
	updateParams();
}

// Set the absolute conductance for this channel
void IonChannel::gAbsolute(Number gabs) 
{ 
	_gParam=gabs;
	_gIsSpecific = false;
	updateParams();
}

// Apply neuromodulation parameters. Subclass must
// handle any ids other than id=token("gModulator")
void IonChannel::setModParams(
	TokenId		id,
	int			numValues,
	Number*		values)
{
	static TokenId gModulatorId = token("gModulator");

	if (id==gModulatorId) {
		if (numValues<1) {
			FatalError("(IonChannel::setModParams) Too few values provided");
		}
		gModulator(values[0]);
	}
}

// Convenience function for invoking setModParams
void IonChannel::setModParam(TokenId id, Number value)
{
	setModParams(id,1,&value);
}


// Convenience function for invoking setModParams
void IonChannel::setModParam(char* modTypeName, Number value)
{
	setModParam(token(modTypeName),value);
}


// Locate the identified component within this ion channel.
// If the component is not found, return NULL.
IonChannel* IonChannel::findIonChannel(TokenId compId)
{
	return compId == componentId() ? this : NULL;
}

// Get the Q10Factor from a cached value
// Use lazy initialization to allow Q10 to be supplied
// by the subclass after base class construction.
Number IonChannel::Q10Factor()
{
	if (_cachedQ10Factor<0) {
		updateParams();
	}
	return _cachedQ10Factor;
}

// Do calculations based on current parameter values.
void IonChannel::updateParams()
{
	// Set the conductance value for use in computing currents.
	// If there is not yet an associated compartment, _g is set 
	// to 0 and then reset when the hook-up occurs later. This 
	// allows setting values during construction, before the 
	// compartment hook-up occurs.
	_g = container()==NULL ? 0 : gAbsolute()*gModulator();

	// Set the Q10 factor based on current temp and Q10
	_cachedQ10Factor = pow(Q10(),(currentTempC()-ratedTempC())/10);
}

// Get the current attributable to Ca++ ions.
// This can be overridden if the current is a mixture of ions.
Number IonChannel::ICa() 
{ 
	return Iion(); 
}

// Set the default temperature and make necessary adjustments
void IonChannel::defaultTempC(Number temp) {

	// Set the global value
	_DefaultTempC = temp;

	// Update dependents
	_FoverRT = FoverRT(temp);
}

// Compute F/RT using the default preparation temperature
Number IonChannel::FoverRT()
{
	return _FoverRT;
}

// Compute F/(RT) for a given temperature (C)
Number IonChannel::FoverRT(Number tempC)
{
	using namespace UOM;

	return Faraday/(GasConstant*degKfromC(tempC));
}

// Set default internal Ca++ concentration
void IonChannel::defaultPeakCaXin(Number x)
{
	_DefaultPeakCaXin = x;
	if (pCaGHKTable()!=NULL) {
		loadCaGHKTable();
	}
}

// Set default external Ca++ concentration
void IonChannel::defaultCaXout(Number x)
{
	_DefaultCaXout = x;
	if (pCaGHKTable()!=NULL) {
		loadCaGHKTable();
	}
}

// Compute Ca++ effective potential using GHK formula
Number IonChannel::ghkCaEffectivePotential(
	Number			v,
	Number			CaXin,
	Number			CaXout,
	Number			tempC)
{
	const Number	epsilon = VStepForIndex/128;	// v=0 test range
	const Number	z=2;							// Ca++ charge

	// Check for singularity at v=0
	if (fabs(v)>epsilon) {

		// Compute with full formula
		double ezvFRT = exp(z*v*FoverRT(tempC));
		return v*(1-CaXin/CaXout*ezvFRT)/(1-ezvFRT);
	}
	else {
		// Return limit at V=0
		return (1-CaXin/CaXout)/(-z*FoverRT(tempC));
	}
}

// Compute an effective conductance by differentiating
// the effective potential with respect to V. When
// multiplied by a constant conductance this gives the
// effective conductance for a calcium channel.
// The derivative is extracted numerically.
Number IonChannel::ghkCaEffectiveCond(
	Number			v,
	Number			CaXin,
	Number			CaXout,
	Number			tempC)
{
	const Number	dVm = VStepForIndex/32;
	Number			E1,E2;

	E2 = ghkCaEffectivePotential(v+dVm,CaXin,CaXout,tempC);
	E1 = ghkCaEffectivePotential(v-dVm,CaXin,CaXout,tempC);
	return (E2-E1)/(2*dVm);
}

// Load the GHK table of values by voltage
// Since this might be invoked for subclasses,
// be sure to use appropriate accessors.
void IonChannel::loadCaGHKTable()
{
	GHKTableEntry*	ghkTbl;
	int				k;
	Number			v;
	Number			CaXin = defaultPeakCaXin();
	Number			CaXout = defaultCaXout();
	Number			tempC = currentTempC();

	// Delete any previously allocated table
	if (*pCaGHKTable()!=NULL) {
		delete[] *pCaGHKTable();
	}

	// Allocate a new table and remember where it is
	*pCaGHKTable() = ghkTbl = new GHKTableEntry[VTableSize];

	// Go through the voltages one step at a time
	for (k=0,v=VMinForIndex; k<VTableSize; k++, v+=VStepForIndex ) {
		ghkTbl[k].Veff=ghkCaEffectivePotential(v,CaXin,CaXout,tempC);
		ghkTbl[k].Geff=ghkCaEffectiveCond(v,CaXin,CaXout,tempC);
	}
}

// Look up GHK Ca++ effective potential using defaults
Number IonChannel::ghkCaEffectivePotential()
{
	GHKTableEntry* ent=(*pCaGHKTable())+container()->VmIndex();
	Number s = container()->VmRem()/VStepForIndex;

	return VTableInterp(s,
		(ent-1)->Veff, ent->Veff, (ent+1)->Veff, (ent+2)->Veff);
}

// Look up GHK Ca++ effective conductance using defaults
Number IonChannel::ghkCaEffectiveCond()
{
	GHKTableEntry* ent=(*pCaGHKTable())+container()->VmIndex();
	Number s = container()->VmRem()/VStepForIndex;
	
	return VTableInterp(s,
		(ent-1)->Geff, ent->Geff, (ent+1)->Geff, (ent+2)->Geff);
}

// Update all params as of start of the simulation and load the
// Ca++ GHK table on start of the simulation (if not done already)
void IonChannel::simulationStarted()
{
	updateParams();

	if (*pCaGHKTable()==NULL) {
		loadCaGHKTable();
	}
}

// Delete the Ca++ GHK table when the simulation ends
void IonChannel::simulationEnded()
{
	delete[] *pCaGHKTable();
	*pCaGHKTable() = NULL;
}



// ====================================================================
// DummyIonChannel class body
// ====================================================================



// Constructor and destructor
DummyIonChannel::DummyIonChannel(int numStateVar,CalciumPool* pool)
{
	_numStateVar = numStateVar;
	if (pool!=NULL) {
		calciumPool(pool);
	}
}

DummyIonChannel::~DummyIonChannel() {}



// ====================================================================
// MultiGateIonChannel class body
// ====================================================================



// Constructor and destructor
MultiGateIonChannel::MultiGateIonChannel(Number gSpVal)
: IonChannel(gSpVal) {}

MultiGateIonChannel::~MultiGateIonChannel() 
{
	IonChannelVectorIt		it;

	// Delete all held gates
	for (it=_gates.begin();it!=_gates.end(); it++) {
		delete *it;
	}
}

// make a gate variable known
void MultiGateIonChannel::addGate(IonChannel* gate)
{
	_gates.push_back(gate);

	// Let gate know about its container
	gate->container( container() );

	// Hook to model
	gate->model( model() );
}

// when new container provided, pass along to gates also
void MultiGateIonChannel::container(Compartment* comp)
{
	IonChannelVectorIt		it;

	IonChannel::container(comp);

	for (it=_gates.begin(); it!=_gates.end(); it++) {
		(*it)->container(comp);
	}
}

// When new model is provided, pass along to gates.
void MultiGateIonChannel::model(Model* m)
{
	IonChannelVectorIt		it;

	// Add this object to the model
	IonChannel::model(m);

	// Add any gates to the model
	for (it=_gates.begin(); it!=_gates.end(); it++) {
		(*it)->model(m);
	}

}

// Check that the required number of gates have been added by the
// time the simulation is started.
void MultiGateIonChannel::simulationStarted()
{
	string	msg1("(MultiGateIonChannel::simulationStarted) "
			"Incorrect number of gates for ");

	if (requiredNumGates() != _gates.size() ) {
		FatalError(msg1+componentName());
	}

	// Let superclass do whatever is needed for simulation started
	IonChannel::simulationStarted();
}

// Supply a default component name taking the name of the first
// gate if there is one or just a default if not.
const char* MultiGateIonChannel::componentName()
{
	static char*	defaultName = "MultiGateIonChannel";

	return _gates.size()>0 ? _gates[0]->componentName() : defaultName;
}

// Add to components to probe when reporting.
// This includes the current object and any held gates.
void MultiGateIonChannel::addToComponentsToProbe(ModelComponentVector& comps)
{

	int i;

	// Add this object
	comps.push_back(this);

	// Add all gates held
	for (i=0;i<_gates.size();i++) {
		_gates[i]->addToComponentsToProbe(comps);
	}
}

// Locate the identified component within this ion channel.
// If the component is not found, return NULL.
IonChannel* MultiGateIonChannel::findIonChannel(TokenId compId)
{
	IonChannelVectorIt it;

	// Check this object for a match
	if (componentId()==compId)
		return this;

	// Check gates and return the first match found
	for (it=_gates.begin(); it!=_gates.end(); it++) {
		if ( (*it)->componentId()==compId) 
			return *it;
	}

	// Nothing found here
	return NULL;
}



// ====================================================================
// M1HIonChannel class body
// ====================================================================



// Constructor and destructor
M1HIonChannel::M1HIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M1HIonChannel::~M1HIonChannel() {}

// Get conductance
Number M1HIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*h;
}

// Get both conductance and current via Ohm's law
void M1HIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	Gout = g()*m*h;
	Iout = Gout*(Vm()-Vrev());
}


// ====================================================================
// M2HIonChannel class body
// ====================================================================



// Constructor and destructor
M2HIonChannel::M2HIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M2HIonChannel::~M2HIonChannel() {}

// Compute conductance
Number M2HIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*h;
}

// Get both conductance and current via Ohm's law
void M2HIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	Gout = g()*m*m*h;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// M3HIonChannel class body
// ====================================================================



// Constructor and destructor
M3HIonChannel::M3HIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M3HIonChannel::~M3HIonChannel() {}

// Compute conductance
Number M3HIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*m*h;
}

// Get both conductance and current via Ohm's law
void M3HIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	Gout = g()*m*m*m*h;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// M4HIonChannel class body
// ====================================================================



// Constructor and destructor
M4HIonChannel::M4HIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M4HIonChannel::~M4HIonChannel() {}

// Compute conductance
Number M4HIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*m*m*h;
}

// Get both conductance and current via Ohm's law
void M4HIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	Gout = g()*m*m*m*m*h;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// M1HCaIonChannel class body
// ====================================================================



// Constructor and destructor
M1HCaIonChannel::M1HCaIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M1HCaIonChannel::~M1HCaIonChannel() {}

// Compute conductance
Number M1HCaIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*h*ghkCaEffectiveCond();
}

// Compute current
Number M1HCaIonChannel::Iion()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*h*ghkCaEffectivePotential();
}




// ====================================================================
// M2HCaIonChannel class body
// ====================================================================



// Constructor and destructor
M2HCaIonChannel::M2HCaIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M2HCaIonChannel::~M2HCaIonChannel() {}

// Compute conductance
Number M2HCaIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*h*ghkCaEffectiveCond();
}

// Compute current
Number M2HCaIonChannel::Iion()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*h*ghkCaEffectivePotential();
}



// ====================================================================
// M3HCaIonChannel class body
// ====================================================================



// Constructor and destructor
M3HCaIonChannel::M3HCaIonChannel(Number gSpVal) : MnHIonChannel(gSpVal) {}

M3HCaIonChannel::~M3HCaIonChannel() {}

// Compute conductance
Number M3HCaIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*m*h*ghkCaEffectiveCond();
}

// Compute current
Number M3HCaIonChannel::Iion()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();

	return g()*m*m*m*h*ghkCaEffectivePotential();
}



// ====================================================================
// M1HSIonChannel class body
// ====================================================================



// Constructor and destructor
M1HSIonChannel::M1HSIonChannel(Number gSpVal) : MnHSIonChannel(gSpVal) {}

M1HSIonChannel::~M1HSIonChannel() {}

// Compute conductance
Number M1HSIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	return g()*m*h*s;
}

// Get both conductance and current via Ohm's law
void M1HSIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	Gout = g()*m*h*s;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// M2HSIonChannel class body
// ====================================================================



// Constructor and destructor
M2HSIonChannel::M2HSIonChannel(Number gSpVal) : MnHSIonChannel(gSpVal) {}

M2HSIonChannel::~M2HSIonChannel() {}

// Compute conductance
Number M2HSIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	return g()*m*m*h*s;
}

// Get both conductance and current via Ohm's law
void M2HSIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	Gout = g()*m*m*h*s;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// M3HSIonChannel class body
// ====================================================================



// Constructor and destructor
M3HSIonChannel::M3HSIonChannel(Number gSpVal) : MnHSIonChannel(gSpVal) {}

M3HSIonChannel::~M3HSIonChannel() {}

// Compute conductance
Number M3HSIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	return g()*m*m*m*h*s;
}

// Get both conductance and current via Ohm's law
void M3HSIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	Gout = g()*m*m*m*h*s;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// M4HSIonChannel class body
// ====================================================================



// Constructor and destructor
M4HSIonChannel::M4HSIonChannel(Number gSpVal) : MnHSIonChannel(gSpVal) {}

M4HSIonChannel::~M4HSIonChannel() {}

// Compute conductance
Number M4HSIonChannel::conductance()
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	return g()*m*m*m*m*h*s;
}

// Get both conductance and current via Ohm's law
void M4HSIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number m=_gates[0]->value();
	Number h=_gates[1]->value();
	Number s=_gates[2]->value();

	Gout = g()*m*m*m*m*h*s;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// HHIonChannel class body
// ====================================================================



// Constructor and destructor
HHIonChannel::HHIonChannel(Number gSpVal) : IonChannel(gSpVal) {}

HHIonChannel::~HHIonChannel() {}

// Set initial state assuming one state variable
// and t=infinity limit computed from current
// alpha and beta values.
void HHIonChannel::setInitialState()
{
	Number		a = alpha();
	Number		b = beta();

	stateValue(0) = a/(a+b);
}


// Calculate the derivatives of a state vector
// from the membrane voltage and current state.
void HHIonChannel::computeDerivatives()
{
	Number a=alpha();
	Number b=beta();
	Number x=stateValue(0);

	derivValue(0) = a-(a+b)*x;
}

// Update the state vector for a time step of size h.
void HHIonChannel::localStateUpdate(SimTime h, CNStepType stepType)
{
	// This logic closely follows that of compute derivatives.
	// If the local ODE is dy/dt=f(y), the local update is:
	// y(t+h)=y(t)+h/(1-h/2*df/dy)*f(y)

	Number a=alpha();
	Number b=beta();
	Number x=stateValue(0);

	// Apply the semi-implicit trapezoid rule
	if (stepType==CNStartingHalfStep) {
		stateValue(0) += h/(1+h*(a+b))*(a-(a+b)*x);
	}
	else {
		stateValue(0) += h/(1+h/2*(a+b))*(a-(a+b)*x);
	}
}

// Print alpha-beta values for debugging.
// If no path name give, print to stdout.
// This sets up a temporary model and compartment
// only for the duration of the print.
void HHIonChannel::printAlphaBeta(char* pathName)
{

	FILE*				out;
	Number				v;
	float				a,b;
	int					k;

	Model*				oldModel;
	Compartment*		oldContainer;
	int					oldSVOffset;

	ClockSolver			newSolver;		// temporary solver - referenced but not used
	Model				newModel;		// temporary model
	Compartment			newContainer;	// temporary compartment

	// Set up a separate environment for holding state and switch to using
	// a new container object where the voltage can be set.

	oldContainer = container();
	oldModel = model();
	oldSVOffset = svOffset();

	// Hook up with new model and container without unhooking from the old ones
	_model = NULL;
	_container = NULL;
	newModel.solver(&newSolver);
	newContainer.model(&newModel);
	newContainer.addIonChannel(this);	

	// Make sure states are set at the beginning
	newModel.simulationStarted();
	newContainer.simulationStarted();
	simulationStarted();

	newContainer.setInitialState();
	setInitialState();

	// Open the output file
	if (pathName==NULL) {
		out = stdout;
	}
	else {
		out = fopen(pathName,"w");
		if (out==NULL) {
			FatalError("HHIonChannel::printAlphaBetaTable) File open failed");
		}
	}

	// Print the table to file
	for (k=0,v=VMinForIndex; k<VTableSize; v+=VStepForIndex,k++) {
		newContainer.Vm(v);
		a=alpha();
		b=beta();
		fprintf(out,"%g,%g,%g\n",v,a,b);
	}

	// Close the file and wrap up any simulation activities
	if (out!=stdout) {
		fclose(out);
	}
	simulationEnded();

	// Restore the old containing compartment and model
	newModel.solver(NULL);
	newContainer.removeIonChannel(this);
	_container = oldContainer;
	_model = oldModel;
	if (oldModel!=NULL) {
		svOffset(oldSVOffset);	// restore state vector array pointer
	}
}

// Function for computing rates of the form:
// r*v/(1-exp(-v/s)) where v can equal zero.
// Use of double avoids warnings when using constants as parameters. 
double HHIonChannel::linoidRate(double r, double v, double s)
{
	const double	epsilon = 1e-6;
	double			ratio = v/s;

	if (fabs(ratio)>epsilon) {
		return r*v/(1-exp(-ratio));
	}
	else {
 		return r*s;
	}
}



// ====================================================================
// BlendedIonChannel class body
// ====================================================================



// Constructors and destructor
BlendedIonChannel::BlendedIonChannel(
	Number gSpVal, 
	HHIonChannel* g1, 
	HHIonChannel* g2, 
	Number r) : HHIonChannel(gSpVal)
{
	_gate1 = g1;
	_gate2 = g2;
	_blendRatio = r;
}

BlendedIonChannel::~BlendedIonChannel()
{
	delete _gate1;
	delete _gate2;
}

// Set gate value and connect with compartment
void BlendedIonChannel::gate1(HHIonChannel* gate)
{
	// Delete any old gate assigment
	delete _gate1;

	// Assign new gate
	_gate1 = gate;
	_gate1->container( container() );
}

// Set gate value and connect with compartment
void BlendedIonChannel::gate2(HHIonChannel* gate)
{
	// Delete any old gate assigment
	delete _gate2;

	// Assign new gate
	_gate2 = gate;
	_gate2->container( container() );
}

// Set blend ratio
void BlendedIonChannel::blendRatio(Number ratio)
{
	// Make sure ratio is in bounds
	if (ratio<0 || ratio>1) {
		FatalError("(BlendedIonChannel::blendRatio) ratio is not in interval 0 through 1.");
	}

	_blendRatio = ratio;
}

// Set container and pass along to constituents
void BlendedIonChannel::container(Compartment* comp)
{
	// Set superclass do hookup for this object
	HHIonChannel::container(comp);

	// Pass container along to constituents
	if (_gate1!=NULL) {
		_gate1->container(comp);
	}
	if (_gate2!=NULL) {
		_gate2->container(comp);
	}
}

// Pass simulation start event along to constituents
void BlendedIonChannel::simulationStarted()
{
	// Let superclass do whatever is needed for this object
	HHIonChannel::simulationStarted();

	// Notify constituents
	_gate1->simulationStarted();
	_gate2->simulationStarted();
}

// Pass simulation end event along to constituents
void BlendedIonChannel::simulationEnded()
{
	// Let superclass do whatever is needed for this object
	HHIonChannel::simulationEnded();

	// Notify constituents
	_gate1->simulationEnded();
	_gate2->simulationEnded();
}

// Return a blended alpha value
Number BlendedIonChannel::alpha()
{
	Number r = blendRatio();
	Number a1 = gate1()->alpha();
	Number a2 = gate2()->alpha();

	return (1-r)*a1+r*a2;
}

// Return a blended beta value
Number BlendedIonChannel::beta()
{
	Number r = blendRatio();
	Number b1 = gate1()->beta();
	Number b2 = gate2()->beta();

	return (1-r)*b1+r*b2;
}




// ====================================================================
// Order1BlendedIonChannel class body
// ====================================================================



// Constructors and destructor
Order1BlendedIonChannel::Order1BlendedIonChannel(
	Number gSpVal, 
	HHIonChannel* g1, 
	HHIonChannel* g2, 
	Number r) :	BlendedIonChannel(gSpVal,g1,g2,r) {}

Order1BlendedIonChannel::~Order1BlendedIonChannel() {}

// Get conductance
Number Order1BlendedIonChannel::conductance()
{
	return g()*value();
}

// Get current via Ohm's law
Number Order1BlendedIonChannel::Iion()
{
	return g()*value()*(Vm()-Vrev());
}

// Get conductance and current together
void Order1BlendedIonChannel::condAndIion(Number& Gout, Number& Iout)
{
	Gout = g()*value();
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// VoltageDepTabChannel class body
// ====================================================================



// Static Variables
bool VoltageDepTabChannel::_alphaBetaLoadInProg = false;
VoltageDepTabChannel* VoltageDepTabChannel::_alphaBetaLoadObject = NULL;

// Constructors and Destructor
VoltageDepTabChannel::VoltageDepTabChannel(Number g)
: HHIonChannel(g) 
{
	// Set cache values to known state
	_pABTable = NULL;
}

VoltageDepTabChannel::~VoltageDepTabChannel() {}

// Return true if the alpha beta tabel has already
// been loaded for relevant voltages.
bool VoltageDepTabChannel::alphaBetaTableLoaded() 
{
	return *pAlphaBetaTable()!=NULL;
}

// Return true if loading of an alpha beta table
// is in progress. This obviously assumes only
// load can be ongoing at a time.
bool VoltageDepTabChannel::alphaBetaLoadInProg()
{
	return _alphaBetaLoadInProg;
}

// Load the alpha beta table one entry at a time.
// Only do this once per class.
void VoltageDepTabChannel::loadAlphaBetaTable()
{
	// Only load the table once per class.
	if (alphaBetaTableLoaded() ) 
		return;

	AlphaBetaEntry*	abTab;
	int				k;
	Number			v;
	Number			alpha,beta,xinf,tau;
	Number			xmin = xinfAbsMin();

	// Provide way for subclasses to know load is in progress
	_alphaBetaLoadInProg = true;
	_alphaBetaLoadObject = this;

	// Allocate the table and remember where it is
	*(pAlphaBetaTable()) = abTab = new AlphaBetaEntry[VTableSize];

	// Go through the voltages one step at a time and load
	// values for alpha, beta, xinf, and tau.
	for (k=0,v=VMinForIndex; k<VTableSize; k++, v+=VStepForIndex ) {

		// See whether alpha/beta or xinf/tau is being used
		alpha = alphaForTable(v);
		if (alpha>=0) {
			beta = betaForTable(v);
			xinf=alpha/(alpha+beta);
			tau=1/(alpha+beta);
		}
		else {

			// Get xinf and tau values and check for validity
			xinf = xinfForTable(v);
			if (xinf<0) {
				FatalError("(VoltageDepTabChannel::loadAlphaBetaTable) "
					"Subclass missing both alpha and xinf functions");
			}
			tau = tauForTable(v);
			if (tau<0) {
				FatalError("(VoltageDepTabChannel::loadAlphaBetaTable) "
					"Invalid value returned for tau");
			}

			// Apply an exponent adjustment to xinf if provided.
			if (xinfExponent()!=1) {
				xinf = pow(xinf,xinfExponent());
			}

			// Apply any multiplication by a sigmoidal function to the xinf to be loaded
			if (xinfMultK()!=0) {
				xinf /= 1+exp(-(v-xinfMultVhalf())/xinfMultK());
			}

			// Also impose a minimum value on xinf (useful for leaky gates)
			xinf = xmin+(1-xmin)*xinf;


			// Quietly impose a lower limit on tau to keep impossibly low
			// value from adversely affecting the simulation time step.
			if (tau<tauAbsMin() ) {
				tau = tauAbsMin();
			}

			// Get alpha and beta equivalent to xinf and tau
			alpha=xinf/tau;
			beta=(1.0-xinf)/tau;
		}

		// Update the table entry 
		abTab[k].alphaValue=alpha;
		abTab[k].betaValue=beta;
		abTab[k].xinfValue=xinf;
		abTab[k].tauValue=tau;
	}

	// All done, tell interested subclasses
	_alphaBetaLoadInProg = false;
	_alphaBetaLoadObject = NULL;
}

// Delete the alpha beta table.
void VoltageDepTabChannel::deleteAlphaBetaTable()
{
	if (alphaBetaTableLoaded()) {
		delete[] *pAlphaBetaTable();
		*pAlphaBetaTable() = NULL;
	}
}

// Initialize when the simulation starts
void VoltageDepTabChannel::simulationStarted()
{
	// Let superclass do its initializations
	HHIonChannel::simulationStarted();

	// Load the alpha beta table with current values
	loadAlphaBetaTable();

	// Save cached address of table to speed up processing later
	_pABTable = pAlphaBetaTable();
}

// Clean up when the simulation ends
void VoltageDepTabChannel::simulationEnded()
{
	// Let superclass do its clean-up also
	HHIonChannel::simulationEnded();

	// Delete the alpha beta table
	if (alphaBetaTableLoaded() ) {
		deleteAlphaBetaTable();
	}
}

// Compute alpha function from table using linear interpolation
Number VoltageDepTabChannel::alpha() 
{ 
	AlphaBetaEntry* abEnt = *_pABTable+container()->VmIndex();
	Number s = container()->VmRem()/VStepForIndex;
	return VTableInterp(s, 
		(abEnt-1)->alphaValue, abEnt->alphaValue, 
		(abEnt+1)->alphaValue, (abEnt+2)->alphaValue);
}

// Compute beta function from table using linear interpolation
Number VoltageDepTabChannel::beta() 
{ 
	AlphaBetaEntry* abEnt = *_pABTable+container()->VmIndex();
	Number s = container()->VmRem()/VStepForIndex;
	return VTableInterp(s, 
		(abEnt-1)->betaValue, abEnt->betaValue, 
		(abEnt+1)->betaValue, (abEnt+2)->betaValue);
}

// Compute xinf function from table using linear interpolation
Number VoltageDepTabChannel::xinf() 
{ 
	AlphaBetaEntry* abEnt = *_pABTable+container()->VmIndex();
	Number s = container()->VmRem()/VStepForIndex;
	return VTableInterp(s, 
		(abEnt-1)->xinfValue, abEnt->xinfValue, 
		(abEnt+1)->xinfValue, (abEnt+2)->xinfValue);
}

// Compute tau function from table using linear interpolation
Number VoltageDepTabChannel::tau() 
{ 
	AlphaBetaEntry* abEnt = *_pABTable+container()->VmIndex();
	Number s = container()->VmRem()/VStepForIndex;
	return VTableInterp(s,
		(abEnt-1)->tauValue, abEnt->tauValue, 
		(abEnt+1)->tauValue, (abEnt+2)->tauValue);
}

// Fast path computation of derivatives.
// Note that xinf() and tau() are bypassed
// and if they are overridden, this function
// should be overridden as well.
void VoltageDepTabChannel::computeDerivatives()
{
	// Interpolate values for xinf and tau
	AlphaBetaEntry*		abEnt = *_pABTable+container()->VmIndex();

	Number		x = stateValue(0);
	Number		s = container()->VmRem()/VStepForIndex;
	Number		xinf,tau;

	xinf= VTableInterp(s, 
			(abEnt-1)->xinfValue, abEnt->xinfValue, 
			(abEnt+1)->xinfValue, (abEnt+2)->xinfValue);
	tau = VTableInterp(s,
			(abEnt-1)->tauValue, abEnt->tauValue, 
			(abEnt+1)->tauValue, (abEnt+2)->tauValue);

	// Apply xinf and tau to get derivative
	// loadAlphaBetaTable prevents tau=0 case from occurring
	derivValue(0) = (xinf-x)/tau;
}

// Fast path computation of local state update via implicit rule.
// Note that xinf() and tau() are bypassed and if they are overridden,
// this function should be overriden as well.
void VoltageDepTabChannel::localStateUpdate(SimTime h, CNStepType stepType)
{
	// This logic closely follows that of compute derivatives.
	// If the local ODE is dy/dt=f(y), the local update is:
	// y(t+h)=y(t)+h/(1-h/2*df/dy)*f(y)

	// Interpolate values for xinf and tau
	AlphaBetaEntry*		abEnt = *_pABTable+container()->VmIndex();

	Number		x = stateValue(0);
	Number		s = container()->VmRem()/VStepForIndex;
	Number		xinf,tau;

	xinf= VTableInterp(s, 
			(abEnt-1)->xinfValue, abEnt->xinfValue,
			(abEnt+1)->xinfValue, (abEnt+2)->xinfValue);
	tau = VTableInterp(s,
			(abEnt-1)->tauValue, abEnt->tauValue,
			(abEnt+1)->tauValue, (abEnt+2)->tauValue);
	
	// Apply the semi-implicit trapezoid rule
	if (stepType==CNStartingHalfStep) {
		stateValue(0) += h/(tau+h)*(xinf-x);
	}
	else {
		stateValue(0) += h/(tau+h/2)*(xinf-x);
	}	
}

// Print alpha-beta table for debugging
// If no path name is given, print to stdout.
// The alpha-beta table is loaded as a side-effect.
void VoltageDepTabChannel::printAlphaBetaTable(char* pathName)
{
	FILE*		out;
	Number		v;
	int			k;

	if (pathName==NULL) {
		out = stdout;
	}
	else {
		out = fopen(pathName,"w+");
		if (out==NULL) {
			FatalError("(VoltageDepTabChannel::printAlphaBetaTable) File open failed");
		}
	}

	// Make sure table is loaded
	loadAlphaBetaTable();
	AlphaBetaEntry*	abt = *pAlphaBetaTable();

	// Print the table to file
	for (k=0,v=VMinForIndex; k<VTableSize; v+=VStepForIndex,k++) {
		fprintf(out,"%g,%g,%g,%g,%g\n",
			double(v),
			double(abt[k].alphaValue),
			double(abt[k].betaValue),
			double(abt[k].xinfValue),
			double(abt[k].tauValue) );
	}

	// Close the file
	if (out!=stdout) {
		fclose(out);
	}
}



// ====================================================================
// EnergyBarrierTabChannel class body
// ====================================================================



// Static variables
Number EnergyBarrierTabChannel::_cachedFoverRT = 0;
Number EnergyBarrierTabChannel::_cachedRate = 0;
Number EnergyBarrierTabChannel::_cachedZeta = 0;

// Constructors and destructor
EnergyBarrierTabChannel::EnergyBarrierTabChannel(Number gSpVal)
: VoltageDepTabChannel(gSpVal) 
{
	// Clear cached values
	_cachedQ10FactorForTauMin = -1;		// no value set
}

EnergyBarrierTabChannel::~EnergyBarrierTabChannel() {}

// Clear cached values before loading the alpha beta table
void EnergyBarrierTabChannel::loadAlphaBetaTable()
{
	// Load the alpha beta table (actual load is one time only)
	_cachedFoverRT = 0;
	_cachedRate = -1;
	_cachedZeta = 0;

	// Let the superclass do its thing
	VoltageDepTabChannel::loadAlphaBetaTable();
}

// Compute a value of Q10 adjustment factor and cache
// it for use during table loading.
Number EnergyBarrierTabChannel::Q10FactorForTauMin()
{
	if (_cachedQ10FactorForTauMin<=0) {
		_cachedQ10FactorForTauMin = 
			pow(Q10ForTauMin(), (currentTempC()-ratedTempC())/10);
	}
	return _cachedQ10FactorForTauMin;
}

// Compute a value of zeta from the slope parameter
// This is done at the rated temperature.
Number EnergyBarrierTabChannel::zeta()
{
	Number zetaValue;

	// See if a cached value is available
	if (_alphaBetaLoadObject == this &&
		_cachedZeta != 0) {
		return _cachedZeta;
	}

	// Make sure either slope or zeta was overridden by subclass
	if (slope()==0) {
		FatalError("(EnergyBarrierTabChannel::zeta) slope parameter not provided");
	}

	// Do the computation so that zeta*F/RT = 1/slope
	zetaValue = 1/(FoverRTAtRatedTemp()*slope());

	if (_alphaBetaLoadObject == this) {
		_cachedZeta = zetaValue;
	}
	return zetaValue;
}

// Compute a default rate based on tauMax at rated temperature.
// We assume that gamma is not voltage dependent (or else this
// function should be overridden by the subclass).
Number EnergyBarrierTabChannel::rate()
{
	Number rateValue;

	// See if a cached value is available
	if (_alphaBetaLoadObject == this &&
		_cachedRate >= 0) {
		return _cachedRate;
	}

	Number		tauDiff = tauMax()-tauMin();
	Number		g = gamma();
	Number		zFRT = zeta()*FoverRTAtRatedTemp();
	Number		vmin,denom;

	// If tauDiff == 0, the time constant is fixed at one value.
	// Return the special value of 0, which would otherwise be invalid.
	if ( fabs(tauDiff) < 1e-8 * UOM::sec) { // Test for approx zero.
											// This should probably be
											// EpsilonTime instead.
											// Fix in the future.
		if (_alphaBetaLoadObject == this) {
			_cachedRate = 0;
		}
		return 0;
	}

	// We want to find the minimum value for the denominator in the
	// expression for tau, but gamma = 0 or 1 is a special case.
	if (g==0 ||g==1) {
		rateValue = 1/tauDiff;
		if (_alphaBetaLoadObject == this) {
			_cachedRate = rateValue;
		}
		return rateValue;
	}

	// Handle other special cases in which the result can be had 
	// with less effort. Here v=vhalf gives the maximum tau value.
	if (zFRT==0 || g==0.5) {
		rateValue = 1/(2*tauDiff);
		if (_alphaBetaLoadObject == this) {
			_cachedRate = rateValue;
		}
		return rateValue;
	}

	// Determine the voltage offset from Vhalf at which the 
	// denominator in the formula for tau is a minimum.
	vmin = log((1-g)/g)/zFRT;

	// Plug in the value for vmin = v-vhalf and find the denominator value
	denom = exp(zFRT*g*vmin)+exp(-zFRT*(1-g)*vmin);

	// Get the rate value that gets the right tauMax
	rateValue = 1/(tauDiff*denom);
	if (_alphaBetaLoadObject == this) {
		_cachedRate = rateValue;
	}
	return rateValue;
}

// Return the value of F/RT for the rated temperature
Number EnergyBarrierTabChannel::FoverRTAtRatedTemp()
{
	if (_alphaBetaLoadObject == this &&
		_cachedFoverRT != 0) {
		return _cachedFoverRT;
	}

	Number FoverRTValue = FoverRT( ratedTempC() );
	if (_alphaBetaLoadObject == this) {
		_cachedFoverRT = FoverRTValue;
	}
	return FoverRTValue;
}

// Compute equilibrium value of state variable
// at the current simulation temperature.
Number EnergyBarrierTabChannel::xinfForTable(Number v)
{
	Number	zFRT = zetaAtV(v)*FoverRT();

	return 1/( 1+exp( -zFRT*(v-Vhalf()) ));
}

// Compute the time constant for a voltage adjusting for
// the current simulation temperature.
Number EnergyBarrierTabChannel::tauForTable(Number v)
{
	// Rate = 0 is a special case indicating that only tauMin is used.
	Number rate = rateAtV(v);
	if (rate==0) {
		return tauMinAtV(v)/Q10FactorForTauMin();
	}

	// Otherwise, do the full computation
	Number g = gammaAtV(v);
	Number zvFRT = zetaAtV(v)*(v - Vhalf())*FoverRT();

	// Compute tau from 1/(alpha+beta)
	Number tau0=tauMinAtV(v);
	Number tau1=1/( rate*( exp(g*zvFRT)+exp(-(1-g)*zvFRT) ));

	// Return a temperature adjusted tau
	return tau1/Q10Factor()+tau0/Q10FactorForTauMin();
}



// ====================================================================
// Order1EnergyBarrierTabChannel class body
// ====================================================================



// Constructors and destructor
Order1EnergyBarrierTabChannel::Order1EnergyBarrierTabChannel(Number gSpVal)
: EnergyBarrierTabChannel(gSpVal) {}

Order1EnergyBarrierTabChannel::~Order1EnergyBarrierTabChannel() {}

// Get conductance
Number Order1EnergyBarrierTabChannel::conductance()
{
	Number x = value();

	return g()*x;
}

// Get current
Number Order1EnergyBarrierTabChannel::Iion()
{
	Number x = value();

	return g()*x*(Vm()-Vrev());
}

// Get conductance and current together
void Order1EnergyBarrierTabChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number x = value();

	Gout = g()*x;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// Order2EnergyBarrierTabChannel class body
// ====================================================================



// Constructors and destructor
Order2EnergyBarrierTabChannel::Order2EnergyBarrierTabChannel(Number gSpVal)
: EnergyBarrierTabChannel(gSpVal) {}

Order2EnergyBarrierTabChannel::~Order2EnergyBarrierTabChannel() {}

// Get conductance
Number Order2EnergyBarrierTabChannel::conductance()
{
	Number x = value();

	return g()*x*x;
}

// Get current
Number Order2EnergyBarrierTabChannel::Iion()
{
	Number x = value();

	return g()*x*x*(Vm()-Vrev());
}

// Get conductance and current together
void Order2EnergyBarrierTabChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number x = value();

	Gout = g()*x*x;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// Order3EnergyBarrierTabChannel class body
// ====================================================================



// Constructors and destructor
Order3EnergyBarrierTabChannel::Order3EnergyBarrierTabChannel(Number gSpVal)
: EnergyBarrierTabChannel(gSpVal) {}

Order3EnergyBarrierTabChannel::~Order3EnergyBarrierTabChannel() {}

// Get conductance
Number Order3EnergyBarrierTabChannel::conductance()
{
	Number x = value();

	return g()*x*x*x;
}

// Get current
Number Order3EnergyBarrierTabChannel::Iion()
{
	Number x = value();

	return g()*x*x*x*(Vm()-Vrev());
}

// Get conductance and current together
void Order3EnergyBarrierTabChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number x = value();

	Gout = g()*x*x*x;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// Order4EnergyBarrierTabChannel class body
// ====================================================================



// Constructors and destructor
Order4EnergyBarrierTabChannel::Order4EnergyBarrierTabChannel(Number gSpVal)
: EnergyBarrierTabChannel(gSpVal) {}

Order4EnergyBarrierTabChannel::~Order4EnergyBarrierTabChannel() {}

// Get conductance
Number Order4EnergyBarrierTabChannel::conductance()
{
	Number x = value();

	return g()*x*x*x*x;
}

// Get current
Number Order4EnergyBarrierTabChannel::Iion()
{
	Number x = value();

	return g()*x*x*x*x*(Vm()-Vrev());
}

// Get conductance and current together
void Order4EnergyBarrierTabChannel::condAndIion(Number& Gout, Number& Iout)
{
	Number x = value();

	Gout = g()*x*x*x*x;
	Iout = Gout*(Vm()-Vrev());
}



// ====================================================================
// Order1CaEnergyBarrierTabChannel class body
// ====================================================================



// Constructors and destructor
Order1CaEnergyBarrierTabChannel::Order1CaEnergyBarrierTabChannel(Number gSpVal)
: EnergyBarrierTabChannel(gSpVal) {}

Order1CaEnergyBarrierTabChannel::~Order1CaEnergyBarrierTabChannel() {}

// Get conductance
Number Order1CaEnergyBarrierTabChannel::conductance()
{
	Number x = value();

	return g()*x*ghkCaEffectiveCond();
}

// Get current
Number Order1CaEnergyBarrierTabChannel::Iion()
{
	Number x = value();

	return g()*x*ghkCaEffectivePotential();
}



// ====================================================================
// Order2CaEnergyBarrierTabChannel class body
// ====================================================================



// Constructors and destructor
Order2CaEnergyBarrierTabChannel::Order2CaEnergyBarrierTabChannel(Number gSpVal)
: EnergyBarrierTabChannel(gSpVal) {}

Order2CaEnergyBarrierTabChannel::~Order2CaEnergyBarrierTabChannel() {}

// Get conductance
Number Order2CaEnergyBarrierTabChannel::conductance()
{
	Number x = value();

	return g()*x*x*ghkCaEffectiveCond();
}

// Get current
Number Order2CaEnergyBarrierTabChannel::Iion()
{
	Number x = value();

	return g()*x*x*ghkCaEffectivePotential();
}



// ====================================================================
// CalciumPool class body
// ====================================================================



// Constructor and destructor
CalciumPool::CalciumPool() 
{
	// Clear the component name (use default)
	_componentName = NULL;
}

CalciumPool::~CalciumPool() 
{
	// Free any allocated name string
	delete _componentName;
}

// Get compartment
void CalciumPool::container(Compartment* comp)
{
	// Save new containing compartment
	_container = comp;
}

// Get the component name
const char* CalciumPool::componentName()
{
	return _componentName==NULL ? 
		"CaPool" : _componentName;
}

// Set the component name
void CalciumPool::componentName(char* name)
{
	_componentName = new char[strlen(name)+1];
	strcpy(_componentName,name);
}

// Add this to the compartment
void CalciumPool::addTo(Compartment* comp)
{
	comp->addCalciumPool(this);
}

// Remove this from the compartment
void CalciumPool::removeFromCompartment()
{
	if (container()!=NULL) {
		container()->removeCalciumPool(this);
	}
}

// Add a channel to those that supply ions
void CalciumPool::addSourceChannel(IonChannel* chan)
{
	_sourceChannels.push_back(chan);
}

// Remove a channel from those that supply ions
void CalciumPool::removeSourceChannel(IonChannel* chan)
{
	IonChannelVectorIt		last;

	last=remove(_sourceChannels.begin(),_sourceChannels.end(),chan);
	_sourceChannels.resize(last - _sourceChannels.begin());
}

// Find all calcium channels in the compartment and
// add them to this pool as calcium sources.
// This would typically be used during processing of
// startSimulation notification after the compartment
// contents are set up.
void CalciumPool::addAllCalciumChannels()
{
	IonChannelVector		chan( container()->getCalciumChannels() );
	IonChannelVectorIt		it;

	for (it=chan.begin();it!=chan.end();it++) {

		// Add the channel to those that contribute to the pool
		addSourceChannel(*it);

		// Hook channel to pool as the channels source of internal
		// calcium ion concentration.
		(*it)->calciumPool(this);
	}
}

// Sum up the current from all channels supplying Ca++ ions
Number CalciumPool::ICa()
{
	IonChannelVectorIt		it;
	Number					sum = 0;

	for (it=_sourceChannels.begin();it!=_sourceChannels.end();it++) {
		sum += (*it)->ICa();
	}
	return sum;
}



// ====================================================================
// SimpleCalciumPool class body
// ====================================================================



// Constructor and destructor
SimpleCalciumPool::SimpleCalciumPool() 
{
	using namespace UOM;

	// Since there is no compartment at this time, initialize
	// values to provide a known starting state.
	_container = NULL;
	_phi = 0;
	_beta = 0;

	// Set somewhat arbitrary default values just to get started.
	// Some of the values are taken from table 3 of
	// Jaffe et al. 1994, J. Neurophysiol. 71, 1065-1077.
	CaXrest			(Number( 50*nanoM ));
	CaXinit			(-1);						// no value set -- use CaXrest
	unboundRatio	(Number( 0.001 ));	
	vmax			(Number( 6e-14*mMole/msec/cm_2 ));
	Kd				(Number( 1*microM ));
	shellDepth		( numeric_limits<Number>::infinity() ); // no shell

	// Set default weighting for ODE error tolerance
	weight(Number( 1.0/(50*nanoM) ));	
}

SimpleCalciumPool::SimpleCalciumPool(
	Number				rest,		// resting concentration
	Number				init,		// initial concentration
	Number				ubr,		// unbound Ca++ ratio
	Number				vmx,		// vmax pump rate
	Number				kd,			// pump half rate concentration
	Number				shell)		// depth of shell 

{
	using namespace UOM;

	// Since there is no compartment at this time, initialize
	// values to provide a known starting state.
	_container = NULL;
	_phi = 0;
	_beta = 0;

	// Initialize with the values provided
	CaXrest			(rest);
	CaXinit			(init);
	unboundRatio	(ubr);	
	vmax			(vmx);
	Kd				(kd);
	shellDepth		(shell);

	// Set default weighting for ODE error tolerance
	weight(Number( 1.0/(100*nanoM) ));	
}

SimpleCalciumPool::SimpleCalciumPool(
	Number				givenPhi,		// fixed phi
	Number				givenBeta,		// fixed beta
	Number				givenRest)		// fixed CaXrest

{
	using namespace UOM;

	// There is no container so far, so set to NULL
	_container = NULL;
	_CaXinit = -1;

	// Set fixed values
	phi(givenPhi);
	beta(givenBeta);
	CaXrest(givenRest);

	// Set default weighting for ODE error tolerance
	weight(Number( 1.0/(100*nanoM) ));	
}

SimpleCalciumPool::~SimpleCalciumPool() {}

// Capture change in setting of container
void SimpleCalciumPool::container(Compartment* comp)
{
	// Let superclass do its version first
	CalciumPool::container(comp);

	// Update values computed from compartment sizes
	updateParams();
}

// Override any calculated value of phi
void SimpleCalciumPool::phi(Number x)
{
	_phi = x;
	_calcPhi = false;
	updateParams();
}

// Override any calculated value of beta
void SimpleCalciumPool::beta(Number x)
{
	_beta = x;
	_calcBeta = false;
	updateParams();
}

// Set phi based on a shell depth within a cylindrical compartment
void SimpleCalciumPool::shellDepth(Number d)
{
	_shellDepth = d;
	_calcPhi = true;
	updateParams();
}

// Set phi based on a constant ratio of unbound and bound forms of Ca
void SimpleCalciumPool::unboundRatio(Number ubr)
{
	_unboundRatio = ubr;
	_calcPhi = true;
	updateParams();
}

// Save the specific value for Vmax and force calculation
// of beta based on pump dynamics.
void SimpleCalciumPool::vmax(Number spvm)
{
	_vmax = spvm;
	_calcBeta = true;
	updateParams();
}

// Set a value for the pump half activation concentration
void SimpleCalciumPool::Kd(Number kd)
{
	_Kd = kd;
	_calcBeta = true;
	updateParams();
}

// Set a value for resting calcium concentration
void SimpleCalciumPool::CaXrest(Number x)
{
	_CaXrest = x;
	updateParams();
}

// Rescale as needed based on compartment size.
void SimpleCalciumPool::updateParams()
{
	// Only do this is there is an associated compartment.
	// This function is invoked when a compartment size
	// changes or when the compartment is initially associated
	// with the pool.

	if (container()!=NULL ) {

		const Number z = 2;	// Ion charge

		Number	area = container()->membraneArea() / container()->areaAdjustment();
		Number	poolVol = container()->subshellVolume( shellDepth() );

		// Update phi and beta as needed

		if (_calcPhi) {
			_phi = unboundRatio()/(z*UOM::Faraday*poolVol);
		}

		if (_calcBeta) {
			// Assign a beta based on scaled vmax and Kd
			// This is scaled to be in the same range
			// and fit into the formula applicable with a
			// non-saturating pump. Note that vmax/kd is the
			// extrusion rate per area unit at [Ca++]-i = 0.
			_beta = area / poolVol * vmax() / Kd();
			
			// Compute a resting concentration that results
			// in the derivative being zero at CaXrest.
			_prest = CaXrest()/(1+CaXrest()/Kd());
		}
		else {
			_prest = CaXrest();
		}
	}
}

// Compute derivatives based on simple diffusion model
void SimpleCalciumPool::computeDerivatives()
{
	Number cax = stateValue(0);		// avoid extra call to CaX
	Number pca;						// effective concentration

	// If a value for beta was calculated, incorporate pump
	// dynamics into effective concentrations.
	if (_calcBeta ) {
		pca = cax/(1+cax/Kd());
	}
	else {
		pca = cax;
	}

	// Note that ICa is an inward current and therefore has
	// a negative value. Hence the leading minus sign below.

	derivValue(0)= -phi()*ICa() - beta()*(pca - _prest);
}

// Make an implicit update to the state advancing by time h.
// An implicit update rule is used.
void SimpleCalciumPool::localStateUpdate(SimTime h, CNStepType stepType)
{

	// This logic closely follows that of compute derivatives.
	// If the local ODE is dy/dt=f(y), the local update is:
	// y(t+h)=y(t)+h/(1-h/2*df/dy)*f(y(t))

	// Even though this update method is second order accurate,
	// the observed effect is that calcium concentrations are
	// only first order accurate when voltage and channel states
	// change rapidly.

	Number cax = stateValue(0);		// avoid extra call to CaX
	Number pca;						// effective concentration
	Number ydot;					// f(y) in above
	Number dfdy;					// df/dy in above
	Number pdd;						// pump dynamics denominator below

	// If a value for beta was calculated, incorporate pump
	// dynamics into effective concentrations.
	if (_calcBeta) {
		pdd = 1+cax/Kd();
		pca = cax/pdd;
		dfdy = -beta()/(pdd*pdd);	// deriv of f wrt cax
	}
	else {
		pca = cax;
		dfdy = -beta();				// deriv of f wrt cax
	}

	// Note that ICa is typically an inward current
	// and therefore has a negative value. Hence the
	// leading minus sign below.
	ydot = -phi()*ICa() - beta()*(pca - _prest);

	// Apply the semi-implicit second order update rule.
	// where the ending state is the midpoint (y(t)+y(t+2h))/2.
	if (stepType == CNFullStep) {
		stateValue(0) += h/(1-h/2*dfdy)*ydot;
	}
	else {
		stateValue(0) += h/(1-h*dfdy)*ydot;
	}
}

// Get calcium channels at simulation start
void SimpleCalciumPool::simulationStarted()
{
	// See if the source channels have been identified.
	// If not, add all calcium channels in the compartment.
	if (_sourceChannels.empty()) {
		addAllCalciumChannels();
	}

	// Initialize the state vector early so that
	// channels can access calcium concentrations
	// during their setInitialState processing
	setInitialState();
}

// Set initial state values
void SimpleCalciumPool::setInitialState()
{
	// Initialize to the resting potential.
	// Use the value in CaXinit if there is one.
	// Otherwise, start at nominal resting potential.
	if (CaXinit() >= 0) {
		stateValue(0) = CaXinit();
	}
	else {
		stateValue(0) = CaXrest();
	}
}

// Set weight values somewhat arbitrarily based on typical
// peak values for pyramidal cell dendrites
void SimpleCalciumPool::setWeightValues()
{
	weightValue(0)= weight();
}



// ====================================================================
// SynapticResponse class body
// ====================================================================



// Initialize static data
int SynapticResponse::_LateArrivingEventCount = 0;
SimTime SynapticResponse::_TotalLateEventTime = 0;

// Constructors and destructor
SynapticResponse::SynapticResponse()
{
	// Initialize data elements
	_groupOwner = NULL;
	_synapses = NULL;
	_synapseCount = 0;
	_synapseStateOffset = 0;
}

SynapticResponse::~SynapticResponse()
{
	// See if  this object owns a list of synapses or not
	if ( !isGroupMember() ) {
		// Remove all postsynaptic data
		destroyAllSynapses();
	}
	else {
		// Remove this response from the associated group
		groupOwner()->remove(this);
	}
}

// Set the group owner while making sure that there are
// no synapses currently associated with this object
void SynapticResponse::groupOwner(SynapticGroupResponse* owner)
{
	// Check for a non-empty synapse list
	if (isActive() && owner!=NULL) {
		FatalError("(SynapticResponse::groupOwner) Synapse list is non-empty while setting new owner.");
	}

	// Set the new owner
	_groupOwner = owner;
}

// Access synapses via group owner
Synapse* SynapticResponse::groupSynapses()
{
	return _groupOwner->synapses();
}

// Access synapse count via group owner
int SynapticResponse::groupSynapseCount()
{
	return _groupOwner->synapseCount();
}

// Access AP queue via group owner
ActionPotentialEventQueue& SynapticResponse::groupAPQueue()
{
	return _groupOwner->apQueue();
}

// Create a synapse to go with this channel.
// Memory is allocated as a block into which a Synapse object is placed
// followed by an instance of an appropriate subclass of SynapseState, 
// the size of which is determined by the implementing subclass.
Synapse* SynapticResponse::createSynapse(
	AxonProcess*			axon,			// presynaptic axon process
	Number					wght,			// initial weight value
	Number					dist)			// axonal distance
{
	const unsigned int		synapseSize		= sizeInBytesRounded(sizeof(Synapse));
	const unsigned int		stateSize		= synapseStateSize();

	unsigned int			allocatedSize;
	Byte*					allocatedMemory;
	Synapse*				allocatedSynapse;

	// This cannot be done by this object if a member of a group.
	// Only the group owner is permitted in that case.
	if ( isGroupMember() ) {
		FatalError("(SynapticResponse::createSynapse) Cannot create synapse for group member.");
	}

	// Determine the size of the memory to allocate
	allocatedSize = synapseSize + stateSize;

	// Allocate the memory. Note that the std template library provides
	// coverage of allocation errors by raising an exception.
	allocatedMemory = new Byte[allocatedSize];

	// Even though it is trivial in this case, fill in offset for state data.
	// For subclasses with multiple responses this is not so simple.
	// An assumption here is that all synapses are create equal and the same.
	synapseStateOffset(synapseSize);

	// Fill in the synapse data. This variation on new is implemented
	// as part of the standard template library. It just runs the constructor.
	allocatedSynapse = new(reinterpret_cast<void*>(allocatedMemory)) Synapse(axon,this,dist);
	addSynapse(allocatedSynapse);

	// Allow subclass to fill in its part using the already known offset.
	// The subclass should fill in the weight as needed.
	createSynapseState(reinterpret_cast<Synapse*>(allocatedMemory),wght);

	// Return the synapse created
	return allocatedSynapse;
}

// Destroy a synapse by invoking associated state destructors 
// and unhooking it from the postsynaptic list. When both presynaptic
// and postsynaptic sides are unhooked, the synapse self destructs.
void SynapticResponse::destroySynapse(Synapse* syn)
{
	// This cannot be done by this object if a member of a group.
	// Only the group owner is permitted in that case.
	if (isGroupMember()) {
		FatalError("(SynapticResponse::destroySynapse) Cannot do this while member of a group.");
	}

	// Destroy synapse state object(s), if any.
	destroySynapseState(syn);

	// Unhook the synapse from the postsynaptic list.
	removeSynapse(syn);
}

// Create synapse state data. This default is for a SynapseState object only,
// but a subclass can override to create a more complex state object.
// This general pattern can be followed, but since the class of the state object
// unknown here, this code cannot be directly extended.
void SynapticResponse::createSynapseState(Synapse* syn, Number wght)
{
	// Locate the state date
	SynapseState*	state = reinterpret_cast<SynapseState*>
		(reinterpret_cast<Byte*>(syn)+synapseStateOffset());

	// Set the associated weight value. This variation on new just runs
	// the constructor for an object at the address given.
	new(reinterpret_cast<void*>(state)) SynapseState;

	// Set the initial weight
	state->weight=wght;
}

// Add an associated synapse
void SynapticResponse::addSynapse(Synapse* syn)
{
	// Make sure all adds/deletes are via common owner
	if (isGroupMember() ) {
		FatalError("(SynapticResponse::addSynapse) " 
			"Group member cannot add a synapse");
	}

	// Add to the beginning of the list and adjust count
	syn->addBeforeInPostsynapticList(_synapses);
	_synapses = syn;
	_synapseCount++;

	// Indicate that a synapse was added
	synapseAdded(syn);
}

// Remove an associated synapse
void SynapticResponse::removeSynapse(Synapse* syn)
{
	Synapse*		next;

	// Make sure all adds/deletes are via common owner
	if (isGroupMember() ) {
		FatalError("(SynapticResponse::removeSynapse) " 
			"Group member cannot remove a synapse");
	}

	// See if the synapse to be removed is the first one
	if (_synapses==syn) {
		_synapses = syn->nextPostsynaptic();
	}
	else {
		// Scan the list looking for the synapse before the one to be removed
		// and adjust the next link to skip the deleted synapse.
		next = _synapses;
		while (next!=NULL && next->nextPostsynaptic()!=syn) {
			next=next->nextPostsynaptic();
		}
		if (next!=NULL) {
			next->removeNextFromPostsynapticList();
		}
		else {
			FatalError("(SynapticResponse::removeSynapse) Synapse to remove not found");
		}
	}

	// Clear the presynaptic link in the synapse removed
	syn->clearPostsynaptic();
	_synapseCount--;

	// Indicate that a synapse was removed
	synapseRemoved(syn);
}

// Destroy all current synapses associated with this object
void SynapticResponse::destroyAllSynapses()
{
	Synapse*	syn;
	Synapse*	next;

	syn=_synapses;
	while (syn!=NULL) {

		// Get the next pointer before destroying the one in hand.
		next = syn->nextPostsynaptic();

		// Destroy this synapse and any state data with it.
		destroySynapse(syn);

		// Move on to the next on the list
		syn=next;
	}
}

// Add a new action potential to the queue
void SynapticResponse::signalActionPotential(
	SimTime presynTime, 
	ActionPotentialEvent* apEvent)
{
	ActionPotentialEvent apev = *apEvent; // make private copy

	// Update the event time to reflect delays between AP arrival and
	// release of neurotransmitter for this type of synapse.
	apev.eventTime(presynTime+synapticDelay());

	// See if the action potential has arrived late.
	// If so allow for special handling to catch up.
	if (apev.eventTime() < currentTime() ) {

		// Keep count of late arriving events
		_LateArrivingEventCount++;
		_TotalLateEventTime += currentTime()-apev.eventTime();

		// Handle a late arrival (subclass may customize)
		handleLateAPEvent(&apev);
	}

	// Queue the event for later processing when its time comes.
	addAPEventToQueue(&apev);
}

// Add an AP event to the queue. Subclass may override to customize.
void SynapticResponse::addAPEventToQueue(ActionPotentialEvent* apEvent)
{
	apQueue().add(apEvent);
}

// Handle start of simulation by initializing (in subclass)
// any relevant class cache objects. Note that every component
// instance gets this individually and must check for cases
// where initialization is done already.
void SynapticResponse::simulationStarted()
{
	// Let superclass do its initialization
	IonChannel::simulationStarted();

	// Force initialization of any class cache (via subclass)
	initializeClassCache();
}

// When the time step is over, update state as needed
void SynapticResponse::timeStepEnded()
{
	SimTime					now = currentTime();
	ActionPotentialEvent*	apPtr;

	// Let a subclass do its own end to step processing
	applyEndOfTimeStep();

	// Do AP purge if this object owns the queue
	if (!isGroupMember() ) {

		// Purge AP queue entries up to the current time
		for (apPtr=_apQueue.begin();
			apPtr!=_apQueue.end() && apPtr->eventTime()<now; 
			apPtr=_apQueue.next(apPtr) ) {

			_apQueue.removeFirst();
		}
	}
}


// Static accessors
int SynapticResponse::lateArrivingEventCount()
{
	return _LateArrivingEventCount;
}

SimTime SynapticResponse::meanLateEventTime()
{
	return (_LateArrivingEventCount > 0) ?
		_TotalLateEventTime / _LateArrivingEventCount : 0;
}

void SynapticResponse::resetLateArrivingEventStatistics()
{
	_LateArrivingEventCount = 0;
	_TotalLateEventTime = 0;
}

// Access the size of the default synapse weight instance
unsigned int SynapticResponse::synapseStateSize() 
{ 
	return sizeInBytesRounded(sizeof(SynapseState)); 
}	

// Get the weight value for a synapse for this resposne
Number SynapticResponse::synapseWeight(Synapse* syn)
{
	return (reinterpret_cast<SynapseState*>
		(reinterpret_cast<Byte*>(syn)+synapseStateOffset()))->weight;
}

// Set the weight value for a synapse for this resposne
void SynapticResponse::synapseWeight(Synapse* syn, Number w)
{
	(reinterpret_cast<SynapseState*>
	(reinterpret_cast<Byte*>(syn)+synapseStateOffset()))->weight = w;
}

// Return the total of all synapses weights for this response type
Number SynapticResponse::totalSynapticWeight(int* pCount)
{
	Synapse*	syn = synapses();
	Number		totalWeight = 0;
	int			count = 0;

	// Sum up total weight
	while (syn!=NULL) {
		totalWeight += 
			(reinterpret_cast<SynapseState*>
			(reinterpret_cast<Byte*>(syn)+synapseStateOffset()))->weight;
		count++;
		syn = syn->nextPostsynaptic();
	}

	// If requested, return the count also
	if (pCount!=NULL) {
		*pCount = count;
	}

	return totalWeight;
}



// ====================================================================
// SynapticGroupResponse class body
// ====================================================================



// Constructors and destructor
SynapticGroupResponse::SynapticGroupResponse() {}

SynapticGroupResponse::~SynapticGroupResponse() 
{
	SynapticResponseVectorIt it;

	// Delete any associated members after explicitly unhooking them
	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->groupOwner(NULL);
		delete *it;
	}
}

// Add a member to the list of responses in the group
void SynapticGroupResponse::addResponse(SynapticResponse* resp) 
{
	// Once synapses exist, then changes can no longer be made
	if (isActive() ) {
		FatalError("(SynapticGroupResponse::addResponse) "
			"Updates cannot be done once synapses are allocated.");
	}

	// Add the response and set ownership
	_members.push_back(resp);
	resp->groupOwner(this);

	// Set model and containing compartment in new member
	resp->model( model() );
	resp->container( container() );

}

// Remove a member from the list of responses in the group
void SynapticGroupResponse::removeResponse(SynapticResponse* resp) 
{
	SynapticResponseVectorIt last;

	// Once synapses exist, then changes can no longer be made
	if (isActive() ) {
		FatalError("(SynapticGroupResponse::removeResponse) "
			"Updates cannot be done once synapses are allocated.");
	}

	// Remove the ion channel here
	last=std::remove(_members.begin(),_members.end(),resp);
	_members.resize( last - _members.begin() );

	// Clear model and container in member
	resp->model( NULL );
	resp->container( NULL );

	// Clear group ownership
	resp->groupOwner( NULL );
}

// Set model and pass along to members
void SynapticGroupResponse::model(Model* m)
{
	SynapticResponseVectorIt it;

	// Set model for this object
	SynapticResponse::model(m);

	// Pass along to members
	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->model(m);
	}
}

// Set container and pass along to members
// Need to also maintain relationship with the
// containing compartment for group members.
void SynapticGroupResponse::container(Compartment* comp)
{
	SynapticResponseVectorIt it;

	// Set the containing compartment for this object
	SynapticResponse::container(comp);

	// Pass along to members
	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->container( comp );
	}
}

// Provide the size to allocate for state data.
unsigned int SynapticGroupResponse::synapseStateSize()
{
	SynapticResponseVectorIt it;

	unsigned int totalSize = 0;

	// Total sizes for needed for all members
	for (it=_members.begin(); it!=_members.end(); it++) {
		totalSize += (*it)->synapseStateSize();
	}

	return totalSize;
}

// Apply a constructor to synapse state data.
void SynapticGroupResponse::createSynapseState(Synapse* syn, Number wght)
{
	SynapticResponseVectorIt it;
	int	offset = synapseStateOffset();
	Number w =wght;

	// Have each member create its own state data
	// adjusting offsets along the way for allocated states.
	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->synapseStateOffset(offset);
		(*it)->createSynapseState(syn,w);
		offset += (*it)->synapseStateSize();

		// See if same weight is applied to all members
		if (!applyInitialWeightToAll() ) {
			w=1;
		}
	}
}

// Apply a destructor to synapse state data.
void SynapticGroupResponse::destroySynapseState(Synapse* syn)
{
	SynapticResponseVectorIt it;

	// Have each member destory its own state data
	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->destroySynapseState(syn);
	}
}

// Pass changes in the synapse list along to members
void SynapticGroupResponse::synapseAdded(Synapse* syn)
{
	SynapticResponseVectorIt	it;

	SynapticResponse::synapseAdded(syn);

	for (it=_members.begin();it!=_members.end();it++) {
		(*it)->synapseAdded(syn);
	}
}

// Pass changes in the synapse list along to members
void SynapticGroupResponse::synapseRemoved(Synapse* syn)
{
	SynapticResponseVectorIt	it;

	SynapticResponse::synapseRemoved(syn);

	for (it=_members.begin();it!=_members.end();it++) {
		(*it)->synapseRemoved(syn);
	}
}

// Notify members of a late arriving action potential event
void SynapticGroupResponse::handleLateAPEvent(ActionPotentialEvent* apEvent)
{
	SynapticResponseVectorIt	it;

	for (it=_members.begin();it!=_members.end();it++) {
		(*it)->handleLateAPEvent(apEvent);
	}
}

// Notify members when time step is over
void SynapticGroupResponse::applyEndOfTimeStep()
{
	SynapticResponseVectorIt	it;

	for (it=_members.begin();it!=_members.end();it++) {
		(*it)->applyEndOfTimeStep();
	}
}

// Response delay for associated synapses.
// By default, value is taken from the first response member.
SimTime SynapticGroupResponse::synapticDelay()
{
	return _members[0]->synapticDelay();
}

// Return total conductance for the group.
Number SynapticGroupResponse::conductance()
{
	SynapticResponseVectorIt it;
	Number sum = 0;

	for (it=_members.begin(); it!=_members.end(); it++) {
		sum += (*it)->conductance();
	}

	return sum;
}

// Return total ionic current for the group.
Number SynapticGroupResponse::Iion()
{
	SynapticResponseVectorIt it;
	Number sum = 0;

	for (it=_members.begin(); it!=_members.end(); it++) {
		sum += (*it)->Iion();
	}

	return sum;
}

// Return total calcium current for the group.
Number SynapticGroupResponse::ICa()
{
	SynapticResponseVectorIt it;
	Number sum = 0;

	for (it=_members.begin(); it!=_members.end(); it++) {
		sum += (*it)->ICa();
	}

	return sum;
}

// Return both conductance and current at once
// by summing over group members.
void SynapticGroupResponse::condAndIion(
	Number&	Gout,			// conductance
	Number& Iout)			// current
{
	// Reset conductance and current
	Gout = 0;
	Iout = 0;

	// Stop now if there are no synapses
	if (isInactive() )
		return;

	// Total conductance and current over group
	SynapticResponseVectorIt it;

	for (it=_members.begin(); it!=_members.end(); it++) {

		Number	Gtemp, Itemp;

		(*it)->condAndIion(Gtemp,Itemp);
		Gout += Gtemp;
		Iout += Itemp;
	}
}

// Add to components to probe when reporting.
// This includes the current object and any member responses.
void SynapticGroupResponse::addToComponentsToProbe(ModelComponentVector& comps)
{

	int i;

	// Add this object
	comps.push_back(this);

	// Add all members
	for (i=0;i<_members.size();i++) {
		_members[i]->addToComponentsToProbe(comps);
	}
}

// Locate the identified component within this ion channel.
// If the component is not found, return NULL.
IonChannel* SynapticGroupResponse::findIonChannel(TokenId compId)
{
	SynapticResponseVectorIt it;


	// Check this object for a match
	if (componentId()==compId)
		return this;

	// Check members and return the first found
	for (it=_members.begin(); it!=_members.end(); it++) {
		if ( (*it)->componentId()==compId) 
			return *it;
	}

	// Nothing found here
	return NULL;
}

// Pass any neuromodulation parameters to members
void SynapticGroupResponse::setModParams(TokenId modId,int numParams,Number* params)
{
	SynapticResponseVectorIt it;

	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->setModParams(modId,numParams,params);
	}
}

// For testing, disable plasticity in all members
void SynapticGroupResponse::disablePlasticity()
{
	SynapticResponseVectorIt it;

	for (it=_members.begin(); it!=_members.end(); it++) {
		(*it)->disablePlasticity();
	}
}



// ====================================================================
// SynapticConductance class body
// ====================================================================



// Construct a new instance
SynapticConductance::SynapticConductance()
{

	// Indicate that there are no plasticity rules
	_presynapticRule = NULL;
	_postsynapticRule = NULL;
}

// Destroy the current instance
SynapticConductance::~SynapticConductance()
{
	// Clear the list of synapses.
	// Note that we have to do this here even though exactly
	// the same function is invoked in SynapticResponse just
	// after this because during the SynapticResponse destructor
	// processing function overrides are no longer applied.
	destroyAllSynapses();

	// Discard any plasticity rules
	if (_presynapticRule != NULL) {
		delete _presynapticRule;
	}
	if (_postsynapticRule != NULL) {
		delete _postsynapticRule;
	}
}

// Indicate an error if gSpecific is invoked
void SynapticConductance::gSpecific(Number gval)
{
	FatalError("(SynapticConductance::gSpecific) Function not supported for this class.");
}

// Indicate an error if gAbsolute is invoked
void SynapticConductance::gAbsolute(Number gval)
{
	FatalError("(SynapticConductance::gAbsolute) Function not supported for this class.");
}

// Set the model and pass along to any plasticity rule
void SynapticConductance::model(Model* m)
{
	// Let superclass hook up this object with model
	SynapticResponse::model(m);

	// Inform any plasticity rules
	if (_presynapticRule != NULL) {
		_presynapticRule->model(m);
	}
	if (_postsynapticRule != NULL) {
		_postsynapticRule->model(m);
	}
}

// Set the presynaptic rule
void SynapticConductance::presynapticRule(PlasticityRule* pr)
{
	// If this is a no-op we can stop now
	if (pr==presynapticRule() )
		return;

	// If there are already synapses, plasticity rule cannot be changed
	if ( isActive() ) {
		FatalError("(SynapticConductance::presynapticRule) "
			"Cannot make changes when synapses already exist");
	}

	// Otherwise, the rule can still be changed,
	// but first get rid of any old rule.
	if (_presynapticRule != NULL) {
		delete _presynapticRule;
	}

	// Hook up with the new rule
	if (pr!=NULL) {
		pr->synapticCond(this);
	}
	_presynapticRule = pr;

	// Force recalculation of offsets
	synapseStateOffset( synapseStateOffset() );
}

// Set the postsynaptic rule
void SynapticConductance::postsynapticRule(PlasticityRule* pr)
{
	// If this is a no-op we can stop now
	if (pr==postsynapticRule() )
		return;

	// If there are already synapses, plasticity rule cannot be changed
	if ( isActive() ) {
		FatalError("(SynapticConductance::plasticityRule) "
			"Cannot make changes when synapses already exist");
		return;
	}

	// Otherwise, the rule can still be changed,
	// but first get rid of any old rule.
	if (_postsynapticRule != NULL) {
		delete _postsynapticRule;
	}

	// Hook up with the new rule
	if (pr!=NULL) {
		pr->synapticCond(this);
	}
	_postsynapticRule = pr;

	// Force recalculation of offsets
	synapseStateOffset( synapseStateOffset() );
}

// Unhook from either a presynaptic or postsynaptic plasticity rule
void SynapticConductance::unhookFromRule(PlasticityRule* pr)
{
	if (presynapticRule()==pr)		_presynapticRule = NULL;
	if (postsynapticRule()==pr)		_postsynapticRule = NULL;
}

// As a testing convenience, unhook and delete plasticity rules.
// Any plasticity data is abandoned in place.
void SynapticConductance::disablePlasticity()
{
	delete _presynapticRule;
	_presynapticRule = NULL;

	delete _postsynapticRule;
	_postsynapticRule = NULL;
}

// Apply any neuromodulation parameters
void SynapticConductance::setModParams(TokenId modId,int numParams,Number* params)
{
	// Defer basic response to superclass
	SynapticResponse::setModParams(modId,numParams,params);

	// Inform any plasticity rules
	if (presynapticRule()!=NULL) {
		presynapticRule()->setModParams(modId,numParams,params);
	}
	if (postsynapticRule()!=NULL) {
		postsynapticRule()->setModParams(modId,numParams,params);
	}
}

// Handle changes in synapse list including via a group owner.
void SynapticConductance::synapseRemoved(Synapse* syn)
{
	// Clear the state if there are no synapses left
	if (isInactive() ) {
		clearState();
	}
}

// Set the offset for synapse state data plus
void SynapticConductance::synapseStateOffset(unsigned int nbytes)
{
	// Set the offsets to allow for receptor specific state, 
	// plasticity state, and presynaptic state in that order.

	unsigned int		presynapticStateOffset;
	unsigned int		postsynapticStateOffset;

	_synapseStateOffset = nbytes;

	presynapticStateOffset = nbytes + receptorSpecificStateSize();
	if (presynapticRule()!=NULL) {
		presynapticRule()->stateOffset(presynapticStateOffset);
	}

	postsynapticStateOffset = presynapticStateOffset;
	if (presynapticRule()!=NULL) {
		postsynapticStateOffset += presynapticRule()->plasticityStateSize();
	}
	if (postsynapticRule()!=NULL) {
		postsynapticRule()->stateOffset(postsynapticStateOffset);
	}
}

// Return the size to allocate for state data including any needed
// to support a plasticity rule.
unsigned int SynapticConductance::synapseStateSize()
{
	unsigned int nbytes = receptorSpecificStateSize();

	if (presynapticRule()!=NULL) {
		nbytes += presynapticRule()->plasticityStateSize();
	}
	if (postsynapticRule()!=NULL) {
		nbytes += postsynapticRule()->plasticityStateSize();
	}

	return nbytes;
}

// Apply a constructor to synapse state data.
// This creates other state data plus any plasticity state.
void SynapticConductance::createSynapseState(Synapse* syn, Number wght)
{
	// Create any other state data
	createReceptorSpecificState(syn,wght);

	// Create any plasticity state data
	if (presynapticRule() != NULL) {
		presynapticRule()->createPlasticityState(syn,wght);
	}
	if (postsynapticRule() != NULL) {
		postsynapticRule()->createPlasticityState(syn,wght);
	}
}

// Apply a destructor to any plasticity and synapse state data.
void SynapticConductance::destroySynapseState(Synapse* syn)
{
	// Destroy the plasticity state data
	if (postsynapticRule() != NULL) {
		postsynapticRule()->destroyPlasticityState(syn);
	}
	if (presynapticRule() != NULL) {
		presynapticRule()->destroyPlasticityState(syn);
	}

	// Destroy "other" state data
	destroyReceptorSpecificState(syn);
}

// Return the size to allocate for other state data.
// By default this is the size of SynapseState.
unsigned int SynapticConductance::receptorSpecificStateSize()
{
	// Use the superclass to get size
	return SynapticResponse::synapseStateSize();
}

// Create a non-plasticity state object in a synapse buffer.
// By default this creates an instance of SynapseState.
void SynapticConductance::createReceptorSpecificState(Synapse* syn, Number wght)
{
	// Use superclass to do the create of a SynapseState object
	SynapticResponse::createSynapseState(syn,wght);
}

// When the time step is over, update state as needed
void SynapticConductance::applyEndOfTimeStep()
{
	// Update state based on the passage of time
	updateStateForTimeStep();

	// Notify any plasticity rule that the time step is ended 
	// (before AP purge). Rules that are not currently enabled 
	// are not notified unless so indicated by the rule subclass.
	if (presynapticRule()!=NULL && presynapticRule()->notifyRule() ) {
		presynapticRule()->applyEndOfStep(apQueue());
	}
	if (postsynapticRule()!=NULL && postsynapticRule()->notifyRule() ) {
		postsynapticRule()->applyEndOfStep(apQueue());
	}
}

// Default implementation to update state to reflect AP Events.
// Interval for events to included is [fromTime toTime)
void SynapticConductance::updateForQueuedAPEvents(
	SimTime from, 
	SimTime to)
{
	ActionPotentialEventQueue&	apQ = apQueue();
	ActionPotentialEvent*		apPtr;

	// Process any events on the queue up to the current time.
	for (apPtr = apQ.begin(); 
		apPtr != apQ.end() && apPtr->eventTime()<to; 
		apPtr  = apQ.next(apPtr) ) {

		// If the action potential lies in the current interval,
		// add its impact to the current working state.
		if (apPtr->eventTime()>=from) {

			// Do any finalization processing when an event is
			// first processed. This allows finalization to be
			// done in time sequence when current state is available.
			// Obviously this can only be done once, especially if
			// the outcome is stochastic. Do not notify the rule if
			// it currently disabled unless overridden by rule subclass.
			if (!apPtr->isFinal() && presynapticRule()!=NULL
				&& presynapticRule()->notifyRule() ) {
				presynapticRule()->finalizeAPEvent(apPtr);
			}

			// Update state to reflect the event.
			updateForAPEvent(apPtr);
		}
	}
}

// Update the receptor state to reflect a single AP Event.
void SynapticConductance::updateForAPEvent(ActionPotentialEvent* apEvent)
{
	FatalError("(SynapticConductance::updateForAPEvent) Subclass must override");
}


// Return a weight adjusted to compensate for spine neck resistance.
// For purposes of the adjustment, the response is assumed to be
// proportional to synapse weight and AP quantity.
Number SynapticConductance::adjustedSynapticWeight(ActionPotentialEvent* apEvent)
{

	Number		r;			// spine neck resistance
	Number		w;			// current synapse weight
	Number		q;			// current response quantity

	w = synapseWeight(apEvent->synapse());

	// Skip the arithmetic if there is no spine neck
	if ( (r=spineNeckResistance()) == 0) {
		return w;
	}

	// Adjust for spine neck resistance by treating
	// the spine as a unipotential ball with a neck
	// formed of a resister. This adjustment approximates
	// the current integral by assuming the average receptor
	// conductance is half the peak value. Peak value is
	// in turn assumed to be proportional to weight and
	// quantity.

	q = apEvent->quantity();
	return 2*w/(2 + w*q*r*gMax());
}

// Utility to return the peak value for a dual exponent 
// formulation: tau2/(tau2-tau1)*(exp(-t/tau2)-exp(-t/tau1))
// Where tau1==tau2, an alpha function formula is used.
Number SynapticConductance::peakDualExpResp(SimTime tau1, SimTime tau2)
{
	Number		speak;
	SimTime		tpeak;

	// See which case applies - alpha function or not
	if (fabs((tau1-tau2)/tau1)<1e-4) {

		// Alpha function case (note scaling by 1/tau1)
		speak = 1/E;
	}
	else {
		// Dual exponent case
		tpeak = log(tau2/tau1)/(1/tau1 - 1/tau2);
		speak = tau2/(tau2-tau1)*(exp(-tpeak/tau2) - exp(-tpeak/tau1));
	}

	return speak;
}

// Utility to return the peak value for a triple exponent formulation.
Number SynapticConductance::peakTripleExpResp(
	SimTime				tau1, 
	SimTime				tau2, 
	SimTime				tau3,
	double				c,
	double				rtol)

{
	// Formula to be maximized is the solution of the underlying ODE.
	// This solution is the superpostion of two dual exponent systems.
	//
	// f(t) = c*tau2/(tau2-tau1)*(exp(-t/tau2)-exp(-t/tau1)) +
	//        (1-c)*tau3/(tau3-tau1)*(exp(-t/tau3)-exp(-t/tau1))
	//
	// where tau1<tau2, tau1<tau3, and 0<=c<=1.

	// Where tau1=tau2 or tau1=tau3, an alpha function limit is used.

	// Note that the alpha function limit is t*exp(-t/tau)/tau because
	// of the way the ODE state variables are scaled (see classes
	// DualExpSynapticCond and TripleExpSynapticCond for a discussion).

	int			niter=0;			// number of iterations
	int			maxIter = 32;		// maximum iterations allowed
	double		f;					// function values
	double		f2,f3;				// temp parts of f relating to tau2 and tau3
	double		t,tL,tH;			// time value and low-high range
	double		t2,t3;				// peak times for tau2 and tau3 components
	double		df;					// function derivative at t
	double		df2,df3;			// parts of df relating to tau2 and tau3
	double		e1,e2,e3;			// exponent values
	double		k2,k3;				// constant multipliers
	double		tol;				// peak time error tolerance

	bool		alpha2 = fabs((tau1-tau2)/tau1)<1e-4;
	bool		alpha3 = fabs((tau1-tau3)/tau1)<1e-4;

	// Precompute useful constants (when alpha function does not apply)
	if (!alpha2)	k2 = tau2/(tau2-tau1);
	if (!alpha3)	k3 = tau3/(tau3-tau1);

	// Get an error tolerance in time based on the fastest time constant.
	tol = rtol*tau1;

	// Locate the bracketing range for t taking advantage of the fact that
	// the function is 0 at t=0 and at t=infinity. Since f can be rewritten
	// as the sum of two dual exponential parts each of which is unimodal,
	// a bounding range can be found from the respective dual exponent parts.

	t2 = alpha2 ? tau1 : log(tau2/tau1)/(1/tau1 - 1/tau2);
	t3 = alpha3 ? tau1 : log(tau3/tau1)/(1/tau1 - 1/tau3);

	tL = minval(t2,t3);
	tH = maxval(t2,t3);	

	// Do a simple line search looking for a zero derivative.
	// Even if the exact maximum point is not found, the peak value
	// should be a good approximation.
	while (tH-tL>tol && niter++<maxIter) {
		
		// Locate the mid point for a bisection.
		t=(tL+tH)/2;

		// Compute the derivative at t
		e1 = exp(-t/tau1);
		e2 = exp(-t/tau2);
		e3 = exp(-t/tau3);
		
		df2 = alpha2 ? (1-t/tau1)/tau1*e1 : k2*(e1/tau1-e2/tau2);
		df3 = alpha3 ? (1-t/tau1)/tau1*e1 : k3*(e1/tau1-e3/tau3);
		df = c*df2+(1-c)*df3;

		// Adjust the range so that df(tL)>0 and df(tH)<0
		if (df>0)		tL=t;
		else			tH=t;
	}

	// Get the value of f at t and return it as the peak value
	f2 = alpha2 ? t/tau1*exp(-t/tau1) : k2*(exp(-t/tau2)-exp(-t/tau1));
	f3 = alpha3 ? t/tau1*exp(-t/tau1) : k3*(exp(-t/tau3)-exp(-t/tau1));
	f = c*f2+(1-c)*f3;

	return f;
}


// ====================================================================
// SingleExpSynapticCond class body
// ====================================================================



// Construct a new instance from a maximum conductance
SingleExpSynapticCond::SingleExpSynapticCond(Number gMaxValue) 
{
	// Save the gMax value now so that it can be accessed
	// during setInitialState which invokes gMax.
	// Note that at this stage of construction, the 
	// gMax(Number x) accessor cannot be used because
	// time constants are currently unknown.
	_gMax = gMaxValue;

	// Clear any cached state at the outset
	_cachedTime = InfiniteFuture;
	_s1Start = 0;
	_s1Now = 0;
}

// Destroy the current instance
SingleExpSynapticCond::~SingleExpSynapticCond() {}

// Set initial state values plus cache values at startup
void SingleExpSynapticCond::setInitialState()
{
	// Set state for superclass
	SynapticConductance::setInitialState();

	// Initialize the class cache
	initializeClassCache();

	// Update values derived from parameters
	gMax(_gMax);

	// Set starting state values to zero
	clearState();
}

// Initialize the class cache. Subclasses may need
// to extend this with their own initializations.
void SingleExpSynapticCond::initializeClassCache()
{
	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	// Stop now if already initialized
	if (pCC->isInitialized) return;

	// Set initial values
	pCC->Tau1 = tau1()/Q10Factor();
	pCC->H = -1; // i.e. no value yet

	// Indicate that cache has been initialized
	pCC->isInitialized = true;
}

// Set conductance based on maximum value attainable
// This is likely to be overridden in subclasses where
// setting the conductance value of is not so degenerate.
void SingleExpSynapticCond::gMax(Number gm)
{
	// Set the scaled maximum conductance value.
	// Of course, in this case, no scaling is needed.
	SynapticResponse::gAbsolute(gm);

	// Save the gMax value in case it is needed later.
	_gMax = gm;

	// Clear any cached state dependent on gMax
	emptyCaches();
}

// Get conductance and current together
void SingleExpSynapticCond::condAndIion(
	Number&	Gout,			// conductance
	Number& Iout)			// current
{
	// Take fast path if there are no synapses
	if (isInactive() ) {
		Gout = 0;
		Iout = 0;
	}
	else {
		Gout = conductance();
		Iout = Gout*(Vm()-Vrev()); 
	}
}

// Empty caches for time step and current time
void SingleExpSynapticCond::emptyCaches()
{
	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	// Set state cached time to an impossible time
	// and class cache step size to an also impossible value.
	_cachedTime = InfiniteFuture;
	pCC->H = -1;
}

// Clear state value to zero
void SingleExpSynapticCond::clearState()
{
	_s1Now = 0;
	saveStartingState();
	emptyCaches();
}

// Save the current state as the starting state for a time step
void SingleExpSynapticCond::saveStartingState() 
{
	_s1Start = _s1Now;
}

// Restore the current internal state to time step start values
void SingleExpSynapticCond::restoreToStartingState() 
{
	_s1Now = _s1Start;
}

// Update cached values to reflect change in time step
void SingleExpSynapticCond::updateCacheForStepSize() 
{
	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	pCC->Exp1 = exp(- pCC->H/pCC->Tau1);
}

// Advance state in time by the amount of the current step size in _H
void SingleExpSynapticCond::advanceState()
{
	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	_s1Now *= pCC->Exp1;
}

// Update the current state value to reflect an AP Event
void SingleExpSynapticCond::updateForAPEvent(ActionPotentialEvent* apEvent)
{
	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	SimTime		hSyn;
	Number		wSyn;

	// If there was no release, then there is no effect
	if (apEvent->quantity()==0)
		return;

	// Compute time since arrival of the AP and associated exponentials.
	hSyn = currentTime() - apEvent->eventTime();

	// Get effective weight for this event
	wSyn = adjustedSynapticWeight(apEvent)*apEvent->quantity();

	// Get the contribution of one AP and add it to the others.
	// Use QDExp to speed up processing time per spike.
	_s1Now += wSyn*qdexp(-hSyn/pCC->Tau1);
}

// Estimate the mean value of s1 over an interval.
Number SingleExpSynapticCond::s1Mean(SimTime h)
{
	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	SimTime		tEnd=currentTime();
	SimTime		tStart = tEnd-h;
	SimTime		tc1 = pCC->Tau1;

	double		s1int;

	updateCachedValues();

	// Just in case h is degenerate, stop now with the current value
	if (h==0) {
		return _s1Now;
	}

	// Get the integral of s1 over the time interval [tStart tEnd)
	// ignoring any action potential events that occurred over this 
	// interval. The contribution of such events is adjusted below.
	s1int=tc1*(qdexp(h/tc1)-1)*_s1Now;

	// Peek at AP events on the queue to see if any
	// fall into the interval and make adjustments
	ActionPotentialEventQueue&	apQ = apQueue();
	ActionPotentialEvent*		apPtr;

	// Process any events on the queue up to the current time.
	for (apPtr = apQ.begin(); 
		apPtr != apQ.end() && apPtr->eventTime()<tEnd; 
		apPtr  = apQ.next(apPtr) ) {

		// If the AP event lies in the target interval,
		// back out the contribution up to the start of the event
		// that was already included in the above s1int computation.
		if (apPtr->eventTime()>=tStart) {

			SimTime	hb = apPtr->eventTime()-tStart;
			Number	wSyn = adjustedSynapticWeight(apPtr)*apPtr->quantity();

			s1int -= wSyn*tc1*(qdexp(hb/tc1)-1);
		}
	}

	// Return the mean value over the interval
	return s1int/h;
}

// Update the cached values as of the current time
void SingleExpSynapticCond::updateCachedValues()
{
	const SimTime		tNow = currentTime();

	// If the cache is already up to date, stop now
	if (fabs(_cachedTime-tNow)<=EpsilonTime)
		return;

	// Access the class cache
	SingleExpSynapticCondClassCache* pCC=pSingleExpClassCache();

	SimTime				h; // time step size so far

	// Recompute state values when needed.

	// If cached time is beyond where we are now, go
	// back to the state at the start of the time step
	// and move forward from there. Allow processing of
	// any late arriving events currently on the queue
	// by setting the cached time to the remote past.
	if (_cachedTime > tNow) {
		_cachedTime = InfinitePast;
		restoreToStartingState();
		h=tNow-timeStepStart();
	}
	else {
		// Set time step based on time since last good step
		h = tNow - _cachedTime;
	}

	// Recompute factors derived from step size, but only if changed
	// Allow a small amount of room for error in judging when the
	// time step has changed since the step size is calculated and
	// may be subject to roundoff error.
	if (fabs(h - pCC->H) > EpsilonTime) {
		pCC->H = h;
		updateCacheForStepSize();
	}

	// Advance the local state to the current time. This will
	// typically be overridden in subclasses as needed.
	if (h>EpsilonTime) {
		advanceState();
	}

	// Process any events since the old cache time up to the current time.
	updateForQueuedAPEvents(_cachedTime, tNow);

	// Our work here is done, update the cache time stamp.
	_cachedTime = tNow;
}


// Update the state at the end of a time step to be
// ready for the next time step.
void SingleExpSynapticCond::updateStateForTimeStep()
{
	// Make sure cache is up to date. For some ODE solvers
	// the previous evaluation time might not be the end of
	// the time step as a whole.
	updateCachedValues();

	// Save any subclass state values
	saveStartingState();
}

// Handle an action potential arriving after the time step
// in which it would be processed is already completed.
void SingleExpSynapticCond::handleLateAPEvent(ActionPotentialEvent* apEvent)
{
	// Take default action from the superclass
	SynapticConductance::handleLateAPEvent(apEvent);

	// Clear the cache time to force reprocessing of the AP queue
	_cachedTime = InfiniteFuture;
}

// Return the internal state value for reporting
Number SingleExpSynapticCond::internalStateValue(int k)
{
	switch (k) {
	case 0:
		return s1();

	default:
		FatalError("(SingleExpSynapticCond::internalStateValue) Invalid state index.");
		return 0;
	}
}



// ====================================================================
// DualExpSynapticCondWithVarTau class body
// ====================================================================



// Construct a new instance from a maximum conductance
DualExpSynapticCondWithVarTau::DualExpSynapticCondWithVarTau(Number gMaxValue) 
{
	// Save gMax value now so that it can be accessed
	// during setInitialProcessing which invokes gMax.
	_gMax = gMaxValue;
} 

// Destroy the current instance
DualExpSynapticCondWithVarTau::~DualExpSynapticCondWithVarTau() {}

// Clear the state including the state vector value
void DualExpSynapticCondWithVarTau::clearState()
{
	// Allow the superclass to clear its state
	SingleExpSynapticCond::clearState();

	// Set starting state values to zero
	stateValue(0) = 0;
}

// Set conductance based on maximum value attainable.
// Note that this requires tau1 and tau2 values via the subclass
// and is normally done only during the start-up sequence.
void DualExpSynapticCondWithVarTau::gMax(Number gm)
{
	Number q10 = Q10Factor();

	// Save the gMax value in case it is needed later,
	// for example, if gMax is changed before setInitialState.
	_gMax = gm;

	// Set the scaled synapse conductance so that the peak value
	// to an impulse response comes out to be gMax.
	SynapticResponse::gAbsolute( gm/peakDualExpResp( tau1()/q10, nominalTau2()/q10));

	// Clear any values based on gmax
	emptyCaches();
}

// Compute derivatives by evaluating the ODE
void DualExpSynapticCondWithVarTau::computeDerivatives()
{
	// Skip the update if there are no synapses.
	// This avoids computing tau2 which may be expensive.
	if ( isInactive() ) {
		derivValue(0) = 0;
		return;
	}

	// Get current values of time constants
	SimTime tc2 = tau2()/Q10Factor();
	SimTime tc1 = pSingleExpClassCache()->Tau1;

	// Apply the ODE to get the derivative
	derivValue(0) = -s2()/tc2 + s1()/tc1;
}

// Apply an implicit trapezoid rule to update the s2 state
void DualExpSynapticCondWithVarTau::localStateUpdate(SimTime h, CNStepType stepType)
{
	// Skip the update if there are no synapses.
	// This avoids computing tau2 which may be expensive.
	if ( isInactive() ) {
		stateValue(0) = 0;
		return;
	}

	// Get current values of time constants
	SimTime tc2 = tau2()/Q10Factor();
	SimTime tc1 = pSingleExpClassCache()->Tau1;

	// Apply the update rule. No special handling of step type is done.
	stateValue(0) += h/(1+h/tc2)*(-s2()/tc2+s1Mean(h)/tc1);
}



// ====================================================================
// DualExpSynapticCond class body
// ====================================================================



// Construct a new instance from a maximum conductance
DualExpSynapticCond::DualExpSynapticCond(Number gMaxValue) 
{
	// Save the gMax value now so that it can be accessed
	// during setInitialState which invokes gMax.
	// Note that at this stage of construction, the 
	// gMax(Number x) accessor cannot be used because
	// time constants are currently unknown.
	_gMax = gMaxValue;
} 

// Destroy the current instance
DualExpSynapticCond::~DualExpSynapticCond() {}

// Initialize the class cache by adding to supperclass.
void DualExpSynapticCond::initializeClassCache()
{
	// Access the class cache
	DualExpSynapticCondClassCache* pCC = pDualExpClassCache();

	// Stop now if already initialized
	if (pCC->isInitialized) return;

	// Cache parameters values. Doing this here ensures that
	// the cache is updated at least once while permitting some
	// flexibility up until the ODE solver starts running.
	Number q10=Q10Factor();

	pCC->Tau1 = tau1()/q10;
	pCC->Tau2 = tau2()/q10;

	pCC->useAlpha = fabs((pCC->Tau1-pCC->Tau2)/pCC->Tau1)<1e-4;
	pCC->H = -1; // i.e. no value yet

	// Indicate that the cache has now been initialized
	pCC->isInitialized = true;
}

// Set conductance based on maximum value attainable.
// Note that this requires tau1 and tau2 values via subclass.
void DualExpSynapticCond::gMax(Number gm)
{
	// Save the gMax value in case it is needed later,
	// for example, if gMax is changed before setInitialState.
	_gMax = gm;

	// Get the temperature adjustment factor
	Number q10 = Q10Factor();

	// Set the scaled synapse conductance so that the peak value
	// to an impulse response comes out to be gMax.
	SynapticResponse::gAbsolute( gm/peakDualExpResp(tau1()/q10,tau2()/q10) );

	// Force recomputation of any cached values
	emptyCaches();
}

// Advance state in time by the amount of the current step size in _H
void DualExpSynapticCond::advanceState()
{
	// Access the class cache for the current time step
	DualExpSynapticCondClassCache* pCC = pDualExpClassCache();

	// Update the state values by solving the ODE explicitly
	_s2Now = pCC->Exp2 * _s2Now + pCC->ImpulseResp * _s1Now;
	_s1Now = pCC->Exp1 * _s1Now;
}

// Clear state values to zero
void DualExpSynapticCond::clearState()
{
	_s1Now = 0;
	_s2Now = 0;
	saveStartingState();
	emptyCaches();
}

// Save the current state as the starting state for a time step
void DualExpSynapticCond::saveStartingState()
{
	_s1Start = _s1Now;
	_s2Start = _s2Now;
}

// Restore the current internal state to time step start values
void DualExpSynapticCond::restoreToStartingState()
{
	_s1Now = _s1Start;
	_s2Now = _s2Start;
}

// Update cached values to reflect change in time step
void DualExpSynapticCond::updateCacheForStepSize()
{
	// Access the class cache for the current time step
	DualExpSynapticCondClassCache* pCC = pDualExpClassCache();

	// Simplify notation
	SimTime		h = pCC->H;
	double		tc1 = pCC->Tau1;
	double		tc2 = pCC->Tau2;

	// Test for alpha function case (or not)
	if (pCC->useAlpha) {

		// Set values for the alpha function case
		pCC->Exp1 = pCC->Exp2 = exp(-h/tc1);
		pCC->ImpulseResp = h /tc1 * pCC->Exp1;
	}
	else {
		// Otherwise, need to do both exponent values
		pCC->Exp1 = exp(-h/tc1);
		pCC->Exp2 = exp(-h/tc2);
		pCC->ImpulseResp = tc2/(tc2-tc1)*(pCC->Exp2 - pCC->Exp1);
	}
}

// Update current conductance to reflect an AP event
void DualExpSynapticCond::updateForAPEvent(ActionPotentialEvent* apEvent)
{
	// If there was no release, then there is no effect
	if (apEvent->quantity()==0)
		return;

	// Access the class cache for the current time step
	DualExpSynapticCondClassCache* pCC = pDualExpClassCache();

	double		e1,e2;
	SimTime		hSyn;
	Number		wSyn;

	// Simplify notation
	double	tc1 = pCC->Tau1;
	double	tc2 = pCC->Tau2;

	// Compute time since arrival of AP and associated exponentials.
	hSyn = currentTime() - apEvent->eventTime();

	// Get the weight for this particular synapse and action potential.
	wSyn = adjustedSynapticWeight(apEvent)*apEvent->quantity();

	// Get the contribution of one AP and add it to the state.
	if (pCC->useAlpha) {
		// Use alpha function
		e1 = qdexp(-hSyn/tc1);
		_s2Now += wSyn * hSyn/tc1 * e1;
		_s1Now += wSyn * e1;
	}
	else {
		// Use dual exponent form.
//		e1 = qdexp(-hSyn/tc1);
//		e2 = qdexp(-hSyn/tc2);
		e1 = exp(-hSyn/tc1);
		e2 = exp(-hSyn/tc2);
		_s2Now += wSyn * tc2/(tc2-tc1) * (e2-e1);
		_s1Now += wSyn * e1;
	}
}

// Return the internal state value for reporting.
Number DualExpSynapticCond::internalStateValue(int k)
{
	switch (k) {
	case 0:
		return s1();

	case 1:
		return s2();

	default:
		FatalError("(DualExpSynapticCond::internalStateValue) Invalid state index.");
		return 0;
	}
}



// ====================================================================
// TripleExpSynapticCond class body
// ====================================================================



// Construct a new instance from a maximum conductance
TripleExpSynapticCond::TripleExpSynapticCond(Number gMaxValue) 
{
	// Save the gMax value now so that it can be accessed
	// during setInitialState which invokes gMax.
	// Note that at this stage of construction, the 
	// gMax(Number x) accessor cannot be used because
	// time constants are currently unknown.
	_gMax = gMaxValue;
} 

// Destroy the current instance
TripleExpSynapticCond::~TripleExpSynapticCond(){}

// Initialize the class cache]
void TripleExpSynapticCond::initializeClassCache()
{
	// Save the class cache address
	TripleExpSynapticCondClassCache* pCC = pTripleExpClassCache();

	// Stop now if already initialized
	if (pCC->isInitialized) return;

	// Cache parameters values. Doing this here ensures that
	// the cache is updated at least once while permitting some
	// flexibility up until the ODE solver starts running.
	Number q10=Q10Factor();

	pCC->C2 = c2();
	pCC->Tau1 = tau1()/q10;
	pCC->Tau2 = tau2()/q10;
	pCC->Tau3 = tau3()/q10;

	pCC->useAlpha2 = fabs((pCC->Tau1-pCC->Tau2)/pCC->Tau1)<1e-4;
	pCC->useAlpha3 = fabs((pCC->Tau1-pCC->Tau2)/pCC->Tau1)<1e-4;

	pCC->H = -1; // i.e. no value yet

	// Indicate that the cache has now been initialized
	pCC->isInitialized = true;
}

// Set conductance based on maximum value attainable
// Note that this uses tau values from the subclass.
void TripleExpSynapticCond::gMax(Number gm)
{
	// Save the gMax value in case it is needed later,
	// for example, if gMax is changed before setInitialState.
	_gMax = gm;

	// Get the temperature factor
	Number q10=Q10Factor();

	// Set the scaled synapse conductance so that the peak value
	// to an impulse response comes out to be gMax.
	SynapticResponse::gAbsolute( gm/peakTripleExpResp(tau1()/q10,tau2()/q10,tau3()/q10,c2() ) );

	// Force recomputation of any cached values	
	emptyCaches();
}

// Advance state in time by the amount of the current step size in _H
void TripleExpSynapticCond::advanceState()
{
	// Access the class cache
	TripleExpSynapticCondClassCache* pCC=pTripleExpClassCache();

	// Update the state values by solving the ODE explicitly
	_s2Now = pCC->Exp2 * _s2Now + pCC->ImpulseResp2 * _s1Now;
	_s3Now = pCC->Exp3 * _s3Now + pCC->ImpulseResp3 * _s1Now;
	_s1Now = pCC->Exp1 * _s1Now;
}

// Clear state values to zero
void TripleExpSynapticCond::clearState()
{
	_s1Now = 0;
	_s2Now = 0;
	_s3Now = 0;
	saveStartingState();
	emptyCaches();
}

// Save the current state as the starting state for a time step
void TripleExpSynapticCond::saveStartingState()
{
	_s1Start = _s1Now;
	_s2Start = _s2Now;
	_s3Start = _s3Now;
}

// Restore the current internal state to time step start values
void TripleExpSynapticCond::restoreToStartingState()
{
	_s1Now = _s1Start;
	_s2Now = _s2Start;
	_s3Now = _s3Start;
}

// Update cached values to reflect change in time step
void TripleExpSynapticCond::updateCacheForStepSize()
{
	// Access the class cache
	TripleExpSynapticCondClassCache* pCC=pTripleExpClassCache();

	// Simplify notation
	SimTime h	= pCC->H;
	SimTime tc1	= pCC->Tau1;
	SimTime tc2	= pCC->Tau2;
	SimTime tc3	= pCC->Tau3;

	// Set values for tau1
	pCC->Exp1 = exp(-h/tc1);

	// Set ImpulseResp for tau2 case
	if (pCC->useAlpha2) {
		pCC->Exp2 = pCC->Exp1;
		pCC->ImpulseResp2 = h/tc1 * pCC->Exp1;
	}
	else {
		pCC->Exp2 = exp(-h/tc2);
		pCC->ImpulseResp2 = tc2/(tc2-tc1)*(pCC->Exp2 - pCC->Exp1);
	}

	// Set ImpulseResp for tau3 case
	if (pCC->useAlpha3) {
		pCC->Exp3 = pCC->Exp1;
		pCC->ImpulseResp2 = h/tc1 * pCC->Exp1;
	}
	else {
		pCC->Exp3 = exp(-h/tc3);
		pCC->ImpulseResp3 = tc3/(tc3-tc1)*(pCC->Exp3 - pCC->Exp1);
	}
}

// Update current conductance to reflect an AP event
void TripleExpSynapticCond::updateForAPEvent(ActionPotentialEvent* apEvent)
{
	// If there was no release, then there is no effect
	if (apEvent->quantity()==0)
		return;

	// Access the class cache
	TripleExpSynapticCondClassCache* pCC=pTripleExpClassCache();

	double		wSyn,hSyn,e1,e2,e3;

	// Simplify notation
	SimTime tc1	= pCC->Tau1;
	SimTime tc2	= pCC->Tau2;
	SimTime tc3	= pCC->Tau3;

	// Compute time since arrival of AP and associated exponentials.
	hSyn = currentTime() - apEvent->eventTime();

	// Get the weight for this particular synapse and action potential.
	wSyn = adjustedSynapticWeight(apEvent);

	// Get the contribution of one AP and add it to the state.
	e1 = qdexp(-hSyn/tc1);

	if (pCC->useAlpha2) {
		_s2Now += wSyn * hSyn/tc1 * e1;
	}
	else {
		e2 = qdexp(-hSyn/tc2);
		_s2Now += wSyn * tc2/(tc2-tc1) * (e2-e1);
	}
	if (pCC->useAlpha3) {
		_s3Now += wSyn * hSyn/tc1 * e1;
	}
	else {
		e3 = qdexp(-hSyn/tc3);
		_s3Now += wSyn * tc3/(tc3-tc1) * (e3-e1);
	}
	_s1Now += wSyn * e1;
}

// Return the net conductance combining both s2 and s3 values.
Number TripleExpSynapticCond::conductance()
{
	// Bring the cache up to date
	updateCachedValues();

	// Access the class cache
	TripleExpSynapticCondClassCache* pCC=pTripleExpClassCache();

	// Return the weighted conductance
	return g()*(pCC->C2*_s2Now + (1-pCC->C2)*_s3Now); 
}

// Return the internal state value for reporting.
Number TripleExpSynapticCond::internalStateValue(int k)
{
	switch (k) {
	case 0:
		return s1();

	case 1:
		return s2();

	case 2:
		return s3();

	default:
		FatalError("(TripleExpSynapticCond::internalStateValue) Invalid state index.");
		return 0;
	}
}



// ====================================================================
// Synapse class body
// ====================================================================



// Construct a new instance. This is always done by the postsynaptic side,
// which is responsible for adding the synapse to its list after object creation.
Synapse::Synapse(
	AxonProcess*			axon,	// presynaptic axon for this synapse
	SynapticResponse*		resp,	// postsynaptic process
	Number					dist)	// axonal distance
{
	// Clear link pointers
	_nextPresynaptic = NULL;
	_nextPostsynaptic = NULL;

	// Set values from passed parameters
	_presynapticProcess = axon;
	_postsynapticProcess = resp;
	
	// Delay is the sum of propagation plus synaptic delays
	_delay = axon->propRate() * dist;

	// Link up with presynaptic axon process
	axon->addSynapse(this);
}

// Destroy this instance
Synapse::~Synapse() 
{
	// Remove from the associated axon list if present
	if (_presynapticProcess!=NULL) {
		_presynapticProcess->removeSynapse(this);
	}

	// Also destroy any associated state data if not
	// done already by the postysnaptic side.
	if (_postsynapticProcess!=NULL) {
		_postsynapticProcess->destroySynapse(this);
	}
}

// Add to presynaptic list before the synapse provided (which may be NULL)
void Synapse::addBeforeInPresynapticList(Synapse* synapseToFollow)
{

	_nextPresynaptic = synapseToFollow;
}

// Add to postsynaptic list before the synapse provided (which may be NULL)
void Synapse::addBeforeInPostsynapticList(Synapse* synapseToFollow)
{

	_nextPostsynaptic = synapseToFollow;
}

// Remove the next synapse after this one in the presynaptic list
void Synapse::removeNextFromPresynapticList()
{
	if (_nextPresynaptic!=NULL) {
		_nextPresynaptic = _nextPresynaptic->nextPresynaptic();
	}
}

// Remove the next synapse after this one in the postsynaptic list
void Synapse::removeNextFromPostsynapticList()
{
	if (_nextPostsynaptic!=NULL) {
		_nextPostsynaptic = _nextPostsynaptic->nextPostsynaptic();
	}
}

// Clear the presynaptic process in preparation for list deletion
Synapse* Synapse::clearPresynaptic()
{
	Synapse* next = _nextPresynaptic;

	// Clear the link
	_presynapticProcess = NULL;
	_nextPresynaptic = NULL;

	// See if the other link is also NULL, in which case the
	// buffer containing this synapse can be freed.
	// The object deconstructor is bypassed.
	if (_postsynapticProcess == NULL) {
		delete this;
	}

	return next;
}

// Clear the postsynaptic process in preparation for list deletion
Synapse* Synapse::clearPostsynaptic()
{
	Synapse* next = _nextPostsynaptic;

	// Clear the link
	_postsynapticProcess = NULL;
	_nextPostsynaptic = NULL;

	// See if the other link is also NULL, in which case the
	// buffer containing this synapse can be freed.
	// The object deconstructor is bypassed.
	if (_presynapticProcess == NULL) {
		delete this;
	}

	return next;
}

// Handle a new action potential by passing it along to the postsynaptic process.
// If there is no postsynapticProcess, this is a dead entry awaiting cleanup.
void Synapse::signalActionPotential(ActionPotentialEvent* apev)
{
	// Add in propagation delay and pass along the event
	if (postsynapticProcess()!=NULL) {
		postsynapticProcess()->signalActionPotential(apev->eventTime()+delay(),apev);
	}
}



// ====================================================================
// PlasticityRule class body
// ====================================================================



// Constructors and destructor
PlasticityRule::PlasticityRule() 
{
	_synapticCond = NULL;
	_stateOffset = 0;
	_isEnabled = true;
}

PlasticityRule::~PlasticityRule() 
{
	// Unhook from any owning synaptic conductance
	if (_synapticCond!=NULL) {
		_synapticCond->unhookFromRule(this);
	}
}

// Set the associated synaptic conductance
void PlasticityRule::synapticCond(SynapticConductance* sc)
{
	// Enforce limits on changing relationship
	if (_synapticCond!=NULL && _synapticCond!=sc && sc!=NULL) {
		FatalError("(PlasticityRule) Cannot change associated synaptic conductance");
	}

	// Set the synaptic conductance and inherit its model
	_synapticCond = sc;
	if (sc!=NULL) {
		model( sc->model() );
	}
}

// Take action at the end of the time step before AP purge.
// Default is to apply the rule one event at a time in order.
void PlasticityRule::applyEndOfStep(ActionPotentialEventQueue& apQueue)
{
	int						idx;
	ActionPotentialEvent*	apPtr;
	SimTime					now=currentTime();

	// Process current events one at a time in order
	// up to the end of the current step time.
	for (idx=0;idx<apQueue.size();idx++) {
		apPtr=apQueue.peek(idx);
		if (apPtr->eventTime()>=now) {
			break;
		}
		applyRule(apPtr);
	}
}

// Make the subclass provides a creation function if
// plasticity state size is not zero.
void PlasticityRule::createPlasticityState(Synapse* syn, Number wght)
{
	if (plasticityStateSize()>0 ) {
		FatalError("(PlasticityRule::createPlasticityState) subclass responsibility");
	}
}



// ====================================================================
// RandomReleaseRule class body
// ====================================================================



// Constructors and destructor
RandomReleaseRule::RandomReleaseRule(Number probRel) 
{
	_releaseProbability = probRel;
}

RandomReleaseRule::~RandomReleaseRule() {}

// Use the release probability to do random release
void RandomReleaseRule::finalizeAPEvent(ActionPotentialEvent* apEvent)
{
	Number pr = releaseProbability()<0 ?
		synapticCond()->releaseProbability() :
		releaseProbability();

	// Set event quantity to 1 or 0 randomly
	if (pr<1) {
		apEvent->quantity(synapticCond()->runif()<pr ? 1 : 0);
	}
	else {
		apEvent->quantity(1);
	}
	
	// Tag the event as final
	apEvent->isFinal(true);
}



// ====================================================================
// ActionPotentialEvent class body
// ====================================================================



// Construct a new event
ActionPotentialEvent::ActionPotentialEvent(
	SimTime t,					// Time of the action potential
	AxonProcess* ax,			// Axon where it originated
	Synapse* syn)				// Synapse receiving this event
:	Event(t)
{
	// Fill in event data
	_axon = ax;
	_synapse = syn;
	_quantity = 1;
	_isi = 0;
	_firingRate = 0;
	_isFinal = false;
}

// Destroy the current instance
ActionPotentialEvent::~ActionPotentialEvent() {}



// ====================================================================
// ActionPotentialEventQueue class body
// ====================================================================



// Constructors and destructor
ActionPotentialEventQueue::ActionPotentialEventQueue(int initialCapacity)
{
	// Initialize an empty queue with the capacity indicated
	_capacity = initialCapacity;
	if (_capacity==0) {
		_queue = NULL;
	}
	else {
		_queue = new ActionPotentialEvent[_capacity];
	}
	clear();
}

ActionPotentialEventQueue::~ActionPotentialEventQueue()
{
	// Note that ActionPotentialEvent destructor does nothing.
	// hence it is not invoked before deleting the queue storage.

	delete[] _queue;
}

// Set the queue to an empty state.
void ActionPotentialEventQueue::clear()
{
	// Since ActionPotentialEvent destructor does nothing,
	// it is not invoked. This would obviously need
	// to be changed if the destructor was otherwise.
	
	_size = 0;
	_begin = _queue;
	_end = _queue;
}

// Add a new event to the queue in time order
void ActionPotentialEventQueue::add(ActionPotentialEvent* apEvent)
{
	// Expand the queue if necessary to hold the new entry
	if (_size==_capacity) {

		ActionPotentialEvent*	newQueue;
		ActionPotentialEvent*	oldPtr;
		int						oldCap,newCap,k;

		// Allocate new expanded storage for the queue
		oldCap = capacity();
		newCap = oldCap==0 ? 2 : 2*oldCap;
		newQueue = new ActionPotentialEvent[newCap];

		// Copy the old contents
		oldPtr=_begin;
		for (k=0;k<_size;k++) {
			newQueue[k] = *oldPtr;
			oldPtr = privateNext(oldPtr);
		}

		// Install the new queue pointer
		delete[] _queue;
		_capacity = newCap;
		_queue = newQueue;
		_begin = &_queue[0];
		_end = &_queue[_size];
	}

	// If the queue is currently empty, quickly add an entry.
	// Reset the entry to the start of the queue array to improve
	// locality of reference later.
	if (_size == 0) {
		_queue[0] = *apEvent;
		_begin = &_queue[0];
		_end = &_queue[1];
		_size=1;
		return;
	}

	// Insert the event into the queue in order by time
	ActionPotentialEvent*	ptr2 = _end;
	ActionPotentialEvent*	ptr1 = privatePrev(_end);

	SimTime					t = apEvent->eventTime();
	int						k;

	// Locate the insertion point for the new event
	for (k=_size-1; k>=0 && ptr1->eventTime()>t;k--) {
		*ptr2 = *ptr1;
		ptr2 = ptr1;
		ptr1 = privatePrev(ptr1);
	}

	// Insert the new event and increment _end and size
	*ptr2 = *apEvent;
	_end = privateNext(_end); 
	_size++;
}


// Remove the first entry from the queue
void ActionPotentialEventQueue::removeFirst()
{
	_begin = privateNext(_begin);
	_size--;
}

// Access the next entry after the pointer provided
ActionPotentialEvent* ActionPotentialEventQueue::next(ActionPotentialEvent* ptr)
{
	ActionPotentialEvent* nextPtr;

	// If the queue is empty or unallocated, do not return a next entry
	if (_size==0) {
		return NULL;
	}
	// If the supplied pointer is NULL, return the first entry
	else if (ptr==NULL) {
		return _begin;
	}
	else {
		// Otherwise, get next entry in circular queue, if any
		nextPtr = privateNext(ptr);
		return nextPtr==_end ? NULL : nextPtr;
	}
}

// Access the previous entry before the pointer provided
ActionPotentialEvent* ActionPotentialEventQueue::previous(ActionPotentialEvent* ptr)
{
	// Check for special cases and handle them individually

	if (_size==0 || ptr==_begin) {
		return NULL;
	}
	else if (ptr==NULL) {
		return privatePrev(_end);
	}
	else {
		return privatePrev(ptr);
	}
}

// Peek at an entry given the offset from the first entry
ActionPotentialEvent* ActionPotentialEventQueue::peek(int n)
{
	int		k = n + (_begin-_queue);	// desired entry number except for wrap-around

	return &_queue[k<_capacity ? k : k-_capacity];
}
