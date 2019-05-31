// Provide classes for simulating a mouse moving in a maze
//
// Copyright 2007 John L Baker. All rights reserved.
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: placecell_baker_2003.cpp
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		6 March 2007
//
// Description:
//
// Place cells simulated here are purely phenomenological. They are simulated
// as PoissonNeurons with random firing patterns, which can be modulated based
// on gamma and theta rhythms. Place field firing probability is determined by
// by location within a maze.
//
// Interneurons are similarly phenomenological. Their firing is modulate by
// gamma and theta rhythms but is not location dependent.
//
// Place cells and interneurons are organized into groups (layers) based
// on similar properties or connectivity. Different connection policies
// can be applied to create synapses onto a target neuron.
//
// See the header file for references.

#include "placecell_baker_2003.h"
#include <iostream>

using namespace std;
using namespace BNSF;
using namespace BAKER_2003;



// --------------------------------------------------------------------
// PlaceCell class body
// --------------------------------------------------------------------



// Constructors and destructor
PlaceCell::PlaceCell(Model* m, PlaceCellLayer* pclayer)
: PoissonNeuron(m)
{
	using namespace UOM;

	// Save the associated place cell layer provided
	_layer = pclayer;

	// Set initial default values (examples only).
	// Note that the superclass has an initial
	// firing rate of 0, which is left as is.
	_meanFiringRate = 1*Hz;
	_distanceEstSD = 1*cm;
	_propDistanceEstSD = 0.2f;
	_thetaPhase = 0*radiansPerDegree;
	_precessionRate = 0*radiansPerDegree/meter;
	_isActive = true;
	_inactiveFiringRate = 0*Hz;

	// Set values indicating that no center has yet been set
	_centerX = numeric_limits<Number>::infinity();
	_centerY = numeric_limits<Number>::infinity();
	_Xvar = 0;
	_Yvar = 0;
	_hx = 0;
	_hy = 0;
	_densityMultiplier = 0;

}
PlaceCell::~PlaceCell() {}

// Access the maze subject
MazeSubject* PlaceCell::subject() 
{ 
	return layer()->subject(); 
}

// Access the maze associated with the subject
Maze* PlaceCell::maze() 
{ 
	return layer()->subject()->maze(); 
}

// Return the peak firing rate
Number PlaceCell::peakRate()
{
	return meanFiringRate()*_densityMultiplier;
}

// Set the place field center location and update fi
void PlaceCell::setCenter(Number x, Number y)
{
	// Save the values
	_centerX = x;
	_centerY = y;

	// Update field properties if possible
	updateProperties();
}

// Set the factors for computing distance estimate standard deviation
void PlaceCell::setDistanceEstSD(Number sdall, Number sdprop)
{
	_distanceEstSD = sdall;
	_propDistanceEstSD = sdprop;
	updateProperties();
}

// Set the flag indicating that this cell is active
// in the current environment.
void PlaceCell::isActive(bool b)
{
	_isActive = b;
}

// Return a location in local coordinates. See also the
// implementation of this in Maze. Subclasses may want to
// override this for a different place field scheme.
void PlaceCell::localCoordinates(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				northDist,		// Boundary dist along (0,1)
	Number&				southDist,		// Boundary dist along (0,-1)
	Number&				eastDist,		// Boundary dist along (1,0)
	Number&				westDist)		// Boundary dist along (-1,0)
{
	maze()->localCoordinates(locX, locY, northDist,southDist,eastDist,westDist);
}

// Provide a heading defined by local conditions. See also the
// implementation of this in Maze. Subclasses may want to
// override this for a different place field scheme.
void PlaceCell::localOrientation(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				hx,				// Unit vector x coord (output)
	Number&				hy)				// Unit vector y coord (output)
{
	maze()->localOrientation(locX, locY, hx, hy);
}

// Get the location of the nearest place field center. By default
// there is only one center, but subclasses may support more than one.
void PlaceCell::locateFieldCenter(
	Number					locX,			// Current location X coord
	Number					locY,			// Current location Y coord
	Number&					ctrX,			// Center location X coord
	Number&					ctrY)			// Center location X coord
{
	ctrX = centerX();
	ctrY = centerY();
}

// Check that everything is set before starting the simulation
void PlaceCell::simulationStarted()
{
	if (_densityMultiplier==0) {
		FatalError("(PlaceCell::simulationStarted) Initialization incomplete.");
		return;
	}

	// Let superclass know simulation is starting
	PoissonNeuron::simulationStarted();
}

// Update the firing rate based on current location
Number PlaceCell::updateFiringRate()
{
	MazeSubject*		sub = subject();
	Number				locx = subject()->locX();
	Number				locy = subject()->locY();

	Number				ctrx,ctry;
	Number				dx,dy;
	Number				dxr,dyr;
	Number				den,rate;

	// Get the nearest place field center for all subsequent processing.
	locateFieldCenter(locx,locy,ctrx,ctry);

	// Handle the case of an inactive cell.
	// Similarly, if the cell has a place field located
	// in the current inactive region, treat it as inactive.
	if ( isInactive() || (
			layer()->inactiveRegion()!=NULL &&
			layer()->inactiveRegion()->containsPoint(ctrx,ctry))) {
		firingRate( inactiveFiringRate() );
		return inactiveFiringRate();
	}

	// Otherwise proceed with setting the rate based
	// on the current location.

	// Get the difference between the current location and field center
	dx = locx - ctrx;
	dy = locy - ctry;

	// Rotate the difference vector based on the local orientation
	dxr =  dx*_hx + dy*_hy;
	dyr = -dx*_hy + dy*_hx;

	// Get the new firing rate from the spatial probability density
	den = _densityMultiplier*exp(-((dxr*dxr)/_Xvar+(dyr*dyr)/_Yvar)/2);
	firingRate( rate = meanFiringRate()*den );

	return rate;
}

// Update field properties if possible (otherwise ignore)
void PlaceCell::updateProperties()
{
	Number		xplus,xminus,yplus,yminus;
	Number		xsd,ysd;
	Number		ctrx,ctry;

	// Check for cases in which the update is not yet possible
	if (layer()==NULL || subject()==NULL || maze()==NULL) {
		return;
	}

	// Get the local orientation for the place field center
	// as well as local coordinates of the center in NSEW format.
	locateFieldCenter(subject()->locX(),subject()->locY(),ctrx,ctry);
	localOrientation(ctrx, ctry, _hx, _hy);
	localCoordinates(ctrx, ctry, yplus,yminus,xplus,xminus);

	// For each coordinate pick the smallest distance for the error model
	// and compute an associated standard deviation.
	xsd = distanceEstSD()+propDistanceEstSD()*minval(xplus,xminus);
	ysd = distanceEstSD()+propDistanceEstSD()*minval(yplus,yminus);

	// Save variances (avoids squaring std deviation each time)
	_Xvar = xsd*xsd;
	_Yvar = ysd*ysd;

	// Compute the density multiplier to scale the probability
	// density to unity assuming a uniform spatial occupancy.
	_densityMultiplier = 1/(2*Pi*xsd*ysd);
}

// Return the spike selection probability applying gamma and
// theta modulation as provided by the place cell layer.
Number PlaceCell::spikeSelectionProbability()
{
	// Check for a zero rate and handle directly if so
	if (firingRate()==0 ) {
		return 0;
	}

	// Get values from the subject
	const Number		hdx	= subject()->headingX();
	const Number		hdy	= subject()->headingY();
	const Number		locx = subject()->locX();
	const Number		locy = subject()->locY();

	// Get values from the layer
	const Number		ga	= layer()->gammaAmplitude();
	const Number		gf	= layer()->gammaFrequency();
	const Number		go	= layer()->gammaOffset();

	const Number		ta	= layer()->thetaAmplitude();
	const Number		tf	= layer()->thetaFrequency();
	const Number		to	= layer()->thetaOffset();

	// Working variables
	Number				ctrx,ctry;				// relevant place field center
	Number				dx,dy;					// distance from center
	Number				s;						// dist to center along heading vector
	Number				a;						// theta phase precession
	Number				gammaMod,thetaMod;		// oscillatory modulation values

	// Get theta phase precession, but only if active in the current context
	if (isActive() && precessionRate()!=0 ) {

		// Based on current position and heading determine whether the
		// the distance to the field center in the current heading direction.
		// This can be positive or negative depending on whether the direction
		// is towards or away from the field center.
		locateFieldCenter(locx,locy,ctrx,ctry);
		dx = ctrx - locx;
		dy = ctrx - locy;
		s = hdx*dx+hdy*dy;

		// Phase precession is a linear function of distance.
		// The amount of precession is arbitrarily limited to [-Pi,Pi].
		a = precessionRate() * s;
		if (a<-Pi)		a =-Pi;
		else if (a>Pi)	a = Pi;
	}
	else {
		a = 0;
	}

	// Get the gamma and theta modulations.
	gammaMod = plusval(1+ga*cos(2*Pi*gf*(_nextSpikeTime-go)));
	thetaMod = plusval(1+ta*cos(2*Pi*tf*(_nextSpikeTime-to)-(thetaPhase()+a)));
	
	// From all this, get a combined spiking probability.
	return isiAdjustedSelProb( gammaMod*thetaMod*firingRate() );
}



// --------------------------------------------------------------------
// Interneuron class body
// --------------------------------------------------------------------



// Constructors and destructor
Interneuron::Interneuron(Model* m, InterneuronLayer* inlayer) 
: PoissonNeuron(m)
{
	// Save the associated interneuron layer provided
	_layer = inlayer;

	// Initialize variables
	_thetaPhase = 0; // arbitary default
}

Interneuron::~Interneuron() {}

// Return the spike selection probability applying gamma-theta modulation
// as provided by the place cell layer
Number Interneuron::spikeSelectionProbability()
{
	// Get parameters from the layer
	const Number gf = layer()->gammaFrequency();
	const Number ga = layer()->gammaAmplitude();
	const Number go = layer()->gammaOffset();

	const Number tf = layer()->thetaFrequency();
	const Number ta = layer()->thetaAmplitude();
	const Number to = layer()->thetaOffset();

	const Number tph = thetaPhase();

	MazeSubject* sub = layer()->subject();

	// Working variables
	Number	oscMod,frate;

	// Compute the oscillatory modulation value as of the next spike time
	oscMod = plusval(1+ga*cos(2*Pi*gf*(_nextSpikeTime-go))) *
		plusval(1+ta*cos(2*Pi*tf*(_nextSpikeTime-to) - thetaPhase()));

	// Check to see if the subject is in a novelty region.
	// NOTE: performance would be improved in this case if the
	// layer were a subscriber of the subject and detected whether
	// the current location was in a novel region. For simple mazes
	// there may not be much improvement though.
	if (layer()->inactiveRegion() != NULL &&
		layer()->inactiveRegion()->containsPoint(sub->locX(),sub->locY()) ) {
		frate = layer()->inactiveRate();
	}
	else {
		frate = layer()->firingRate();
	}

	// Get a combined spiking probability for the next candidate spike time.
	// The overall firing rate comes from the layer rather than this object.
	return isiAdjustedSelProb(oscMod*frate);
}

// Get firing rate directly from the layer
Number Interneuron::firingRate()
{
	return layer()->firingRate();
}



// --------------------------------------------------------------------
// SawToothInterneuron class body
// --------------------------------------------------------------------



// Constructors and destructor
SawToothInterneuron::SawToothInterneuron(Model* m, InterneuronLayer* inlayer) 
: Interneuron(m, inlayer) 
{
	// Initialize variables
	_contraThetaPhase = 0; // default that must be changed before use.
}

SawToothInterneuron::~SawToothInterneuron() {}

// Return the spike selection probability applying gamma-theta modulation
// as provided by the place cell layer. Instead of a cosine wave, a
// saw-tooth wave is used with a maximum and minimum phase as specified.
Number SawToothInterneuron::spikeSelectionProbability()
{
	// Get parameters from the layer

	const Number gf = layer()->gammaFrequency();
	const Number ga = layer()->gammaAmplitude();
	const Number go = layer()->gammaOffset();

	const Number tf = layer()->thetaFrequency();
	const Number ta = layer()->thetaAmplitude();
	const Number to = layer()->thetaOffset();

	const Number twoPi = 2*Pi;

	MazeSubject* sub = layer()->subject();

	// Working variables
	Number	gammaMod,thetaMod,frate,omega,peakPhase,w;

	// Compute the oscillatory modulation value as of the next spike time
	// starting with the gamma modulation, which is sinusoidal.
	gammaMod = plusval(1+ga*cos(twoPi*gf*(_nextSpikeTime-go)));

	// Get the phase relative to peak phase and work out
	// a modulator by cases. The value derived (w) is in the
	// range [0 1] where w=0 at phase ctp and w=1 at tph.
	// 2w-1 then has range [-1 1] and replaces the cosine function
	// in the normal formulation of theta modulation.
	
	omega = fmod(Number(twoPi*tf*(_nextSpikeTime-to)-contraThetaPhase()),twoPi);
	peakPhase = fmod(Number(thetaPhase()-contraThetaPhase()),twoPi);

	// Fix the C-defined bug in fmod when handling negative numbers
	omega = omega>=0 ? omega : omega+twoPi;
	peakPhase = peakPhase>=0 ? peakPhase : peakPhase + twoPi;

	// The saw-tooth function is now a minimum at omega=0 and 2*pi.
	// Pick the rising or falling regions based on omega vs peak phase.
	// Try to handle the peakPhase==0 case correctly.
	w = omega<peakPhase ? omega/peakPhase : (twoPi-omega)/(twoPi-peakPhase);
	thetaMod = plusval(1+ta*(2*w-1));

	// Check to see if the subject is in a novelty region.
	// NOTE: performance would be improved in this case if the
	// layer were a subscriber of the subject and detected whether
	// the current location was in a novel region. For simple mazes
	// there may not be much improvement though.
	if (layer()->inactiveRegion() != NULL &&
		layer()->inactiveRegion()->containsPoint(sub->locX(),sub->locY()) ) {
		frate = layer()->inactiveRate();
	}
	else {
		frate = layer()->firingRate();
	}

	// Get a combined spiking probability for the next candidate spike time.
	// The overall firing rate comes from the layer rather than this object.
	return isiAdjustedSelProb(gammaMod*thetaMod*frate);
}



// --------------------------------------------------------------------
// CellLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
CellLayer::CellLayer(MazeSubject* sub, int numericId, UniformRandom* spkunif)
{
	using namespace UOM;
	
	// Save the values provided
	_subject = sub;
	_numericIdentifier = numericId;

	// Set initial values
	_gammaFrequency = 0;
	_gammaAmplitude = 0;
	_gammaOffset = 0;
	_thetaFrequency = 0;
	_thetaAmplitude = 0;
	_thetaOffset = 0;
	_inactiveRegion = NULL;

	// Create a default model and clock solver
	// with an arbitrary default time step.
	ClockSolver* mySolver = new ClockSolver;
	Model* myModel = new Model;
	model(myModel);
	mySolver->model(myModel);
	mySolver->timeStep(5*msec);

	// Set the spike randomizer in the model
	model()->uniformRandom(spkunif);
}

CellLayer::~CellLayer()
{
	// Tell any subscribers that this object is terminating
	changed(terminatedChange);

	// Remove any active change subscriptions
	if (_subject != NULL) {
		_subject->removeSubscriber(this);
	}

	// Delete the model and solver if not
	// already done by the subclass.
	if (model()!=NULL) {
		delete model()->solver();
		delete model();
	}
}

// Add a probe to all cells
void CellLayer::addProbeToAll(Probe* pr)
{
	int k;

	for (k=1;k<=numCells();k++) {
		cell(k)->addProbe(pr);
	}
}

// Remove a probe from all cells
void CellLayer::removeProbeFromAll(Probe* pr)
{
	int k;

	for (k=1;k<=numCells();k++) {
		cell(k)->removeProbe(pr);
	}
}


// Assign a numeric identifier to the layer and associated place cells.
void CellLayer::numericIdentifier(int n)
{
	int		k;
	int		nextId = n;
	int		ncell = numCells();

	// Save the identifier value
	_numericIdentifier = nextId++;

	// Assign identifier to place cells
	for (k=1;k<=ncell;k++) {
		cell(k)->numericIdentifier(nextId++);
	}
}


// --------------------------------------------------------------------
// PlaceCellLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
PlaceCellLayer::PlaceCellLayer(
	MazeSubject* sub, 
	int numCells, 
	int numericId, 
	UniformRandom* spkunif)
:	CellLayer(sub, numericId, spkunif)
{
	using namespace UOM;

	// Set initial values
	_numPlaceCells = numCells;
	_placeCells = NULL;
	_totalFiringRate = 0;
	_fractionActive = 1;

	// Subscribe to changes from the subject
	subject()->addSubscriber(this);
}

PlaceCellLayer::~PlaceCellLayer()
{
	int			k;

	// Remove any dependencies
	if (_subject!=NULL) {
		_subject->removeSubscriber(this);
	}

	// Delete the model here to speed up processing
	delete model()->solver();
	delete model();
	_model = NULL;

	// Delete allocated place cells
	if (_placeCells!=NULL) {
		for (k=1;k<=numCells();k++) {
			delete placeCell(k);
		}
	}
	delete[] _placeCells;
}

// Access an place cell by number (n=1 is the first)
PlaceCell* PlaceCellLayer::placeCell(int n)
{
	// Make sure n is valid
	if (n<1 || n>numCells() ) {
		FatalError("(PlaceCellLayer::placeCell) Invalid cell number requested.");
		return NULL;
	}

	return _placeCells[n-1];
}

// Receive an update notification from the subject
void PlaceCellLayer::updateFrom(ModelComponent* sub, int reason)
{
	int k;

	// If the subject is terminating, do not reference it further
	if (reason==terminatedChange) {
		_subject = NULL;
	}

	// If this is a state change, update all firing frequencies
	else if (reason==stateChange) {

		// Let each place cell update its firing rate
		// and get the new total.
		_totalFiringRate = 0;
		for (k=1;k<=numCells();k++) {
			_totalFiringRate += 
				placeCell(k)->updateFiringRate();
		}

		// Tell any interested parties about the change in rates
		changed(stateChange);
	}

	// Check for a maze being changed
	else if (reason==parameterChange) {
		if (maze()!=NULL) {

			// Check for a change of maze or the first setting
			// of the maze value. If so, update place cells
			// using the new maze coordinates.

			if (subject()->maze()!=subject()->previousMaze() ) {

				// Reset place cell locations based on the new maze.
				setPlaceFieldCenters();
			}
		}
	}

	// Any other change types are ignored -- hence no-op
	else {}
}

// Find place field centers with the maze. The results are returned
// in the preallocated arrays X and Y, which are assumed to be large enough.
void PlaceCellLayer::setPlaceFieldCenters()
{
	int			k;
	Number		x,y;

	// Set the new centers for each cell
	for (k=1;k<=numCells();k++) {

		// Have the maze select a random point using the 
		// random number stream reserved for this purpose.
		maze()->randomPointInMaze(x,y,subject()->locationRandomizer() );

		// Set center to the point found
		placeCell(k)->setCenter(x,y);
	}

	// Set a new active set (if needed)
	if (fractionActive()>0 && fractionActive()<1) {
		setFractionActive( fractionActive() );
	}
}

// Set the activity status of range of cells
void PlaceCellLayer::setFractionActive(Number activeFraction)
{
	int				k,n;
	bool			defaultActive;
	Number			prob;
	UniformRandom*	unif=subject()->networkRandomizer();

	// Save the fraction active
	_fractionActive = activeFraction;

	// Efficiently handle the special case when
	// active fraction is 1 or 0.
	if (activeFraction==0 || activeFraction==1.0) {
		for (k=1;k<=numCells();k++) {
			placeCell(k)->isActive( activeFraction!=0 );
		}
		return;
	}

	// Decide whether to select active or inactive cells
	// depending on which involves less random selection.
	if (activeFraction>0.5) {
		defaultActive = true;
		prob = 1 - activeFraction;
	}
	else {
		defaultActive = false;
		prob = activeFraction;
	}

	// Start by setting everying to the default status
	for (k=1;k<=numCells();k++) {
		placeCell(k)->isActive(defaultActive);
	}

	// Now pick cells at random to set to the opposite status.
	n=int( activeFraction*numCells()+0.5 );
	while(n>0) {
		k=int( unif->next()*numCells()+1 );
		if (placeCell(k)->isActive() == defaultActive) {
			placeCell(k)->isActive( !defaultActive );
			n--;
		}
	}
}

// Set a group of cells active or inactive
void PlaceCellLayer::setActivity(bool activeState, int from, int to)
{
	int k;

	for (k=from;k<=to;k++) {
		placeCell(k)->isActive(activeState);
	}
}

// Disable theta phase precession in all cells
void PlaceCellLayer::disableThetaPhasePrecession()
{
	int k;

	for (k=1;k<=numCells();k++) {
		placeCell(k)->precessionRate(0);
	}
}


// Print field centers to a file
void PlaceCellLayer::printPlaceFieldCenters(char* pathName)
{
	using namespace UOM;

	FILE*			out;
	int				i;
	PlaceCell*		cell;

	// Open the output file
	if (pathName!=NULL) {
		out = fopen(pathName,"w");
		if (out==NULL) {
			FatalError("(PlaceCellLayer::printFieldCenters) "
				"File open failed");
		}
	}
	else {
		out = stdout;
	}

	// Write a header line
	fprintf(out,"id,x,y,peakrate,isactive\n");

	// Write out the field centers in CSV format
	for (i=1;i<=numCells();i++) {
		cell = placeCell(i);
		fprintf(out,"%d,%g,%g,%g,%d\n",
			cell->numericIdentifier(),
			cell->centerX()/cm,
			cell->centerY()/cm,
			cell->peakRate()/Hz,
			cell->isActive() ? 1 : 0 );
	}

	// Close the output file
	if (out!=stdout) {
		fclose(out);
	}
}

// Allocate place cells.
void PlaceCellLayer::allocatePlaceCells()
{
	int k;

	// Make sure everything is ready to proceed
	if (_placeCells != NULL) {
		FatalError("(PlaceCellLayer::allocatePlaceCells) "
			"Trying to allocate place cells twice.");
		return;
	}

	// Allocate an array of place cells and initialize it
	_placeCells = new PlaceCell* [numCells()];
	for (k=1;k<=numCells();k++) {
		_placeCells[k-1] = newPlaceCell(k);
	}
	
	// Locate place field centers at random locations
	// excluding the currently unexplored novel region.
	if (maze()!=NULL) {
		setPlaceFieldCenters();
	}

	// Update all numeric identifiers
	numericIdentifier( numericIdentifier() );
}

// Allocate a new place cell with a center at the indicated location. 
// This can be extended by subclasses to reflect layer-specific properties.
PlaceCell* PlaceCellLayer::newPlaceCell(int k)
{
	return new PlaceCell(model(), this);
}



// --------------------------------------------------------------------
// InteneuronLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
InterneuronLayer::InterneuronLayer(
	MazeSubject* sub, 
	int numCells, 
	int numericId,
	UniformRandom* spkunif)
:	CellLayer(sub, numericId, spkunif)
{
	// Set initial values
	_numInterneurons = numCells;
	_interneurons = NULL;
	_firingRate = 0;
	_inactiveRate = 0;
}

InterneuronLayer::~InterneuronLayer()
{
	int			k;

	// Delete the model here to speed up processing
	delete model()->solver();
	delete model();
	_model = NULL;

	// Delete allocated interneurons
	if (_interneurons!=NULL) {
		for (k=1;k<=numCells();k++) {
			delete interneuron(k);
		}
	}
	delete[] _interneurons;
}

// Allocate interneurons.
void InterneuronLayer::allocateInterneurons()
{
	int k;

	// Make sure everything is ready to proceed
	if (_interneurons != NULL) {
		FatalError("(InterneuronLayer::allocateInterneurons) "
			"Trying to allocate interneurons twice.");
		return;
	}

	// Allocate and initialize interneurons
	_interneurons = new Interneuron* [numCells()];

	// Build the array of place cells
	for (k=1;k<=numCells();k++) {
		_interneurons[k-1] = newInterneuron(k);
	}

	// Set all numeric identifiers
	numericIdentifier( numericIdentifier() );

	// Set all firing rates
	firingRate( firingRate() );
}

// Allocate a new interneuron.. This can be extended by 
// subclasses to reflect layer-specific properties.
Interneuron* InterneuronLayer::newInterneuron(int k)
{
	Interneuron* cell = new Interneuron(model(), this);

	return cell;
}

// Access an interneuron by number (n=1 is the first)
Interneuron* InterneuronLayer::interneuron(int n)
{
	// Make sure n is valid
	if (n<1 || n>numCells() ) {
		FatalError("(InterneuronLayer::interneuron) Invalid cell number requested.");
		return NULL;
	}

	return _interneurons[n-1];
}



// --------------------------------------------------------------------
// ConnectionPolicy class body
// --------------------------------------------------------------------



// Constructors and destructor
ConnectionPolicy::ConnectionPolicy(
	CellLayer*		source, 
	UniformRandom*	randomizer)
{
	using namespace UOM;

	// Clear current target neuron pointer
	_target = NULL;

	// Save values provided
	_layer = source;
	_uniformRandom = randomizer;

	// Set starting values for params.
	// These will be updated as needed in subclasses.
	_minLaminarDist = 0;
	_maxLaminarDist = 0;
	_minBasalDist = 0;
	_maxBasalDist = 0;
	_enpassantDist = 0;
	_minAxonDist = 100*micron;
	_maxAxonDist = 100*micron;
	_dendriteSynDensity = 0;
	_somaSynDensity = 0;
	_ISSynDensity = 0;
	_synapseWeight = 1;
	_synapseWeightAlt = 1;
	_synapseType = NullTokenId;
	_synapseTypeAlt = NullTokenId;
	_connectOrder = sequentialOrder;

	// Zero statistics counters
	_totalSomaSynapses = 0;
	_totalISSynapses = 0;
	_totalAxonSynapses = 0;
	_totalApicalSynapses = 0;
	_totalBasalSynapses = 0;

	// Initialize connection index
	_nextToConnect = 1;
}
ConnectionPolicy::~ConnectionPolicy() {}

// Return total connected synapses
int ConnectionPolicy::totalSynapses()
{
	return 
		_totalSomaSynapses +
		_totalISSynapses +
		_totalAxonSynapses +
		_totalApicalSynapses +
		_totalBasalSynapses;
}

// Print connection stats for debug
void ConnectionPolicy::printSynapseCounts()
{
	cerr<<"totalSomaSynapses = "		<<totalSomaSynapses()<<endl
		<<"totalISSynapses = "			<<totalISSynapses()<<endl
		<<"totalAxonSynapses = "		<<totalAxonSynapses()<<endl
		<<"totalDendriteSynapses = "	
			<<totalApicalSynapses()<<" Apical, "
			<<totalBasalSynapses()<<" Basal, "
			<<totalApicalSynapses()+totalBasalSynapses()<<" Total"<<endl;
}

// Connect the target neuron with the cell layer
void ConnectionPolicy::connectWith(
	MorphologicalNeuron* target, 
	ConnectionSequence connSeq)
{
	// Save target and order for inner functions
	_target = target;
	_connectOrder = connSeq;

	// Connect with various parts of the cell
	connectWithSoma();
	connectWithIS();
	connectWithAxon();
	connectWithDendrites();

	// Clear the target pointer (just to be clean)
	_target = NULL;
}

// Connect with the soma
void ConnectionPolicy::connectWithSoma()
{
	int k,n;

	for (k=1;k<=_target->numSomaComp();k++) {

		Compartment* soma = _target->somaComp(k);
		n= numberToConnect(soma,1,somaSynDensity() );
		createSynapses(soma,n);
		_totalSomaSynapses += n;
	}
}

// Connect with the initial segment
void ConnectionPolicy::connectWithIS()
{
	int k,n;

	for (k=1;k<=_target->numISComp();k++) {

		Compartment* IS = _target->ISComp(k);
		n= numberToConnect(IS,1,ISSynDensity() );
		createSynapses(IS,n);
		_totalISSynapses += n;
	}
}

// Connect with the dendrites
void ConnectionPolicy::connectWithDendrites()
{
	Compartment* pcomp;
	int k;

	// Start sequential connections with the first cell
	_nextToConnect = 1;

	// Make connections for each compartment in turn
	for (k=1;k<=_target->numDendriteComp(); k++) {

		Number	f;
		int		nsyn;
		
		// See how much of a compartment is in range.
		// Optimize performance when out of range.
		f = dendriteFractionInRange(k, minLaminarDist(), maxLaminarDist() );
		f += dendriteFractionInRange(k, minBasalDist(), maxBasalDist() );
		if (f>0) {

			// Create and connect synapses
			pcomp = _target->dendriteComp(k);
			nsyn  = numberToConnect(pcomp, f, dendriteSynDensity() );
			createSynapses(pcomp, nsyn);

			// Count the synapses
			if (_target->isApicalDendrite(k) ) {
				_totalApicalSynapses += nsyn;
			}
			else {
				_totalBasalSynapses += nsyn;
			}
		}
	}
}

// Get an estimated fraction of a target dendrite compartment
// within a pathway defined by a range of distances from the origin.
Number ConnectionPolicy::dendriteFractionInRange(
	int					dendNbr,			// number of target dendrite
	Number				minDist,			// minimum distance
	Number				maxDist)			// maximum distance		
{
	using namespace UOM;

	// If both min and max are the same then there is nothing to do.
	// This tyically arises when both are defaults (0).
	if (minDist==maxDist)
		return 0;

	// Get values from the target neuron
	Compartment*		pDend = _target->dendriteComp(dendNbr);
	MorphologyEntry*	pMorph = _target->dendriteMorphologyEntry(dendNbr);
	Number				orX = _target->orientationX();
	Number				orY = _target->orientationY();
	Number				orZ = _target->orientationZ();

	// Working variables
	Number				s1,s2;				// distance orthogonal to layer
	Number				t1,t2;				// distances trimmed to boundaries
	Number				xprox,yprox,zprox;	// proximal end dendrite location
	Number				xdist,ydist,zdist;	// distal end dendrite location
	Number				temp;				// temporary register

	// Ensure that distance specs are right way around just
	// in case someone missed the lecture on negative numbers.
	if (minDist>maxDist) {
		FatalError("(ConnectionPolicy::dendriteFractionInRange) "
			"Minimum laminar distance greater than maximum");
		return 0;
	}

	// Locate the proximal and distals ends of the compartment 
	// assuming the dendrite is a perfect cylinder. Note that
	// the dx,dy,dz values of the compartment provide the
	// location of end points but do not imply the length of
	// the cylinder (dendrites are not really that straight).
	// Note that the morphology table is not converted for
	// units of measure and entries must be converted to microns.
	xprox = pMorph->x*micron - pMorph->dx*micron/2;
	yprox = pMorph->y*micron - pMorph->dy*micron/2;
	zprox = pMorph->z*micron - pMorph->dz*micron/2;
	xdist = pMorph->x*micron + pMorph->dx*micron/2;
	ydist = pMorph->y*micron + pMorph->dy*micron/2;
	zdist = pMorph->z*micron + pMorph->dz*micron/2;

	// Using the orientation of the layer, get a distance
	// aligned with the layer assuming that the layer is
	// defined by parallel planes (another rough assumption).
	s1 = orX*xprox + orY*yprox + orZ*zprox;
	s2 = orX*xdist + orY*ydist + orZ*zdist;

	// Get the fraction based on the distance along the axis
	// of the compartment that lies within the layer. Separate cases
	// are the easiest way to handle this logic.

	if (s1<minDist && s2<minDist) {
		// Dendrite is entirely below the layer
		return 0;
	}
	else if (s1>maxDist && s2>maxDist) {
		// Dendrite is entirely above the layer
		return 0;
	}
	else {
		// Dendrite is at least partially within the layer
		if (s2<s1) {
			temp = s2;
			s2 = s1;
			s1 = temp;
		}
		t1=maxval(s1,minDist);
		t2=minval(s2,maxDist);

		// It is possible that the dendrite is exactly
		// parallel to the layer boundary. If so, handle
		// as a special case to avoid divide by 0. Otherwise
		// the membrane area is prorated based on the
		// fraction of the dendrite (as a straight line)
		// that is within the layer boundaries.
		if (s1==s2) {
			return 1;
		}
		else {
			return (t2-t1)/(s2-s1);
		}
	}
}

// Find the number of synapses to connect. The number is derived
// from a Poisson distribution with a mean of the value determined 
// by applying synapse density to the implied dendrite area.
int ConnectionPolicy::numberToConnect(
	Compartment*		pcomp,				// target compartment
	Number				fraction,			// fraction of comp in layer
	Number				density)			// synapse density
{
	Number				numExpected;
	int					numActual;

	numExpected = fraction*density*pcomp->extraShellArea( enpassantDist() );

	// Assuming connections occur uniformly, use random round to convert
	// to an integer number of connections. This avoids the possibility
	// of consistently rounding to zero or one for small compartments.
	numActual = int( numExpected+uniformRandom()->next() );

	return numActual;
}

// Create the necessary synapses by selecting at random from 
// source cells in the place cell layer.
void ConnectionPolicy::createSynapses(
	Compartment*		pcomp,				// comparmtent to connect with
	int					numSynapses)		// number of synapses
{
	SynapticResponse*	resp1=NULL;
	SynapticResponse*	resp2=NULL;
	Number				len;

	int k,n;

	// If no connections are to be made stop here.
	// The required synaptic response may not even
	// exist in this compartment.
	if (numSynapses==0)
		return;

	// If connections are requested, make sure that cells 
	// are defined for the source layer.
	if (numSynapses>0 && _layer->numCells()==0) {
		FatalError("(ConnectionPolicy::createSynapses) No cells defined for layer.");
	}

	// Find the synaptic response(s) to be connected with
	if (synapseType()!=NullTokenId) {
		resp1 = pcomp->findSynapticResponse(synapseType() );
	}
	if (synapseTypeAlt()!=NullTokenId) {
		resp2 = pcomp->findSynapticResponse(synapseTypeAlt() );
	}

	// Create synapses
	for (k=1;k<=numSynapses;k++) {

		// Use the appropriate method for picking a source neuron
		switch (_connectOrder) {

		case sequentialOrder:
			if (_nextToConnect >= _layer->numCells() ) {
				_nextToConnect = 1;
			}
			n = _nextToConnect++;
			break;

		case randomOrder:
			n = int( _layer->numCells()*uniformRandom()->next()+1 );
			break;

		default:
			FatalError("(ConnectionPolicy::createSynapses) Invalid connect order");
		}

		// Connect with the synaptic response(s)
		len = axonLength(pcomp);	// get a common length
		if (resp1!=NULL) {
			resp1->createSynapse(_layer->cell(n)->axonProcess(),synapseWeight(),len);
		}
		if (resp2!=NULL) {
			resp2->createSynapse(_layer->cell(n)->axonProcess(),synapseWeightAlt(),len);
		}
	}
}

// Select an axon length for the current connection.
// By default this is a uniformly distributed random number
// between axon min and max distances.
Number ConnectionPolicy::axonLength(Compartment* pcomp)
{
	// If min and max are the same, skip the random number
	// generation overhead. Otherwise generate a distance.
	return minAxonDist()==maxAxonDist() 
		? minAxonDist() 
		: uniformRandom()->next(minAxonDist(), maxAxonDist());
}
