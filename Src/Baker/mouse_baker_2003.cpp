// Provide classes for simulating a mouse moving in a maze
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: mouse_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// To simulate the effects of NMDA knockout in the CA3 portion of the hippocampus,
// it is necessary to generate inputs similar in form to that which might be present
// for a mouse moving in a maze, typically a Morris water maze or its dry equivalent.
//
// The maze subject is nominally a mouse with motion parameters similar to those
// reported by Nakazawa. Only random swimming (or walking) is provided. Goal
// seeking behavior is not included in the simulation,
//
// This mouse model is basically the way a hippocampus specialist would think
// of a mouse, i.e. a hippocampus with some trivial behavior attached. No real
// attempt at even low fidelity simulation of mouse behavior is intended.


#include "mouse_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace BAKER_2003;



// --------------------------------------------------------------------
// Mouse class body
// --------------------------------------------------------------------



// Default random number seed values. These values are
// from a random units table in CRC Standard Mathematical 
// Tables 13th edition (1964 - a vintage year), which uses 
// values from a table prepared by the Bureau of Transport
// Economics and Statistics of the Interstate Commerce 
// Commission, Wash. DC. Using values from this source at 
// least avoids the issue of picking good and bad seed values
// in some preferential fashion. Doubles are used to simplify 
// combining pairs of 5 digit values into a single seed and 
// avoiding the C convention for integer constants in which 
// 010=8 because the constant is taken as octal (really).
// There is some possibility that these seeds could result in
// overlaps in the random streams, but that seems unlikely.

double Mouse::_defaultSeeds[8][6] = {

//	network       PC locations  movements     PC spikes     IN spikes     synapses
	0.1048015011, 0.0153602011, 0.8164791646, 0.6917914194, 0.6259036207, 0.2096999570, 
	0.9129190700, 0.2236846573, 0.2559585393, 0.3099589198, 0.2798253402, 0.9396534095, 
	0.5266619174, 0.3961599505, 0.2413048360, 0.2252797265, 0.7639364809, 0.1517924830, 
	0.4934032081, 0.3068019655, 0.6334858629, 0.4216793093, 0.0624361680, 0.0785616376, 
	0.3944053537, 0.7134157004, 0.0084974917, 0.9775816379, 0.3757039975, 0.8183716656,
	0.0612191782, 0.6046881305, 0.4968460672, 0.1411006927, 0.0126354613, 0.7792106907, 
	0.1100842751, 0.2775653498, 0.1860270659, 0.9065515053, 0.2191681825, 0.4439442880,
	0.9956272905, 0.5642069994, 0.9887231016, 0.7119418738, 0.4401348840, 0.6321321069};

// Constructors and destructor
Mouse::Mouse(
	Maze*				m,
	MovementType		moveType,
	double*				seedValues,
	bool				buildNet, 
	bool				doInit) 
: MazeSubject(m)
{
	using namespace UOM;

	// Do basic initialization of variables
	_ECPlaceCells = NULL;
	_DGPlaceCells = NULL;
	_CA3PlaceCells = NULL;
	_CA3AxoAxonicCells = NULL;
	_CA3BasketCells = NULL;
	_CA3BistratifiedCells = NULL;
	_CA3OLMCells = NULL;

	_PP = NULL;
	_MF = NULL;
	_AC = NULL;
	_INaxon = NULL;
	_INsoma = NULL;
	_INproximal = NULL;
	_INdistal = NULL;

	_targetCell = NULL;

	_buildNetwork = buildNet;
	_ownsTargetCell = false;

	// Set a default speed appropriate for a mouse
	_speed = 10*cm/sec;

	// Save the initial movement policy
	switch (moveType) {

	case wallFollower:
		// Set movement to follow walls such that the
		// target place field center is in the middle
		// of a 5-cm spatial bin.
		movementPolicy(new WallFollower(this, 7.5*cm) );
		break;

	case randomWalk:
		// Set random walk to more or less cover a 50-cm
		// diameter circle in a few minutes.
		movementPolicy(new RandomWalk(this, 4/sec, 2*sec));
		break;

	default:
		FatalError("(Mouse::Mouse) Invalid movement type");
	}

	// Set seeds for random number generator
	setSeedValues(seedValues==NULL ? _defaultSeeds[0] : seedValues);

	// Finish initialization if permitted
	if (doInit) {
		initialize();
	}
}

Mouse::~Mouse()
{
	// Delete any allocated storage
	delete _ECPlaceCells;
	delete _DGPlaceCells;
	delete _CA3PlaceCells;
	delete _CA3AxoAxonicCells;
	delete _CA3BasketCells;
	delete _CA3BistratifiedCells;
	delete _CA3OLMCells;

	delete _PP;
	delete _MF;
	delete _AC;
	delete _INaxon;
	delete _INsoma;
	delete _INproximal;
	delete _INdistal;

	// Since the movement policy is created here it is deleted.
	// Also, must set to NULL to prevent superclass from using ptr.
	delete _movementPolicy;
	_movementPolicy = NULL;

	// Delete the target cell if not supplied
	// externally.
	if (_ownsTargetCell) {
		delete _targetCell;
	}
}

// Set the maze for this subject and set location to its origin
void Mouse::maze(Maze* m)
{

	// Let superclass do the assignment
	MazeSubject::maze(m);

	// If the maze has been set, initialize location
	if (m!=NULL) {
		locX( m->originX() );
		locY( m->originY() );
	}
}

// Set the novelty region here and in affected layers
void Mouse::noveltyRegion(SpatialRegion* novelArea)
{
	_noveltyRegion = novelArea;

	// For now, assume EC is not affected by novelty.
	// Otherwise, pass novelty region to all layers.
	
	DGPlaceCells()->inactiveRegion(novelArea);
	CA3PlaceCells()->inactiveRegion(novelArea);
	
	CA3AxoAxonicCells()->inactiveRegion(novelArea);
	CA3BasketCells()->inactiveRegion(novelArea);
	CA3BistratifiedCells()->inactiveRegion(novelArea);
	CA3OLMCells()->inactiveRegion(novelArea);
}

// Set the target cell to a non-default value
void Mouse::targetCell(CA3PyramidalCell* cell)
{
	// Make sure target cell not changed once network is in place
	if (CA3PlaceCells() != NULL ) {
		FatalError("(Mouse::targetCell) Target cell cannot be changed after initialization.");
	}
	
	// Delete earlier cell if created here
	if (_ownsTargetCell) {
		delete _targetCell;
	}

	// Save the new cell and set ownership flag
	_targetCell = cell;				// save cell provided
	_ownsTargetCell = cell==NULL;	// does not own cell if supplied here
}

// Initialize the instance
void Mouse::initialize()
{
	// Set the starting point as the origin of the maze, if available
	if (maze() != NULL) {
		locX( maze()->originX() );
		locY( maze()->originY() );
	}

	// Start with a random heading
	randomizeHeading();

	// Create the necessary components
	if (_targetCell==NULL) {
		allocateTargetCell();
	}
	allocateCellLayers();
	setLayerActivities();
	if ( buildNetwork() ) {
		connectNetwork();
	}
	
	// Inform cell layers about any novelty region
	noveltyRegion( noveltyRegion() );
}

// Set randomizer seed values (array of 6 doubles)
void Mouse::setSeedValues(double* seedValues)
{
	networkRandomizer() ->setSeed(seedValues[0]);
	locationRandomizer()->setSeed(seedValues[1]);
	movementRandomizer()->setSeed(seedValues[2]);
	pcSpikeRandomizer() ->setSeed(seedValues[3]);
	inSpikeRandomizer() ->setSeed(seedValues[4]);
	synapticRandomizer()->setSeed(seedValues[5]);
}

// Allocate the target cell using a default cell definition
void Mouse::allocateTargetCell()
{
	_targetCell = new L56aPyramidalCell;
	_ownsTargetCell = true;
	_targetCell->model()->uniformRandom(synapticRandomizer());
	_targetCell->numericIdentifier(0);
}

// Allocate cell layers (EC, DG, CA3)
void Mouse::allocateCellLayers()
{
	using namespace UOM;

	UniformRandom*	pcspk = pcSpikeRandomizer();		// for PC Poisson cells
	UniformRandom*	inspk = inSpikeRandomizer();		// for IN Poisson cells

	// Create place cell layers ------------- #cells 1st Id randgen
	_ECPlaceCells		  = new ECLayer (this,  2000, 10000, pcspk);
	_DGPlaceCells		  = new DGLayer (this,  4000, 20000, pcspk);
	_CA3PlaceCells		  = new CA3Layer(this, 12000, 30000, pcspk);

	// Create interneuron cell layers ----------------- #cells 1st Id randgen
	_CA3AxoAxonicCells	  = new CA3AxoAxonicLayer   (this, 100, 50000, inspk);
	_CA3BasketCells		  = new CA3BasketLayer      (this, 100, 51000, inspk);
	_CA3BistratifiedCells = new CA3BistratifiedLayer(this, 500, 52000, inspk);
	_CA3OLMCells		  = new CA3OLMLayer         (this, 100, 53000, inspk);
}

// Set activity fractions for affected networks
void Mouse::setLayerActivities()
{
	_DGPlaceCells		  ->setFractionActive(0.05f);
	_CA3PlaceCells		  ->setFractionActive(0.18f);
}

// Connect the layers with the target
void Mouse::connectNetwork()
{

	UniformRandom*	urnet = networkRandomizer();

	// Initialize connection policy objects
	_PP = new ECPerforantPath(_ECPlaceCells, urnet);
	_AC = new CA3AssocCollaterals(_CA3PlaceCells, urnet);
	_MF = new DGMossyFibers(_DGPlaceCells, urnet);

	// Make the various connections for place cells.
	// Sequential order is used to force one connection per
	// afferent neuron and to allow identification of which
	// cells in the layer are actually afferents. This also
	// means that each cell has only one connection with
	// the target cell as long as sufficient cells were
	// initially allocated in the place cell layer.
	_PP->connectWith(targetCell(), ConnectionPolicy::sequentialOrder);
	_MF->connectWith(targetCell(), ConnectionPolicy::sequentialOrder);
	_AC->connectWith(targetCell(), ConnectionPolicy::sequentialOrder);

	// Now create interneuron connections using network randomizer.
	// Random order is used to permit multiple connections
	// per afferent neuron as is found experimentally.
	_INaxon		= new CA3AxonicInhibition	(_CA3AxoAxonicCells,	urnet);
	_INsoma		= new CA3SomaticInhibition	(_CA3BasketCells,		urnet);
	_INproximal	= new CA3ProximalInhibition	(_CA3BistratifiedCells,	urnet);
	_INdistal	= new CA3DistalInhibition	(_CA3OLMCells,			urnet);

	_INaxon		->connectWith(targetCell(), ConnectionPolicy::randomOrder);
	_INsoma		->connectWith(targetCell(), ConnectionPolicy::randomOrder);
	_INproximal	->connectWith(targetCell(), ConnectionPolicy::randomOrder);
	_INdistal	->connectWith(targetCell(), ConnectionPolicy::randomOrder);
}

// Print synpase creation statistics
void Mouse::printSynapseCounts()
{
	cerr<<endl<<"PP"<<endl;
	_PP->printSynapseCounts();

	cerr<<endl<<"MF"<<endl;
	_MF->printSynapseCounts();

	cerr<<endl<<"AC"<<endl;
	_AC->printSynapseCounts();

	cerr<<endl<<"INaxon"<<endl;
	_INaxon->printSynapseCounts();

	cerr<<endl<<"INsoma"<<endl;
	_INsoma->printSynapseCounts();

	cerr<<endl<<"INproximal"<<endl;
	_INproximal->printSynapseCounts();

	cerr<<endl<<"INdistal"<<endl;
	_INdistal->printSynapseCounts();

	cerr<<endl;
}

// When adding this to a controller, add the place cell layers also
void Mouse::addToController(Controller* cont)
{
	// Let superclasses handle adding this object
	MazeSubject::addToController(cont);

	// Add the place cell layers to the controller
	ECPlaceCells()			->addToController(cont);
	DGPlaceCells()			->addToController(cont);
	CA3PlaceCells()			->addToController(cont);
	CA3AxoAxonicCells()		->addToController(cont);
	CA3BasketCells()		->addToController(cont);
	CA3BistratifiedCells()	->addToController(cont);
	CA3OLMCells()			->addToController(cont);

	// If the network was built during initialization,
	// add the target cell to the controller. Otherwise not
	// since the whole idea of not building the network
	// is to avoid the time to simulate the target cell.
	if (buildNetwork() ) {
		targetCell()->addToController(cont);
	}
}
