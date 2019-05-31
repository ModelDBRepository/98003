// Provide classes for simulating cells of the hippocampal formation
//
// Copyright 2007 John L Baker. All rights reserved.
// This software is provided under the terms of the Open Source MIT License.
// See http://www.opensource.org/licenses/mit-license.php.
//
// File: hippoform_2003.cpp
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		3 March 2007
//
// Description:
//
// To simulate the effects of NMDA knockout in the CA3 portion of the
// hippocampus, it is necessary to generate inputs similar in form to 
// that which might be present for a mouse moving in a maze, typically 
// a Morris water maze or its dry equivalent. A single target cell is 
// simulated in a biologically realistic way. Other cells are simulated 
// phenomenologically based on location modulated firing rates.
//
// This header file declares classes used in the phenomenological 
// simulation. Only EC, DG, and CA3 are included here. Principal cells 
// for EC, DG, and CA3 are included. Interneurons are included for CA3 only.
//
// Interneuron firing rates and phase relationships come largely from
// Klausberger et al. Overall phase is based on the scheme in Fox et al.
// in which peak CA3 pyramidal activity occurs at 19-deg wrt theta.
// Klausberger puts peak CA1 pyramidal activity at 20-deg which is
// pretty much the same thing.
//
// If there is no feedback loop between target cell firing and interneuron
// rates then input spike trains should be reproducible. If interneuron rates
// are affected by target cell firing, then reproduciblity does not hold
// if target cell firing changes in even the slightest fashion. Even changing
// the target cell from 0.5ms time step to 1.0ms is enough to affect 
// afferent spike trains though hopefully not aggregate statistics..

#include "hippoform_baker_2003.h"
#include "mouse_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace BAKER_2003;



// --------------------------------------------------------------------
// CA3InterneuronLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3InterneuronLayer::CA3InterneuronLayer(
	Mouse* sub, int numCells, int numericId, UniformRandom* spkunif)
:	InterneuronLayer(sub, numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Save the maze subject for local reference in the right class
	_mouse = sub;

	// Locate input layers and save them here
	ECPlaceCells( sub->ECPlaceCells() );
	DGPlaceCells( sub->DGPlaceCells() );
	CA3PlaceCells(sub->CA3PlaceCells() );

	// Set layer specific defaults.Subclass will set
	// more meaningful values during construction.

	gammaFrequency(40*Hz);
	gammaAmplitude(0);
	gammaOffset(0);

	thetaFrequency(8*Hz);
	thetaAmplitude(0);
	thetaOffset(0);
	thetaPhase(0);

	ECFeedforward(0);
	DGFeedforward(0);
	CA3Feedforward(0);
	CA3Feedback(0,0*Hz,1*Hz);
	CA3FeedbackOngoingRate(0);
	CA3FeedbackTau(1*minute);
	baselineRate(0);

	// Set as of time to force special start-up handling
	_CA3FeedbackAsOfTime = InfinitePast;

	// The subclass must invoke allocateInterneurons() after
	// any further changes are made during construction.
}

CA3InterneuronLayer::~CA3InterneuronLayer() 
{
	// Remove the dependencies, if any
	if (_ECPlaceCells!=NULL) {
		_ECPlaceCells->removeSubscriber(this);
	}
	if (_DGPlaceCells!=NULL) {
		_DGPlaceCells->removeSubscriber(this);
	}
	if (_CA3PlaceCells!=NULL) {
		_CA3PlaceCells->removeSubscriber(this);
	}
}

// Save the EC layer and hook up as a dependent.
void CA3InterneuronLayer::ECPlaceCells(PlaceCellLayer* layer)
{
	_ECPlaceCells = layer;
	_ECPlaceCells->addSubscriber(this);
}

// Save the DG layer and hook up as a dependent.
void CA3InterneuronLayer::DGPlaceCells(PlaceCellLayer* layer)
{
	_DGPlaceCells = layer;
	_DGPlaceCells->addSubscriber(this);
}

// Save the CA3 layer and hook up as a dependent.
void CA3InterneuronLayer::CA3PlaceCells(PlaceCellLayer* layer)
{
	_CA3PlaceCells = layer;
	_CA3PlaceCells->addSubscriber(this);
}

// Set the baseline rate and get a new total rate
void CA3InterneuronLayer::baselineRate(Number r)
{
	_baselineRate = r;
	updateFrom(this,stateChange);
}

// Act on the end of the timestep to update feedback rate.
void CA3InterneuronLayer::timeStepEnded()
{
	// First let superclass do what it needs to at this point
	InterneuronLayer::timeStepEnded();

	// Force an update of the feedback rate now
	// Timing relations with the target cell are maintained
	// in the update function. We need only ensure that the
	// updates are done relatively frequently.
	updateCA3FeedbackRate();
}

// Respond to an update from an input place cell layer
// indicating that new rates are to be set. This gets done
// more often than really needed (once per input) but should
// have minimal performance impact.
void CA3InterneuronLayer::updateFrom(ModelComponent* mc, int reason)
{
	// Respond to state changes.
	if (reason==stateChange) {

		Number newFiringRate = baselineRate();

		// Sum up the contributions to the firing rate
		if (ECPlaceCells()!=NULL) {
			newFiringRate += ECFeedforward() * 
				ECPlaceCells()->totalFiringRate() / ECPlaceCells()->numCells();
		}

		if (DGPlaceCells()!=NULL) {
			newFiringRate += DGFeedforward() * 
				DGPlaceCells()->totalFiringRate() / DGPlaceCells()->numCells();
		}

		if (CA3PlaceCells()!=NULL) {
			newFiringRate += CA3Feedforward() * 
				CA3PlaceCells()->totalFiringRate() / CA3PlaceCells()->numCells();
		}

		if (mouse()->targetCell()!=NULL) {
			// Make sure rate is up to date
			updateCA3FeedbackRate(); 

			// Add the increment from target cell feedback
			Number fr = CA3FeedbackOngoingRate();
			Number frmax = CA3FeedbackMax();
			Number frh = CA3FeedbackFRh();
			Number frk = CA3FeedbackFRk();

			newFiringRate += frmax/(1+exp(-(fr-frh)/frk));
		}

		// Set the firing rate for this layer.
		// Make sure rate is positive in case activity
		// changes decrease firing rates.
		firingRate(plusval(newFiringRate));
	}

	// Respond to termination of an input layer.
	// Do not reference the input layer further.
	else if (reason==terminatedChange) {
		if (_ECPlaceCells==mc)  _ECPlaceCells = NULL;
		if (_DGPlaceCells==mc)  _DGPlaceCells = NULL;
		if (_CA3PlaceCells==mc) _CA3PlaceCells = NULL;
	}

	// Other change types are ignored.
	else {}
}

// Update the estimate of target cell firing rate
void CA3InterneuronLayer::updateCA3FeedbackRate()
{
	// Get the effective time of the rate which originates
	// in the target cell model and may differ from the
	// time in the model for this object.
	SimTime asOfTime	= mouse()->targetCell()->currentTime();
	SimTime h			= asOfTime - _CA3FeedbackAsOfTime;

	// See if a new update is warranted. If not stop now.
	// Note that we skip the update on the first time step
	// so that any initial value of feedback rate set in
	// subclasses will not be aged out prematurely.
	if (_CA3FeedbackAsOfTime==InfinitePast) {
		_CA3FeedbackAsOfTime = asOfTime;
		return;
	}
	if (h<=EpsilonTime) {
		return;
	}

	// Use an exponential smoother (decaying average)
	double exphtau = exp(-h/CA3FeedbackTau() );

	_CA3FeedbackOngoingRate *= exphtau;
	_CA3FeedbackOngoingRate += (1-exphtau) * mouse()->targetCell()->firingRate();

	// For next time, save the time when this rate was computed.
	_CA3FeedbackAsOfTime = asOfTime;
}

// Create a new interneuron
Interneuron* CA3InterneuronLayer::newInterneuron(int k)
{
	using namespace UOM;

	Interneuron* cell = InterneuronLayer::newInterneuron(k);

	// Set theta phase based on layer defaults
	cell->thetaPhase( thetaPhase() );

	return cell;
}



// --------------------------------------------------------------------
// CA3Inhibition class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3Inhibition::CA3Inhibition(
	CA3InterneuronLayer* source, 
	UniformRandom* randomizer)
:	ConnectionPolicy(source,randomizer)
{
	// Set addditional connection defaults
	numSomaSynapses(0);
	numISSynapses(0);
}

CA3Inhibition::~CA3Inhibition() {}

// Connect with the soma
void CA3Inhibition::connectWithSoma()
{
	Compartment* soma = _target->somaComp();
	int n= numSomaSynapses();
	createSynapses(soma,n);
	_totalSomaSynapses += n;
}

// Connect with the initial segment
void CA3Inhibition::connectWithIS()
{
	Compartment* IS = _target->ISComp();
	int n= numISSynapses();
	createSynapses(IS,n);
	_totalISSynapses += n;
}



// --------------------------------------------------------------------
// ECLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
ECLayer::ECLayer(MazeSubject* sub, int numCells, 
				 int numericId, UniformRandom* spkunif)
: PlaceCellLayer(sub,numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set layer specific defaults
	gammaAmplitude(0.5f);
	gammaFrequency(40*Hz);

	thetaAmplitude(1.0f);
	thetaFrequency(8*Hz);

	// Create all cells
	allocatePlaceCells();
}

ECLayer::~ECLayer() {}

// Allocate a new place cell centered at the given location
PlaceCell* ECLayer::newPlaceCell(int k)
{
	using namespace UOM;

	const bool useCoarseSpatialTuning = false;

	PlaceCell* cell = PlaceCellLayer::newPlaceCell(k);

	// Theta phase is taken from Quirk, Barry & Fox, 1992. 
	// Kamondi et al. has a CSD for CA1 with PP but it is ambiguous
	// because of pairing of sources  and sinks. One interpretation
	// is two population groups in PP, e.g. MEC and LEC or L-II and L-III.
	cell->thetaPhase(5*radiansPerDegree);

	// Set EC place field parameters based on either coarse
	// spatial tuning (Quirk et al.) or based on CA3-like tuning
	// of grid cells (Hafting et al.). The dissertation model uses
	// coarse tuning but subsequent article(s) may not.
	if (useCoarseSpatialTuning) {

		// Assume larger place fields (Quirk et al.)
		cell->meanFiringRate(0.8*Hz);
		cell->setDistanceEstSD(10*cm, 0.5f);

		// Set theta phase precession parameters.
		// There is evidence of phase precession as an input
		// to the hippocampus, but no hard data exists for EC.
		// Because the place fields in EC are larger than CA3,
		// the rate of precession is decreased accordingly.
		// If the precession rate is much higher than this,
		// theta modulation of EC population activity is no
		// longer apparent even when all cells are modulated.
		cell->precessionRate(3/cm*radiansPerDegree);
	}
	else {

		// Assume smaller place fields ala grid cells (Hafting et al.)
		cell->meanFiringRate(0.4*Hz);
		cell->setDistanceEstSD(3.5*cm, 0.15f);

		// Assume same precession rate as CA3 for grid cells.
		// Theta phase precession of grid cells has been reported
		// via abstract but not yet published.
		cell->precessionRate(10/cm*radiansPerDegree);
	}

	return cell;
}



// --------------------------------------------------------------------
// ECPerforantPath class body
// --------------------------------------------------------------------



// Constructors and destructor
ECPerforantPath::ECPerforantPath(
	PlaceCellLayer* source, 
	UniformRandom* randomizer)
:	ConnectionPolicy(source,randomizer)
{
	using namespace UOM;

	minAxonDist			(1*mm);
	maxAxonDist			(4*mm);
	maxLaminarDist		(500*micron);
	minLaminarDist		(350*micron);
	enpassantDist		(1*micron);
	dendriteSynDensity	(0.1/micron_2);
	synapseType			("PP_GluR");

	// Initialize weights.
	synapseWeight		(0.0f);
}
ECPerforantPath::~ECPerforantPath() {}



// --------------------------------------------------------------------
// DGLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
DGLayer::DGLayer(MazeSubject* sub, int numCells, 
				 int numericId, UniformRandom* spkunif)
: PlaceCellLayer(sub,numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set layer specific defaults
	gammaAmplitude(0.5f);
	gammaFrequency(40*Hz);

	thetaAmplitude(1.0f);
	thetaFrequency(8*Hz);

	// Create all cells
	allocatePlaceCells();
}

DGLayer::~DGLayer() {}

// Allocate a new place cell centered at the given location
PlaceCell* DGLayer::newPlaceCell(int k)
{
	using namespace UOM;

	PlaceCell* cell = PlaceCellLayer::newPlaceCell(k);

	// Set params for layer specific place field properties
	cell->meanFiringRate(0.2*Hz);
	cell->inactiveFiringRate(0.01*Hz);
	cell->setDistanceEstSD(2*cm, 0.1f);

	// Set mean theta phase based on Fox, Wolfson, and Ranch 1986.
	cell->thetaPhase(296*radiansPerDegree);		// 83-deg ahead of CA3

	// Set theta phase precession paramters.
	// See Skaggs et al. for data but there is no obviously
	// correct value for theta phase precession.
	// The value from CA3 is used, in which case
	// phase implies the same distance from the
	// center of the place field as for CA3.
	cell->precessionRate(10/cm*radiansPerDegree);

	return cell;
}




// --------------------------------------------------------------------
// DGMossyFiber class body
// --------------------------------------------------------------------



// Constructors and destructor
DGMossyFibers::DGMossyFibers(
	PlaceCellLayer* source, 
	UniformRandom* randomizer)
:	ConnectionPolicy(source,randomizer)
{
	using namespace UOM;

	minAxonDist			(1.0*mm);
	maxAxonDist			(1.1*mm);
	maxLaminarDist		(100*micron);
	minLaminarDist		(20*micron);
	enpassantDist		(1*micron);
	dendriteSynDensity	(1.5e-2/micron_2);
	synapseType			("MF_GluR");
	synapseWeight		(1.0f);
}
DGMossyFibers::~DGMossyFibers() {}



// --------------------------------------------------------------------
// CA3Layer class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3Layer::CA3Layer(MazeSubject* sub, int numCells, 
				   int numericId, UniformRandom* spkunif)
: PlaceCellLayer(sub,numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set layer specific defaults
	gammaFrequency(40*Hz);
	gammaAmplitude(0.20f);

	thetaFrequency(8*Hz);
	thetaAmplitude(0.75f);

	// Create all cells
	allocatePlaceCells();
}

CA3Layer::~CA3Layer() {}

// Allocate a new place cell centered at the given location
PlaceCell* CA3Layer::newPlaceCell(int k)
{
	using namespace UOM;

	PlaceCell* cell = PlaceCellLayer::newPlaceCell(k);

	// Set params for layer specific place field properties
	cell->inactiveFiringRate(0.05*Hz);
	cell->meanFiringRate(0.4*Hz);
	cell->setDistanceEstSD(3.5*cm, 0.15f);

	// Set mean theta phase based on Fox, Wolfson, and Ranch 1986.
	// This is consistent with CA1 results in Klausberger et al.
	cell->thetaPhase(19*radiansPerDegree);

	// Set theta phase precession parameters.
	// This is consistent with CA3 phase precession
	// in O'Keefe and Reece 1993 Hippocampus 3, 317-330.
	// Just in case it's not obvious, 1000 deg/m = 10 deg/cm.
	cell->precessionRate(10/cm*radiansPerDegree);

	return cell;
}



// --------------------------------------------------------------------
// CA3AssocCollateral class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3AssocCollaterals::CA3AssocCollaterals(
	PlaceCellLayer* source, 
	UniformRandom* randomizer)
:	ConnectionPolicy(source,randomizer)
{
	using namespace UOM;

	maxAxonDist			(3*mm);
	maxLaminarDist		(350*micron);
	minLaminarDist		(100*micron);
	maxBasalDist		(0*micron);
	minBasalDist		(-999*micron);
	enpassantDist		(1*micron);
	dendriteSynDensity	(0.1/micron_2);
	synapseType			("AC_GluR");

	// Initialize weights.
	synapseWeight		(0.0f);
}
CA3AssocCollaterals::~CA3AssocCollaterals() {}



// --------------------------------------------------------------------
// CA3AxoAxonicLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3AxoAxonicLayer::CA3AxoAxonicLayer(
	Mouse* sub, int numCells, int numericId, UniformRandom* spkunif)
:	CA3InterneuronLayer(sub,numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set firing rates based on Klausberger et al.
	// Variations involve making rates dependent on
	// activity in other layers.
	Number			nominalRate = 17*Hz;

	baselineRate	(nominalRate);			// Normal initial firing rate	
	inactiveRate	(0.5*nominalRate);		// Novel region firing rate

	// Set amplitude and phase roughly in line with Klausberger et al.
	thetaAmplitude(1.00f);
	thetaPhase(185*radiansPerDegree);		// 194-deg ahead of CA3
	allocateInterneurons();
}
CA3AxoAxonicLayer::~CA3AxoAxonicLayer() {}



// --------------------------------------------------------------------
// CA3AxonicInhibition class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3AxonicInhibition::CA3AxonicInhibition(
	CA3InterneuronLayer* source, 
	UniformRandom* randomizer)
:	CA3Inhibition(source,randomizer)
{
	numISSynapses		(200);

	synapseType			("GABAa");
	synapseWeight		(1.0f);
}
CA3AxonicInhibition::~CA3AxonicInhibition() {}



// --------------------------------------------------------------------
// CA3BasketLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3BasketLayer::CA3BasketLayer(
	Mouse* sub, int numCells, int numericId, UniformRandom* spkunif)
:	CA3InterneuronLayer(sub, numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set firing rates based on Klausberger et al.
	// Variations involve making rates dependent on
	// activity in other layers.
	Number			nominalRate = 8*Hz;

	baselineRate	(nominalRate);			// Normal initial firing rate	
	inactiveRate	(0.5*nominalRate);		// Novel region firing rate

	// Couple with feedback from target cell firing
	// Peak firing rate is estimated based on firing during
	// sharp wave intervals (data from Klausberger et al.)
	CA3Feedback(3*nominalRate,3.0*Hz,0.5*Hz);	// Dynamic rate adjustment
	CA3FeedbackTau(60*sec);						// Time constant to estimate rate

	// Set amplitude and phase roughly in line with Klausberger et al.
	thetaAmplitude(0.65f);
	thetaPhase(271*radiansPerDegree);		// 108-deg ahead of CA3
	allocateInterneurons();
}
CA3BasketLayer::~CA3BasketLayer() {}

// Allocate a new interneuron.. This can be extended by 
// subclasses to reflect layer-specific properties.
Interneuron* CA3BasketLayer::newInterneuron(int k)
{
	using namespace UOM;

	SawToothInterneuron* cell = new SawToothInterneuron(model(), this);

	// Set theta phase based on layer defaults
	cell->thetaPhase( thetaPhase() );

	// Hard code the phase at which firing is minimum
	cell->contraThetaPhase( 0*radiansPerDegree );

	return cell;
}



// --------------------------------------------------------------------
// CA3SomaticInhibition class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3SomaticInhibition::CA3SomaticInhibition(
	CA3InterneuronLayer* source, 
	UniformRandom* randomizer)
:	CA3Inhibition(source,randomizer)
{
	using namespace UOM;

	numSomaSynapses		(200);

	synapseType			("GABAa");
	synapseWeight		(1.0f);
}
CA3SomaticInhibition::~CA3SomaticInhibition() {}



// --------------------------------------------------------------------
// CA3BistratifiedLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3BistratifiedLayer::CA3BistratifiedLayer(
	Mouse* sub, int numCells, int numericId, UniformRandom* spkunif)
:	CA3InterneuronLayer(sub, numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set firing rates based on Klausberger et al.
	// Variations involve making rates dependent on
	// activity in other layers.
	Number			nominalRate = 6*Hz;

	baselineRate	(nominalRate);			// Normal initial firing rate	
	inactiveRate	(0.5*nominalRate);		// Novel region firing rate

	// Couple with feedback from target cell firing
	// Peak firing rate is estimated based on firing during
	// sharp wave intervals (data from Klausberger et al.)
	CA3Feedback(3*nominalRate,3.0*Hz,0.5*Hz);	// Dynamic rate adjustment
	CA3FeedbackTau(60*sec);						// Time constant to estimate rate

	// Set amplitude and phase roughly in line with Klausberger et al.
	// See Klausberger et al. 2003 table 2 for various values by cell.
	thetaAmplitude(0.95f);
	thetaPhase(359*radiansPerDegree);		// 20-deg ahead of CA3	

	// Build the layer
	allocateInterneurons();
}
CA3BistratifiedLayer::~CA3BistratifiedLayer() {}



// --------------------------------------------------------------------
// CA3ProximalInhibition class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3ProximalInhibition::CA3ProximalInhibition(
	CA3InterneuronLayer* source, 
	UniformRandom* randomizer)
:	CA3Inhibition(source,randomizer)
{
	using namespace UOM;

	maxLaminarDist		(350*micron);
	minLaminarDist		(-999*micron);
	enpassantDist		(1*micron);
	dendriteSynDensity	(0.025/micron_2);
	synapseType			("GABAa","GABAb");

	synapseWeight		(0.33f, 0.33f);
}
CA3ProximalInhibition::~CA3ProximalInhibition() {}



// --------------------------------------------------------------------
// CA3OLMLayer class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3OLMLayer::CA3OLMLayer(
	Mouse* sub, int numCells, int numericId, UniformRandom* spkunif)
:	CA3InterneuronLayer(sub, numCells,numericId,spkunif) 
{
	using namespace UOM;

	// Set firing rates based on Klausberger et al.
	// Variations involve making rates dependent on
	// activity in other layers.
	Number			nominalRate = 5*Hz;

	baselineRate	(nominalRate);			// Normal initial firing rate	
	inactiveRate	(0.5*nominalRate);		// Novel region firing rate

	// Set amplitude and phase roughly in line with Klausberger et al.
	thetaAmplitude(0.85f);
	thetaPhase(19*radiansPerDegree);		// Same as CA3
	allocateInterneurons();
}
CA3OLMLayer::~CA3OLMLayer() {}



// --------------------------------------------------------------------
// CA3DistalInhibition class body
// --------------------------------------------------------------------



// Constructors and destructor
CA3DistalInhibition::CA3DistalInhibition(
	CA3InterneuronLayer* source, 
	UniformRandom* randomizer)
:	CA3Inhibition(source,randomizer)
{
	using namespace UOM;

	maxLaminarDist		(999*micron);
	minLaminarDist		(350*micron);
	enpassantDist		(1*micron);
	dendriteSynDensity	(0.025/micron_2);

	synapseType			("GABAas", "GABAb");

	synapseWeight		(0.33f, 0.33f);
}
CA3DistalInhibition::~CA3DistalInhibition() {}
