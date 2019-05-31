// Provide classes for simulating phonemological place cells during maze movement
//
// Copyright 2007 John L Baker. All rights reserved.
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: placecell_baker_2003.h
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		6 March 2007
//
// Description:
//
// Place cells simulated here are purely phenomenological. They are simulated
// using PoissonNeurons with random firing patterns, which can be modulated based
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
// References:
//
// Buzsaki G (2002) Theta oscillations in the hippocampus. Neuron 33: 325-340.
//
// Csicsvari J, Jamieson B, Wise KD, Buzsáki G (2003) Mechanisms of gamma 
// oscillations in the hippocampus of the behaving rat. Neuron 37: 311-322.
//
// Dayan P, Abbott LF (2001) Theoretical Neuroscience: computational and
// mathematical modeling of neural systems. MIT Press, Cambridge MA.
//
// Hartley T, Burgess N, Lever C, Cacucci F, O'Keefe J (2000) Modeling
// place fields in terms of the cortical inputs to the hippocampus.
// Hippocampus 10: 369-379.
//
// Kali S, Dayan P (2000) The involvement of recurrent connections in area
// CA3 in establishing the properties of place fields: a model.
// Journal of Neuroscience 20: 7463-7477.
//
// O'Keefe J, Burgess N (1996) Geometric determinants of the place fields
// of hippocampal neurons. Nature 381: 425-428.
// 
// O'Keefe J, Dostrovsky J (1971) The hippocampus as a spatial map: 
// preliminary evidence from unit activity in the freely moving rat. 
// Brain Research 34: 171-175.
//
// Wallenstein GV, Hasselmo ME (1997) GABAergic modulation of hippocampal
// population activity: sequence learning, place field development, and the
// phase precession effect. Journal of Neurophysiology 78: 393-408.



// --------------------------------------------------------------------
// Only include the definitions in this header once
// --------------------------------------------------------------------

#ifndef __PLACECELL_BAKER_2003_H_
#define __PLACECELL_BAKER_2003_H_


// --------------------------------------------------------------------
// MICROSOFT SPECIFIC DECLARATIONS
// --------------------------------------------------------------------
#ifdef WIN32

// Disable warning C4786: symbol greater than 255 character,
#pragma warning( disable: 4786)

#endif
// --------------------------------------------------------------------
// END OF MICROSOFT SPECIFIC DECLARATIONS
// --------------------------------------------------------------------


#include "bnsf.h"
#include "bnsf_liaf.h"
#include "maze_baker_2003.h"
#include "subject_baker_2003.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace BAKER_2003 {

	// ----------------------------------------------------------------
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ----------------------------------------------------------------

	class PlaceCell;
	class Interneuron;
	class CellLayer;
		class PlaceCellLayer;
		class InterneuronLayer;
	class ConnectionPolicy;

	// ----------------------------------------------------------------
	// CLASS:	PlaceCell
	// EXTENDS:	PoissonNeuron
	// DESC:	Defines a cell with place field response
	//			with spike generation at random at a rate
	//			determined by the local geometry of the maze.
	// RESP:
	//		1.	Set average firing rate based on location.
	//		2.	Modulate firing rate based on gamma and theta rhythms.
	//		3.	Modulate firing phase based on direction and distance.
	//
	// NOTES:	Place field characteristics are determine by a
	//			two-dimensional Gaussian density centered at the
	//			place field center. The distance error model of the
	//			layer is used to provide the standard deviation of
	//			each coordinate in the density. 
	//
	//			A mean firing rate is used to scale the density
	//			to an actual firing rate. The calculation is an
	//			estimate that does not take into account edge effects.
	//			or the finite size of any given maze.
	//
	//			This class supports only a single place field center,
	//			but includes framework support for multiple centers
	//			in which the relevant center is determined from the
	//			current position.
	// ----------------------------------------------------------------

	class PlaceCell : public PoissonNeuron {

	public:

		// Constructors and destructor
		PlaceCell(Model* m, PlaceCellLayer* pclayer);
		virtual ~PlaceCell();

		// Accessors
		virtual PlaceCellLayer*		layer() { return _layer; }
		virtual MazeSubject*		subject();
		virtual Maze*				maze();

		// Access a flag indicating whether or not the cell is active
		// in the current environment/context.
		inline  bool				isActive() { return _isActive; }
		inline  bool				isInactive() { return !_isActive; }
		virtual void				isActive(bool b);

		// Access an idealized mean firing rate for an active cell. 
		// This assumes spatial uniformity over a large (i.e. effectively 
		// infinite) plane. Depending on other parameters, this may differ
		// from the average firing rate in any particular enclosure.
		inline  Number				meanFiringRate() { return _meanFiringRate; }
		virtual void				meanFiringRate(Number r) { _meanFiringRate=r; }

		// Access the firing rate used when the cell is inactive.
		// This rate is not dependent on location.
		inline  Number				inactiveFiringRate() { return _inactiveFiringRate; }
		virtual void				inactiveFiringRate(Number r) { _inactiveFiringRate=r; }

		// Access the coordinates of the field center location
		inline  Number				centerX() { return _centerX; }
		inline  Number				centerY() { return _centerY; }
		virtual void				setCenter(Number x, Number y);

		// Access the theta phase offset (in radians) for this cell.
		// This offset is the phase at which firing probability is maximal.
		inline  Number				thetaPhase() {return _thetaPhase; }
		virtual void				thetaPhase(Number tp) { _thetaPhase = tp; }

		// Access the rate of theta phase precession in radians per unit length
		inline  Number				precessionRate() { return _precessionRate; }
		virtual void				precessionRate(Number r) { _precessionRate=r; }

		// Access parameters used to estimate standard deviation of estimated distances
		inline  Number				distanceEstSD() { return _distanceEstSD; }
		inline	Number				propDistanceEstSD() { return _propDistanceEstSD; }
		virtual void				setDistanceEstSD(Number sdall, Number sdprop);

		// Access the firing rate at the place field center
		inline  Number				peakRate();

		// Recompute firing rate based on current position and return computed value
		virtual Number				updateFiringRate();

		// Verify that all properties are specified prior to startup
		virtual void				simulationStarted();

	protected:
		PlaceCellLayer*				_layer;				// associated place cell layer

		Number						_meanFiringRate;	// desired mean firing rate
		Number						_fieldRadius;		// Nominal radius of best estimate
		Number						_distanceEstSD;		// std dev independent of distance
		Number						_propDistanceEstSD;	// std dev proportional to distance

		Number						_thetaPhase;		// phase offset for this cell
		Number						_precessionRate;	// rate of theta phase precession

		Number						_centerX;			// X coord of place field center
		Number						_centerY;			// Y coord of place field center

		Number						_Xvar;				// variance in X coord
		Number						_Yvar;				// variance in Y coord
		Number						_hx;				// local orientation x component
		Number						_hy;				// local orientation y component

		Number						_densityMultiplier;		// multiplier to normalize density
		Number						_inactiveFiringRate;	// rate used when not active
		bool						_isActive;				// rate is location dependent

		// Return a location in local coordinates. These are distance measures
		// in different directions as appropriate for this location in the maze.
		// By default, distances are returned using the directions north, south, 
		// east, and west corresponding to vectors as rotated by the local orientation.
		// This delegated to the maze, but subclasses can choose otherwise.
		virtual void				localCoordinates(
			Number					locX,			// Current location X coord
			Number					locY,			// Current location Y coord
			Number&					northDist,		// Boundary dist along (0,1)
			Number&					southDist,		// Boundary dist along (0,-1)
			Number&					eastDist,		// Boundary dist along (1,0)
			Number&					westDist);		// Boundary dist along (-1,0)

		// Provide a heading defined by local conditions in the maze such as
		// the direction of the nearest boundary point. The purpose is to provide
		// an orientation for a local coordinate system for this part of the maze.
		// Default is a constant vector (1,0). Subclass should override as needed.
		// This delegated to the maze, but subclasses can choose otherwise.
		virtual void				localOrientation(
			Number					locX,			// Current location X coord
			Number					locY,			// Current location Y coord
			Number&					hx,				// Unit vector x coord (output)
			Number&					hy);			// Unit vector y coord (output)

		// Get the location of the nearest place field center. By default
		// there is only one center, but subclasses may support more than one.
		virtual void				locateFieldCenter(
			Number					locX,			// Current location X coord
			Number					locY,			// Current location Y coord
			Number&					ctrX,			// Center location X coord
			Number&					ctrY);			// Center location X coord

		// Set the properties of the place field based on current values
		// for center, layer properties, and the maze subject.
		virtual void				updateProperties();

		// Return the probability of selecting a spike among those generated
		// at peakFiringRate and at the time of the next time.
		virtual Number				spikeSelectionProbability();
	};

	// ----------------------------------------------------------------
	// CLASS:	Interneuron
	// EXTENDS:	PoissonNeuron
	// DESC:	Defines a cell where spike rate is controlled by
	//			gamma and theta rhythm parameters.
	// RESP:
	//		1.	Get firing rate from common layer value.
	//		2.	Modulate firing rate based on gamma and theta rhythms.
	// ----------------------------------------------------------------

	class Interneuron : public PoissonNeuron {

	public:

		// Constructors and destructor
		Interneuron(Model* m, InterneuronLayer* inlayer);
		virtual ~Interneuron();

		// Accessors
		virtual InterneuronLayer*	layer() { return _layer; }

		// Get current firing rate directly from the layer.
		// Function to set rate is a pass
		virtual Number				firingRate();

		// Access the theta phase offset (in radians) for this cell.
		// This offset is the phase at which firing probability is maximal.
		inline  Number				thetaPhase() {return _thetaPhase; }
		virtual void				thetaPhase(Number tp) { _thetaPhase = tp; }

	protected:
		InterneuronLayer*			_layer;				// associated place cell layer
		Number						_thetaPhase;		// phase offset for this cell

		// Return the probability of selecting a spike among those generated
		// at peakFiringRate and at the time of the next time.
		virtual Number				spikeSelectionProbability();
	};

	// ----------------------------------------------------------------
	// CLASS:	SawToothInterneuron
	// EXTENDS:	Interneuron
	// DESC:	Defines a cell where spike rate is controlled by gamma
	//			and theta rhythm parameters where the firing rate
	//			varies in a saw-tooth wave, that is separate time-
	//			linear rising and falling regions in the theta cycle.
	// RESP:
	//		1.	Modulate firing rate based on gamma and theta rhythms.
	// ----------------------------------------------------------------

	class SawToothInterneuron : public Interneuron {

	public:

		// Constructors and destructor
		SawToothInterneuron(Model* m, InterneuronLayer* inlayer);
		virtual ~SawToothInterneuron();

		// Access the theta phase offset (in radians) for this cell.
		// This offset is the phase at which firing probability is minimal.
		inline  Number				contraThetaPhase() {return _contraThetaPhase; }
		virtual void				contraThetaPhase(Number tp) { _contraThetaPhase = tp; }

	protected:
		Number						_contraThetaPhase; // phase for min firing rate

		// Return the probability of selecting a spike among those generated
		// at peakFiringRate and at the time of the next time.
		virtual Number				spikeSelectionProbability();
	};

	// ----------------------------------------------------------------
	// CLASS:	CellLayer
	// EXTENDS:	ModelComponent
	// DESC:	Abstract class Defining a collection of place cells 
	//			or interneurons.
	// RESP:
	//		1.	Define protocol for common access to layer cells.
	//		2.	Provide vector of subscribers for changed/update
	//		3.	Provide parameters for gamma and theta oscillations
	//		4.	Set numeric identifiers for cells
	//		5.	Create a common model and clock solver for the layer
	// ----------------------------------------------------------------

	class CellLayer : public ModelComponent {

	public:

		// Constructors and destructor
		CellLayer(
			MazeSubject*			sub,				// maze subject
			int						numericId=0,		// starting numeric id
			UniformRandom*			spkunif=NULL);		// spike generation randomizer 

		virtual ~CellLayer();

		// Accessors
		virtual ModelComponentVector* subscribers() { return &_subscribers; }
		virtual ODESolver*			solver() { return model()->solver(); }
		virtual SimTime				timeStep() { return solver()->timeStep(); }
		virtual void				timeStep(SimTime h) { solver()->timeStep(h); }

		inline  MazeSubject*		subject() { return _subject; }
		inline  Maze*				maze() { return subject()->maze(); }

		// Assign numeric identifier associated place cells
		// Place cells are assigned sequentially higher numbers.
		// The cell layer itself is assigned the first identifier.
		virtual int					numericIdentifier() { return _numericIdentifier; }
		virtual void				numericIdentifier(int n);

		// Locate a spatial region overwhich layer cells are inactive
		inline  SpatialRegion*		inactiveRegion() { return _inactiveRegion; }
		virtual void				inactiveRegion(SpatialRegion* r) { _inactiveRegion = r; }

		// Add/remove a probe for all cells
		virtual void				addProbeToAll(Probe* pr);
		virtual void				removeProbeFromAll(Probe* pr);

		// Accessors of gamma and theta rhythm modulation
		// Each rhythm is treated as a multiplier on the firing rate
		// with the formula multiplier(t)=max(0,1+a*cos(2*pi*f*t-o)) where
		// a=amplitude, f=freq, o=offset, and t=time.
		inline	Number				gammaFrequency() { return _gammaFrequency; }
		virtual void				gammaFrequency(Number f) { _gammaFrequency = f; }

		inline	Number				gammaAmplitude() { return _gammaAmplitude; }
		virtual void				gammaAmplitude(Number a) { _gammaAmplitude = a; }

		inline  SimTime				gammaOffset() { return _gammaOffset; }
		virtual void				gammaOffset(SimTime t) { _gammaOffset = t; }

		inline	Number				thetaFrequency() { return _thetaFrequency; }
		virtual void				thetaFrequency(Number f) { _thetaFrequency = f; }

		inline	Number				thetaAmplitude() { return _thetaAmplitude; }
		virtual void				thetaAmplitude(Number a) { _thetaAmplitude = a; }

		inline  SimTime				thetaOffset() { return _thetaOffset; }
		virtual void				thetaOffset(SimTime t) { _thetaOffset = t; }

		// Indicate that this object holds no external state variables.
		virtual int					numStateVar() { return 0; }

		// Subclass responsibilities ----------------------------------

		virtual int					numCells() = 0;
		virtual SpikingNeuron*		cell(int n) = 0;	// return cell n (n=1 is the first)

	protected:
		ModelComponentVector		_subscribers;		// subscribers for state change
		MazeSubject*				_subject;			// associated maze subject
		int							_numericIdentifier;	// starting id for this layer

		SpatialRegion*				_inactiveRegion;	// region of layer inactivity

		Number						_gammaFrequency;	// Frequency of gamma rhythms
		Number						_gammaAmplitude;	// Amplitude of the modulation
		SimTime						_gammaOffset;		// Offset of rhythm in time
		
		Number						_thetaFrequency;	// Frequency of theta rhythms
		Number						_thetaAmplitude;	// Amplitude of the modulation
		SimTime						_thetaOffset;		// Offset of rhythm in time
	};

	// ----------------------------------------------------------------
	// CLASS:	PlaceCellLayer
	// EXTENDS:	CellLayer
	// DESC:	Defines a collection of place cells of
	//			similar characteristics.
	// RESP:
	//		1.	Create a common model and clock solver
	//		2.	Allocate and initialize place cells
	//		3.	Connect place cells with a target neuron
	//		4.	Invoke recalculation of spike rates upon location change
	//
	// NOTES:	The model allocated here is provided to place cells so that
	//			the layer object and place cells will function synchronously.
	//
	//			Because different layers have different place field
	//			characteristics, subclasses should provide any unique
	//			place cell parameters required.
	// ----------------------------------------------------------------

	class PlaceCellLayer : public CellLayer {

	public:

		// Constructors and destructor
		PlaceCellLayer(
			MazeSubject*			sub,			// maze subject
			int						numCells,		// number of cells to allocate
			int						numericId=0,	// starting numeric id
			UniformRandom*			spkunif=NULL);	// spike generation randomizer

		virtual ~PlaceCellLayer();

		// Accessors
		inline  Number				totalFiringRate() { return _totalFiringRate; }
		virtual int					numCells() { return _numPlaceCells; }
		inline  Number				fractionActive() { return _fractionActive; }
	
		// Return a place cell by number starting with the first as n=1
		virtual PlaceCell*			placeCell(int n);
		virtual SpikingNeuron*		cell(int n) { return placeCell(n); }

		// Respond to an update in state of the maze subject
		virtual void				updateFrom(ModelComponent* mc, int reason);

		// Set place cell centers randomly within the maze. This can be 
		// overridden if random assigment of centers is not appropriate.
		virtual void				setPlaceFieldCenters();

		// For a subset of cells, set them active and all others inactive.
		// The subset is selected randomly but the number set active
		// is predetermined by the activeFraction.
		virtual void				setFractionActive(Number activeFraction);

		// Set a group of place cells active or inactive. From and to are numbered
		// from zero to match cells() and placeCell() accessors.
		virtual void				setActivity(bool activeState, int from, int to);

		// Disable theta phase precession in all cells of the layer.
		// Previous phase precession setting are lost when this is done.
		virtual void				disableThetaPhasePrecession();

		// Print field centers to an external file or, if none, stdout.
		// Output is in units of centimeters.
		virtual void				printPlaceFieldCenters(char* pathName=NULL);

	protected:
		int							_numPlaceCells;		// number of place cells
		PlaceCell**					_placeCells;		// array of place cell pointers
		Number						_totalFiringRate;	// total as of last update
		Number						_fractionActive;	// fraction chosen as active

		// Subclass responsibilities ----------------------------------

		// Allocate the place cells and set their parameters based on
		// the current subject and maze. This can only be done once.
		// Since virtual functions are not inherited during construction,
		// subclasses will need to invoke this when all parameters are set.
		virtual void				allocatePlaceCells();

		// Allocate and initialize a place cell. Subclasses 
		// can extend this for customized initializations. 
		// k is the number of the cell in this layer.
		virtual PlaceCell*			newPlaceCell(int k);
	};

	// ----------------------------------------------------------------
	// CLASS:	InterneuronLayer
	// EXTENDS:	CellLayer
	// DESC:	Defines a collection of place cells of
	//			similar characteristics.
	// RESP:
	//		1.	Create a common model and clock solver
	//		2.	Allocate and initialize place cells
	//		3.	Connect place cells with a target neuron
	//		4.	Invoke recalculation of spike rates upon location change
	//
	// NOTES:	The model allocated here is provided to place cells so that
	//			the layer object and place cells will function synchronously.
	//
	//			Because different layers have different place field
	//			characteristics, subclasses should provide any unique
	//			place cell parameters required.
	// ----------------------------------------------------------------

	class InterneuronLayer : public CellLayer {

	public:

		// Constructors and destructor
		InterneuronLayer(
			MazeSubject*			sub,			// maze subject
			int						numCells,		// number of cells to allocate
			int						numericId=0,	// starting numeric id
			UniformRandom*			spkunif=NULL);	// spike generation randomizer
		virtual ~InterneuronLayer();

		// Accessors
		inline  Number				firingRate() { return _firingRate; }
		virtual void				firingRate(Number r) { _firingRate = r; }

		inline  Number				inactiveRate() { return _inactiveRate; }
		virtual void				inactiveRate(Number r) { _inactiveRate = r; }

		virtual int					numCells() { return _numInterneurons; }

		// Return a cell by number starting with the first as n=1
		virtual SpikingNeuron*		cell(int n) { return interneuron(n); }
		virtual Interneuron*		interneuron(int n);	// n=1 returns first cell

	protected:
		int							_numInterneurons;	// number of interneurons
		Interneuron**				_interneurons;		// array of interneurons
		Number						_firingRate;		// interneuron base firing rate
		Number						_inactiveRate;		// rate when layer is inactive

		// Subclass responsibilities ----------------------------------

		// Allocate the interneurons and set their parameters based on
		// the current subject and maze. This can only be done once.
		// Since virtual functions are not inherited during construction,
		// subclasses will need to invoke this when all parameters are set.
		virtual void				allocateInterneurons();

		// Allocate and initialize an interneuron. Subclasses 
		// can extend this for customized initializations. 
		// k is the number of the cell in this layer.
		virtual Interneuron*		newInterneuron(int	k);
	};

	// ----------------------------------------------------------------
	// CLASS:	ConnectionPolicy
	// EXTENDS:	none
	// DESC:	Defines a policy for connecting place cells
	//			interneurons with a target neuron.
	// RESP:
	//		1.	Know afferent place cell layer
	//		2.	Know target neuron
	//		3.	Know boundaries of connection pathway
	//		4.	Know synapse types used in the connection
	//		5.	Connect cell from the place cell layer with the target
	//
	// NOTES:	Parameter values are supplied by subclass constructors.
	//
	//			The place cell layer subject normally supplies the source
	//			of random numbers used in making connections.
	// ----------------------------------------------------------------

	class ConnectionPolicy {

	public:

		// Typedefs and enums for this class
		enum ConnectionSequence { 
			sequentialOrder,			// use source cells in sequence
			randomOrder};				// select source cells at random

		// Constructors and destructor
		ConnectionPolicy(
			CellLayer*				source,			// afferent cell layer
			UniformRandom*			randomizer);	// source of random numbers
		
		~ConnectionPolicy();

		// Accessors for parameters
		inline  UniformRandom*		uniformRandom() { return _uniformRandom; }
		virtual void				uniformRandom(UniformRandom* unif) { _uniformRandom = unif; }

		inline  Number				minLaminarDist() { return _minLaminarDist; }
		virtual void				minLaminarDist(Number x) { _minLaminarDist=x; }

		inline  Number				maxLaminarDist() { return _maxLaminarDist; }
		virtual void				maxLaminarDist(Number x) { _maxLaminarDist=x; }

		inline  Number				minBasalDist() { return _minBasalDist; }
		virtual void				minBasalDist(Number x) { _minBasalDist=x; }

		inline  Number				maxBasalDist() { return _maxBasalDist; }
		virtual void				maxBasalDist(Number x) { _maxBasalDist=x; }

		inline  Number				enpassantDist() { return _enpassantDist; }
		virtual void				enpassantDist(Number x) { _enpassantDist=x; }

		inline  Number				minAxonDist() { return _minAxonDist; }
		virtual void				minAxonDist(Number x) { _minAxonDist=x; }

		inline  Number				maxAxonDist() { return _maxAxonDist; }
		virtual void				maxAxonDist(Number x) { _maxAxonDist=x; }

		inline  Number				somaSynDensity() { return _somaSynDensity; }
		virtual void				dendriteSynDensity(Number x) { _dendriteSynDensity=x; }

		inline  Number				ISSynDensity() { return _ISSynDensity; }
		virtual void				somaSynDensity(Number x) { _somaSynDensity=x; }

		inline  Number				dendriteSynDensity() { return _dendriteSynDensity; }
		virtual void				ISSynDensity(Number x) { _ISSynDensity=x; }

		inline  Number				synapseWeight() { return _synapseWeight; }
		virtual void				synapseWeight(Number x) { _synapseWeight=x; }

		inline  Number				synapseWeightAlt() { return _synapseWeightAlt; }
		virtual void				synapseWeightAlt(Number x) { _synapseWeightAlt=x; }

		inline  TokenId				synapseType() { return _synapseType; }
		virtual void				synapseType(TokenId x) { _synapseType=x; }
		virtual void				synapseType(char* s) { _synapseType=token(s); }

		inline  TokenId				synapseTypeAlt() { return _synapseTypeAlt; }
		virtual void				synapseTypeAlt(TokenId x) { _synapseTypeAlt=x; }
		virtual void				synapseTypeAlt(char* s) { _synapseTypeAlt=token(s); }

		// Set multiple values in one function invocation as a convenience
		virtual void				synapseWeight(Number w1, Number w2)
		{ synapseWeight(w1); synapseWeightAlt(w2); }

		virtual void				synapseType(TokenId st1, TokenId st2)
		{ synapseType(st1); synapseTypeAlt(st2); }

		virtual void				synapseType(char* st1, char* st2)
		{ synapseType(st1); synapseTypeAlt(st2); }

		// Accessors for statistics
		inline  int					totalSomaSynapses() { return _totalSomaSynapses; }
		inline  int					totalISSynapses() { return _totalISSynapses; }
		inline  int					totalAxonSynapses() { return _totalAxonSynapses; }
		inline  int					totalApicalSynapses() { return _totalApicalSynapses; }
		inline  int					totalBasalSynapses() { return _totalBasalSynapses; }

		// Return total of all connected synapses
		virtual int					totalSynapses();

		// Connect the source layer with the target cell
		virtual void				connectWith(
			MorphologicalNeuron*	target,
			ConnectionSequence		connSeq = sequentialOrder);

		// Print current synapse counts for debug
		virtual void				printSynapseCounts();

	protected:

		// Supplied parameters
		CellLayer*					_layer;				// source layer
		UniformRandom*				_uniformRandom;		// connection randomizer
		MorphologicalNeuron*		_target;			// target neuron

		// Connection statistics (for reporting and debugging)
		int							_totalSomaSynapses;
		int							_totalISSynapses;
		int							_totalAxonSynapses;
		int							_totalApicalSynapses;
		int							_totalBasalSynapses;

		// Define the minimum and maximum laminar distance for the
		// layer in terms of the target cell origin and orientation.
		// These specifications apply only to dendrite connections.
		Number						_minLaminarDist;
		Number						_maxLaminarDist;

		// For bistratified cells (including pyramidals) with dendrites
		// in a basal direction, provide a minimum and maximum
		// laminar distance for connections.
		Number						_minBasalDist;
		Number						_maxBasalDist;

		// The effective area for synapse densities is a function of
		// cylinder length from the compartment and radius of the
		// compartment plus a separation between dendrite and axon.
		Number						_enpassantDist;

		// minAxonDist and maxAxonDist give upper and lower bounds
		// for axon lengths used in making connections. By default
		// a uniform distribution is assumed in generating lengths.
		Number						_minAxonDist;
		Number						_maxAxonDist;

		// Specify synapse densities as a function of location.
		Number						_dendriteSynDensity;
		Number						_somaSynDensity;
		Number						_ISSynDensity;

		// Provide weight and type for synapses
		// An alternate is provided for cases where
		// associated synapse types cannot be grouped.
		Number						_synapseWeight;
		TokenId						_synapseType;

		Number						_synapseWeightAlt;
		TokenId						_synapseTypeAlt;

		// Specify how connections are made. For sequential,
		// make the connections by taking afferent neurons
		// in sequential order, wrapping around as needed.
		// Otherwise, choose afferent neurons at random.
		// _nextToConnect keeps track of next in sequence.
		ConnectionSequence			_connectOrder;
		int							_nextToConnect;

		// Utility fuctions and subclass responsibilities -------------

		// Connect with different parts of the cell, (as appropriate)
		virtual void				connectWithSoma();
		virtual void				connectWithIS();
		virtual void				connectWithDendrites();
		virtual void				connectWithAxon() {} // default is no-op

		// Get an estimated fraction of a target dendrite compartment
		// within a pathway defined by a range of distances from the origin.
		virtual Number				dendriteFractionInRange(
			int						dendNbr,			// number of target dendrite
			Number					minDist,			// minimum laminar distance
			Number					maxDist);			// maximum laminar distance

		// Find the number of synapses to connect. The number is derived
		// from a Poisson distribution with a mean of the value determined 
		// by applying synapse density to the implied dendrite area.
		virtual int					numberToConnect(
			Compartment*			pcomp,				// compartment to connect with
			Number					fraction,			// fraction of comp in layer
			Number					density);			// synapse density

		// Create the necessary synapses by selecting at random from 
		// source cells in the place cell layer.
		virtual void				createSynapses(
			Compartment*			pcomp,				// compartment to connect with
			int						nbrSynapses);		// number of synapses

		// Select an axon length for the current connection.
		// By default this is a uniformly distributed random number
		// between axon min and max distances.
		virtual Number				axonLength(
			Compartment*			pcomp);				// compartment to connect with
	};
};

#endif // #ifndef
