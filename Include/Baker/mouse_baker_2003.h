// Provide classes for simulating a mouse moving in a maze
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
// File: mouse_baker_2003.h
//
// Description:
//
// To simulate the effects of NMDA knockout in the CA3 portion of the hippocampus,
// it is necessary to generate inputs similar in form to that which might be present
// for a mouse moving in a maze, typically a Morris water maze or its dry equivalent.
//
// The maze subject is nominally a mouse with motion parameters similar to those
// reported by Nakazawa. Goal seeking behavior is not included in the simulation,
// Only random walk and wall following movements are supported.
//
// References:
//
// Amaral DG, Ishizuka N, Claiborne B (1990) Neurons, numbers and 
// the hippocampal network. Progress in Brain Research 83: 1-11.
//
// O'Keefe J, Dostrovsky J (1971) The hippocampus as a spatial map:
// preliminary evidence from unit activity in the freely moving rat.
// Brain Research 34: 171-175.
//
// Nakazawa K, Quick MC, Chitwood RA, Watanabe M, Yeckel MF, Sun LD, 
// Kata A, Carr CA, Johnston D, Wilson MA, Tonegawa S (2002) 
// Requirement for hippocampal CA3 NMDA receptors in associative memory
// recall. Science 297: 211-218.



// --------------------------------------------------------------------
// Only include the definitions in this header once
// --------------------------------------------------------------------

#ifndef __MOUSE_BAKER_2003_H_
#define __MOUSE_BAKER_2003_H_


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


#include "hippoform_baker_2003.h"
#include "neuron_baker_2003.h"

using namespace std;
using namespace BNSF;

// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace BAKER_2003 {


	// ----------------------------------------------------------------
	// CLASS:	Mouse
	// EXTENDS:	MazeSubject
	// DESC:	Define class to model movements and brain activity
	//			of a mouse in a maze. The objective here is modeling
	//			response of a realistic CA3 pyramidal cell using
	//			phomenologically modeled inputs.
	// RESP:
	//		1.	Set up random walk movement policy
	//		2.	Define EC, DG, and CA3 place cell layers
	//		3.	Define CA3 interneurons
	//		4.	Connect with a CA3 pyramidal cell model
	//
	// NOTE:	This model is basically the way a hippocampus specialist 
	//			would think of a mouse, i.e. a hippocampus with some 
	//			trivial behavior attached. No attempt at even a low 
	//			fidelity simulation of mouse behavior is intended.
	// ----------------------------------------------------------------

	class Mouse : public MazeSubject {

	public:

		// typedefs and enums
		enum MovementType {wallFollower, randomWalk};		

		// Constructors and destructor
		Mouse(
			Maze*				m=NULL,					// maze to use
			MovementType		moveType=randomWalk,	// type of movements generated
			double*				seedValues=NULL,		// array of 6 values in (0,1)
			bool				buildNet=true,			// build network connections or not
			bool				doInit=true);			// do initializations or not
		virtual ~Mouse();

		// Accessor for default random number seed values.
		// The intent here is to allow access so that an alternaive
		// array of seed values can be built to supply to a constructor.
		static double*			defaultSeeds(int n) { return _defaultSeeds[n]; }

		// Accessors for spike generators
		inline  ECLayer*			ECPlaceCells() { return _ECPlaceCells; }
		inline  DGLayer*			DGPlaceCells() { return _DGPlaceCells; }
		inline  CA3Layer*			CA3PlaceCells() { return _CA3PlaceCells; }

		inline  CA3AxoAxonicLayer*	CA3AxoAxonicCells() { return _CA3AxoAxonicCells; }
		inline  CA3BasketLayer*		CA3BasketCells() { return _CA3BasketCells; }
		inline  CA3BistratifiedLayer* CA3BistratifiedCells() { return _CA3BistratifiedCells; }
		inline  CA3OLMLayer*		CA3OLMCells() { return _CA3OLMCells; }

		// Access the maze object. If place cell layers have already
		// been established, place field centers will be set also.
		inline  Maze*			maze() { return _maze; }
		virtual void			maze(Maze* m);

		// Access the novel region. This region is also provided to
		// place cell layers as needed. NULL means no such region.
		inline  SpatialRegion*	noveltyRegion() { return _noveltyRegion; }
		virtual void			noveltyRegion(SpatialRegion* novelArea);

		// Access the target cell. Default is L56a in CA3b.
		// This can be changed prior to initialization.
		// In such a case, the cell object is not deleted
		// when the mouse is destroyed.
		inline  CA3PyramidalCell* targetCell() { return _targetCell; }
		virtual void			targetCell(CA3PyramidalCell* cell);

		// Set randomizer seeds in one function call.
		// Seed values must reference an array of 5 values.
		virtual void			setSeedValues(double* seedValues);

		// Indicate whether synpases with target cell are to be created.
		// This is a debug option that allows all input spike trains
		// to be generated without the overhead of processing them
		// against target cell synapses..
		inline	bool			buildNetwork() { return _buildNetwork; }
		virtual void			buildNetwork(bool x) { _buildNetwork = x; }

		// Print synapses statistics to cerr for debug
		virtual	void			printSynapseCounts();

		// When this is added to a controller, add neurons as well 
		virtual void			addToController(Controller* cont);

		// Initialize this instance once everything is set
		virtual void			initialize();		

		// Functions to build components of the network
		// Subclasses can override as needed for variations
		virtual void			allocateTargetCell();
		virtual void			allocateCellLayers();
		virtual void			setLayerActivities();
		virtual void			connectNetwork();

	protected:

		// Static constants
		static double			_defaultSeeds[8][6];

		// Control parameters
		bool					_buildNetwork;
		bool					_ownsTargetCell;

		// Target cell model
		CA3PyramidalCell*		_targetCell;

		// Place cell spike generators
		ECLayer*				_ECPlaceCells;
		DGLayer*				_DGPlaceCells;
		CA3Layer*				_CA3PlaceCells;

		// Interneuron spike generators
		CA3AxoAxonicLayer*		_CA3AxoAxonicCells;
		CA3BasketLayer*			_CA3BasketCells;
		CA3BistratifiedLayer*	_CA3BistratifiedCells;
		CA3OLMLayer*			_CA3OLMCells;

		// Connection policies for above
		ECPerforantPath*		_PP;
		DGMossyFibers*			_MF;
		CA3AssocCollaterals*	_AC;

		CA3AxonicInhibition*	_INaxon;
		CA3SomaticInhibition*	_INsoma;
		CA3ProximalInhibition*	_INproximal;
		CA3DistalInhibition*	_INdistal;
	};
};

#endif // #ifndef
