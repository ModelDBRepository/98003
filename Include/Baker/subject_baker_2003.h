// Provide classes for simulating a mouse moving in a maze
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: subject_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header provides a framework for defining an anmimal moving in such an
// environment such as a Morris Water Maze. The definition here is abstract
// and concrete implementation involves extension of the MazeSubject class.
//
// A basic movement policy is provided in which the animal moves about the
// maze at random without attempting to seek a goal or explicitly explore
// unfamiliar regions of the maze.
//
// References:
//
// Dudek G, Jenkin M (2000) Computational principles of mobile robotics.
// New York: Cambridge University Press.
//
// Nakazawa K, Quick MC, Chitwood RA, Watanabe M, Yeckel MF, Sun LD, 
// Kata A, Carr CA, Johnston D, Wilson MA, Tonegawa S (2002) Requirement
// for hippocampal CA3 NMDA receptors in associative memory recall.
// Science 297: 211-218.


// --------------------------------------------------------------------
// Only include the definitions in this header once
// --------------------------------------------------------------------

#ifndef __SUBJECT_BAKER_2003_H_
#define __SUBJECT_BAKER_2003_H_


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
#include "maze_baker_2003.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace BAKER_2003 {

	// ----------------------------------------------------------------
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ----------------------------------------------------------------

	class MazeSubject;

	class MovementPolicy;
		class RandomWalk;
		class WallFollower;


	// ----------------------------------------------------------------
	// CLASS:	MazeSubject
	// EXTENDS:	ModelComponent
	// DESC:	Abstract class for defining subjects for maze
	//			experiments. .
	// RESP:
	//		1.	Store relationship with maze.
	//		2.	Know current position and heading.
	//		3.	Allocate private model and solver if none provided
	//		4.	Compute angular relationship between current position
	//			and all visible landmarks.
	//		5.	Know current novelty region	
	//
	// NOTES:	For now all landmarks are assumed to be visible.
	//			Because a clock solver is used, subclasses must
	//			provide motion in terms of discrete change in location
	//			with each time step.
	//
	//			Subclass must establish any behaviors associated with
	//			being within the novelty region, including possibly
	//			passing the region to affected place cell layers.
	// ----------------------------------------------------------------

	class MazeSubject : public ModelComponent {

	public:

		// Constructors and destructor
		MazeSubject(Maze* m=NULL);
		virtual ~MazeSubject();

		// Accessors
		inline  Maze*			maze() { return _maze; }
		virtual void			maze(Maze* m);

		inline  SpatialRegion*	noveltyRegion() { return _noveltyRegion; }
		virtual void			noveltyRegion(SpatialRegion* nr) { _noveltyRegion=nr; }

		inline  MovementPolicy* movementPolicy() { return _movementPolicy; }
		virtual void			movementPolicy(MovementPolicy* policy);

		inline  Number			speed() { return _speed; }
		virtual void			speed(Number v);

		inline	Number			locX() { return _locX; }
		virtual void			locX(Number x);

		inline  Number			locY() { return _locY; }
		virtual void			locY(Number y);

		inline  Number			headingX() { return _headingX; }
		virtual void			headingX(Number x) { _headingX = x; }

		inline  Number			headingY() { return _headingY; }
		virtual void			headingY(Number y) { _headingY = y; }

		virtual Number			headingAngle();
		virtual void			headingAngle(Number angle);

		virtual ModelComponentVector* subscribers() { return &_subscribers; }

		// Provide access to the previous maze, but only while switching over.
		inline  Maze*			previousMaze() { return _previousMaze; }

		// Interface for updating location and heading in one function call.
		// This reduces the overhead of multiple change notifications.
		virtual void			setLocationAndHeading(
			Number				x,			// location x coord
			Number				y,			// location y coord
			Number				hx,			// heading vector x coord
			Number				hy);		// heading vector y coord

		// Set heading to a uniformly distributed random direction
		virtual void			randomizeHeading();

		// Access sensory information available from the current location.
		// For now all landmarks are assumed to be fully visible at all locations.
		virtual LandmarkVector&	landmarks() { return maze()->landmarks(); }
		virtual NumberArray		landmarkAngles() { return _landmarkAngles; }

		// Interface for various sources of random numbers ------------

		inline  UniformRandom*	networkRandomizer()		{ return _networkRandomizer; }
		inline  UniformRandom*	locationRandomizer()	{ return _locationRandomizer; }
		inline  UniformRandom*	movementRandomizer()	{ return _movementRandomizer; }
		inline  UniformRandom*	pcSpikeRandomizer()		{ return _pcSpikeRandomizer; }
		inline  UniformRandom*	inSpikeRandomizer()		{ return _inSpikeRandomizer; }
		inline  UniformRandom*	synapticRandomizer()	{ return _synapticRandomizer; }

		// Interface with Clock Solver --------------------------------

		inline  SimTime			timeStep() { return model()->solver()->timeStep(); }
		virtual void			timeStep(SimTime h) { model()->solver()->timeStep(h); }

		virtual void			timeStepEnded();

		// Accessors for probing state values -------------------------

		virtual int				numStateVar() { return 0; }
		virtual int				numInternalStateVar() { return 4; }
		virtual Number			internalStateValue(int n);

		virtual const char*		componentName() { return "MazeSubject"; }
		virtual const char**	stateLabels() {
			static const char* sl[]={"locX","locY","headingX","heandingY"}; 
			return sl; }

		virtual Number*			unitsOfMeasure() {
			static Number um[] = {UOM::cm, UOM::cm, 1.0, 1.0 }; return &(um[0]); }

	protected:

		ModelComponentVector	_subscribers;	// subscribers to be informed of state changes
		MovementPolicy*			_movementPolicy; // Policy object implementing movements

		// Different random number streams, each of which can be
		// separately controlled, including setting of seeds for
		// reproducible results. The subject instance owns these
		// object and is responsible for creating and deleting them.
		UniformRandom*			_networkRandomizer;		// for establishing the network
		UniformRandom*			_locationRandomizer;	// for setting place cell locations
		UniformRandom*			_movementRandomizer;	// for for generating movements
		UniformRandom*			_pcSpikeRandomizer;		// for place cell spike trains
		UniformRandom*			_inSpikeRandomizer;		// for interneuron spike trains
		UniformRandom*			_synapticRandomizer;	// for random events in synapses

		// Parameters
		Maze*					_maze;			// current maze environment or NULL
		Maze*					_previousMaze;	// old maze before changing to a new one
		SpatialRegion*			_noveltyRegion;	// novelty region or NULL
		Number					_speed;			// walk/swim speed

		// State variables
		Number					_locX;			// current location X coord
		Number					_locY;			// current location Y coord
		Number					_headingX;		// current heading vector X coord
		Number					_headingY;		// current heading vector Y coord

		// Derived information
		NumberArray				_landmarkAngles; // angle to each landmark from current location

		// Compute angles to landmarks from current location.
		virtual	void			updateLandmarkAngles();
	};

	// ----------------------------------------------------------------
	// CLASS:	MovementPolicy
	// EXTENDS:	none
	// DESC:	Abstract class for maze subject movement policies.
	// RESP:
	//		1.	Know maze subject
	//		2.	Provide interface for movement
	//
	// ----------------------------------------------------------------

	class MovementPolicy {

	public:

		// Constructors and destructor
		MovementPolicy(MazeSubject* sub=NULL);
		virtual ~MovementPolicy();

		// Accessors
		inline  MazeSubject*	subject() { return _subject; }
		virtual void			subject(MazeSubject* sub) { _subject = sub; }

		inline  Maze*			maze() { return subject()->maze(); }

		inline  Number			locX() { return subject()->locX(); }
		inline  Number			locY() { return subject()->locY(); }

		inline  Number			headingX() { return subject()->headingX(); }
		inline  Number			headingY() { return subject()->headingY(); }

		inline  Number			speed() { return subject()->speed(); }

		// Pass updates of state on to subject
		virtual void			setLocationAndHeading(
			Number				x,				// location x coord
			Number				y,				// location y coord
			Number				hx,				// heading vector x coord
			Number				hy);			// heading vector y coord

		// Update the maze subject state to take the next movement
		virtual void			takeNextStep()=0;	// subclass responsibility

	protected:
		MazeSubject*			_subject;		// subject using this policy
	};

	// ----------------------------------------------------------------
	// CLASS:	RandomWalk
	// EXTENDS:	MovementPolicy
	// DESC:	Defines movements that wanders at random in a
	//			Maze. Path is a random walk combining monentum,
	//			noisey changes in direction, and random reset
	//			of direction upon encountering a maze wall.
	// RESP:
	//		1.	Know parameters for random walk.
	//		2.	Determine next change in location for a time step.
	//
	// NOTES:	Parameter defaults are intended to create a space-
	//			filling random walk with somewhat mouse-like aspects,
	//			there is no attempt to statistically model specific
	//			observed mouse data. To avoid excessive time spent
	//			bouncing off walks, randomness in the walk is
	//			temporarily set to zero and recovers with a time
	//			constant specified through wanderTau.
	// ----------------------------------------------------------------

	class RandomWalk : public MovementPolicy {

	public:

		// Constructors and destructor
		RandomWalk(MazeSubject* sub=NULL, 
			Number wrate=1/UOM::sec,			// wander rate value
			Number wtau=1*UOM::sec);			// wander tau value
		virtual ~RandomWalk();

		// Accessors
		inline  Number			wanderRate() { return _wanderRate; }
		virtual void			wanderRate(Number wr) { _wanderRate=_currentWanderRate=wr; }

		inline  Number			wanderTau() { return _wanderTau; }
		virtual void			wanderTau(Number tm) { _wanderTau = tm; }

		// Take the next step
		virtual void			takeNextStep();

	protected:
		Number					_wanderRate;		// Initial rate of random direction changes
		Number					_wanderTau;			// Time constant for recovering randomness
		Number					_currentWanderRate;	// Wander rate as of the present
	};

	// ----------------------------------------------------------------
	// CLASS:	WallFollower
	// EXTENDS:	MovementPolicy
	// DESC:	Defines movements that wanders while staying near
	//			a boundary. Movements are somewhat irregular but
	//			at each step there is an objective of staying
	//			a given distance from a boundary.
	// RESP:
	//		1.	Know desired boundary distance.
	//		2.	Determine next change in location for a time step.
	// ----------------------------------------------------------------

	class WallFollower : public MovementPolicy {

	public:

		// Constructors and destructor
		WallFollower(MazeSubject* sub=NULL, Number wd=2*UOM::cm);
		virtual ~WallFollower();

		// Accessors
		inline  Number			wallDistance() { return _wallDistance; }
		virtual void			wallDistance(Number d) { _wallDistance = d; }

		// Take the next step
		virtual void			takeNextStep();

	protected:
		Number					_wallDistance;
	};
};

#endif // #ifndef
