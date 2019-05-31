// Provide classes for simulating a mouse moving in a maze
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: subject_baker_2003.cpp
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		6 March 2007
//
// Description:
//
// To simulate the effects of NMDA knockout in the CA3 portion of the hippocampus,
// it is necessary to generate inputs similar in form to that which might be present
// for a mouse moving in a maze, typically a Morris water maze or its dry equivalent.
//
// Classes implemented here provide a framework for defining a subject animal
// moving in such a maze.


#include "subject_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace BAKER_2003;



// --------------------------------------------------------------------
// MazeSubject class body
// --------------------------------------------------------------------



// Constructors and destructor
MazeSubject::MazeSubject(Maze* m)
{
	ClockSolver*	mySolver = new ClockSolver;
	Model*			myModel = new Model;

	// Set an arbitrary default time step that can be changed later
	SimTime			defaultTimeStep = 25*UOM::msec;

	// Create a default model and clock solver
	model(myModel);
	mySolver->model(myModel);
	mySolver->timeStep(defaultTimeStep);

	// Set initial values
	_maze = NULL;
	_previousMaze = NULL;
	_noveltyRegion = NULL;
	_movementPolicy = NULL;
	_speed = 0;
	_locX = 0;
	_locY = 0;
	_headingX = 1;
	_headingY = 0;

	// Create random number generators
	_networkRandomizer	= new MT19937_UniformRandom();
	_locationRandomizer	= new MT19937_UniformRandom();
	_movementRandomizer	= new MT19937_UniformRandom();
	_pcSpikeRandomizer	= new MT19937_UniformRandom();
	_inSpikeRandomizer	= new MT19937_UniformRandom();
	_synapticRandomizer	= new MT19937_UniformRandom();

	// Hook up the movementRandomizer with the model
	// which is where movement policy looks for its
	// random number stream.
	model()->uniformRandom( movementRandomizer() );

	// Assign the maze if any provided
	if (m!=NULL) {
		maze(m);
	}
}

MazeSubject::~MazeSubject()
{
	// Inform any subscribers that this object is terminated
	changed(terminatedChange);

	// Tell movement policy that this subject is gone
	if (movementPolicy() != NULL) {
		movementPolicy()->subject(NULL);
	}

	// Since solver and model were created for this
	// object, they are now be deleted.
	delete model()->solver();
	delete model();

	// Similarly, random number generators can go
	// back to the bit bucket from which they came.
	delete _networkRandomizer;
	delete _locationRandomizer;
	delete _movementRandomizer;
	delete _pcSpikeRandomizer;
	delete _inSpikeRandomizer;
	delete _synapticRandomizer;
}

// Set the maze in which the subject is located
void MazeSubject::maze(Maze* m) 
{
	// Switch from one maze to another
	_previousMaze = maze();
	_maze = m;

	// Do corresponding updates.
	// Previous maze is available if needed
	// for processing changed/update calls.
	updateLandmarkAngles();
	changed(parameterChange);

	// The old previous maze is no longer relevant.
	// Once updates are processed, the previous and
	// current mazes are the same.
	_previousMaze = maze();
}

// Set the movement policy.
void MazeSubject::movementPolicy(MovementPolicy* policy)
{
	_movementPolicy = policy;
	_movementPolicy->subject(this);
	changed(parameterChange);
}

// Set speed of travel
void MazeSubject::speed(Number v)
{
	_speed = v;
	changed(parameterChange);
}

// Set current location and recompute angles
void MazeSubject::locX(Number x)
{
	_locX = x;
	updateLandmarkAngles();
	changed(stateChange);
}

// Set current location and recompute angles
void MazeSubject::locY(Number y)
{
	_locY = y;
	updateLandmarkAngles();
	changed(stateChange);
}

// Get the current heading angle
Number MazeSubject::headingAngle()
{
	return atan2(_headingY,_headingX);
}

// Set the heading to the supplied angle
void MazeSubject::headingAngle(Number a)
{
	_headingX=cos(a);
	_headingY=sin(a);
	changed(stateChange);
}

// Set location and heading in one function call
void MazeSubject::setLocationAndHeading(
	Number	x,			// location x coord
	Number	y,			// location y coord
	Number	hx,			// heading vector x coord
	Number	hy)			// heading vector y coord
{
	_locX = x;
	_locY = y;
	_headingX = hx;
	_headingY = hy;

	updateLandmarkAngles();
	changed(stateChange);
}

// Set heading to a random direction
void MazeSubject::randomizeHeading()
{
	headingAngle( runif(-Pi,Pi) );
}

// Handle the end of the time step
void MazeSubject::timeStepEnded()
{
	// Take the next step and update landmark angles
	movementPolicy()->takeNextStep();
	updateLandmarkAngles();

	// Inform subscribers that the state has changed
	changed(stateChange);
}

// Compute angles to the various landmarks associated with the maze
void MazeSubject::updateLandmarkAngles()
{
	LandmarkVector&		lm = maze()->landmarks();
	int					nlm = lm.size();
	int					i;

	const double		epsilon = 1e-8;
	double				dx,dy,a;

	// Resize array of angles if needed
	if (_landmarkAngles.size()!=nlm) {
		_landmarkAngles.resize(nlm);
	}

	for (i=0;i<nlm;i++) {
		dx = lm[i]->locX()-locX();
		dy = lm[i]->locY()-locY();
		a = fabs(dx)+fabs(dy)>epsilon ? atan2(dy,dx) : 0;
		_landmarkAngles[i] = a;
	}
}

// Return the value of an internal state.
// See stateLabels for required ordering of variables.
Number MazeSubject::internalStateValue(int n)
{
	// Return the appropriate value
	switch (n) {
	case 0:
		return locX();
	case 1:
		return locY();
	case 2:
		return headingX();
	case 3:
		return headingY();
	default:
		FatalError("(MazeSubject::internalStateValue) Invalid state index.");
		return 0;
	}
}



// --------------------------------------------------------------------
// MovementPolicy class body
// --------------------------------------------------------------------



// Constructors and destructor
MovementPolicy::MovementPolicy(MazeSubject* sub)
{
	_subject = sub;
}

MovementPolicy::~MovementPolicy() {}

// Pass state updates along to subject
void MovementPolicy::setLocationAndHeading(
	Number		x,			// location x coord
	Number		y,			// location y coord
	Number		hx,			// heading vector x coord
	Number		hy)			// heading vector y coord
{
	subject()->setLocationAndHeading(x,y,hx,hy);
}



// --------------------------------------------------------------------
// RandomWalk class body
// --------------------------------------------------------------------



// Constructors and destructor
RandomWalk::RandomWalk(MazeSubject* sub, Number wrate, Number wtau) 
: MovementPolicy(sub)
{
	// Set initial values.

	// Wander rate may not scale exactly with units of measure
	// because of fractal nature of the random walk.
	// These values are only examples.

	wanderRate(wrate);
	wanderTau(wtau);
}

RandomWalk::~RandomWalk() {}

// Take the next random step
void RandomWalk::takeNextStep()
{
	// Provide short names for common variables
	Number		dt = subject()->timeStep();
	Number		stepSize = speed()*dt;

	// Working variables
	Number				x,y,hx,hy,dist;
	Number				hxtmp,hytmp;
	double				r;

	// Set a new trial direction by turning by a random angle from
	// the current heading. Scaling by dt is not entirely accurate
	// but roughly preserves the overall turning rate.
	r=subject()->rnorm() * dt * _currentWanderRate;
	hxtmp = cos(r);
	hytmp = sin(r);

	hx = headingX()*hxtmp-headingY()*hytmp;
	hy = headingX()*hytmp+headingY()*hxtmp;

	// See if a step in the new direction leads outside the maze.
	// If so, try an alternate direction at random.
	dist=maze()->boundaryDistance(locX(),locY(),hx,hy);
	while (dist<=stepSize) {

		// Temporarily decrease randomness in movement
		_currentWanderRate = 0;

		// Set a new heading by rotating a random amount.
		// This is equivalent to picking a random heading but
		// avoids some minor imperfections in the random 
		// number generators.
		double a = 2*Pi*subject()->runif();
		hxtmp=cos(a);
		hytmp=sin(a);
		hx = headingX()*hxtmp-headingY()*hytmp;
		hy = headingX()*hytmp+headingY()*hxtmp;

		dist=maze()->boundaryDistance(locX(),locY(),hx,hy);
	}

	// Use the new heading to get a new position
	x = locX()+speed()*hx*dt;
	y = locY()+speed()*hy*dt;
	setLocationAndHeading(x,y,hx,hy);

	// Update the current wander rate to reflect the passage of time.
	// Precision is not really needed here, but an implicit update is used
	// just in case instability is possible.
	_currentWanderRate += dt*(wanderRate()-_currentWanderRate)/(wanderTau()+dt);
}



// --------------------------------------------------------------------
// WallFollower class body
// --------------------------------------------------------------------



// Constructors and destructor
WallFollower::WallFollower(MazeSubject* sub, Number wd)
: MovementPolicy(sub)
{
	// Set initial defaults.
	wallDistance(wd);
}

WallFollower::~WallFollower() {}

// Take the next random step
void WallFollower::takeNextStep()
{
	using namespace UOM;

	// Provide short names for common variables
	Number			wd = wallDistance();
	Number			dt = subject()->timeStep();
	Number			s = speed()*dt;
	Number			ds = s/100;

	// Get current location and heading
	Number			x = locX();
	Number			y = locY();
	Number			hx = headingX();
	Number			hy = headingY();

	// Remember the old heading
	Number			oldhx = hx;
	Number			oldhy = hy;

	// Working variables
	Number			dx,dy;
	Number			d0,d1;
	double			r,err1,err2;

	// Treat distance from a boundary as a potential
	// function, the differential of which provides 
	// direction of movement to minimize error.
	// The objective function here is (dist-wd)^2;

	err1 = maze()->boundaryDistance(x+ds,y) - wd;
	err2 = maze()->boundaryDistance(x-ds,y) - wd;
	dx = err1*err1-err2*err2;
	err1 = maze()->boundaryDistance(x,y+ds) - wd;
	err2 = maze()->boundaryDistance(x,y-ds) - wd;
	dy = err1*err1-err2*err2;
	r = norm(dx,dy);

	// If it happens that we are at a local minimum
	// then r will be nearly 0. If so, keep the current
	// heading. Otherwise, use the differential direction.
	if (r>0) {
		dx = -dx/r;
		dy = -dy/r;
	}
	else {
		dx = hx;
		dy = hy;
	}

	// If we have reached the desired distance,
	// or more specifically if the next step
	// will cross the desired distance from the wall,
	// we want to move at right angles to the differential
	// to make forward progress along the wall. Other cases
	// are below.

	d0 = maze()->boundaryDistance(x,y);
	d1 = maze()->boundaryDistance(x+s*dx,y+s*dy);

	if (d0>=wd && d1<wd) {

		// Still outside the desired distance but
		// close enough. Move at right angles in
		// a counter clockwise direction.

		hx = -dy;
		hy =  dx;
		d1 = maze()->boundaryDistance(x+s*hx,y+s*hy);
	}
	else if (d1>=wd && d0<wd) {

		// Inside the desired distance.
		// Rotate heading by a -45 deg angle
		// to make some forward progress while
		// improving the distance.

		Number sa = sin(-Pi/4);
		Number ca = cos(-Pi/4);

		hx = dx*ca-dy*sa;
		hy = dy*ca+dx*sa;
		d1 = maze()->boundaryDistance(x+s*hx,y+s*hy);
	}
	else {

		// Otherwise, not yet at the desired
		// distance. Follow the differential.

		hx =  dx;
		hy =  dy;
	}

	// If we have just turned around but the distances
	// do not seem to be improving, try adding a
	// rotation to get out of a prossible loop.
	if (fabs(d1-wd)>=fabs(d0-wd) && hx*oldhx+hy*oldhy<0) {

		Number sa = sin(0.7);
		Number ca = cos(0.7);

		hx = hx*ca-hy*sa;
		hy = hy*ca+hx*sa;
	}

	// Use the new heading to get a new position
	// Only go as far as the boundary will allow.
	// No out of the box experiences allowed here.
	s = minval(s,maze()->boundaryDistance(x,y,hx,hy));
	x = locX()+s*hx;
	y = locY()+s*hy;
	setLocationAndHeading(x,y,hx,hy);
}


