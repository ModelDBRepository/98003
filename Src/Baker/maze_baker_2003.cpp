// Provide classes for simulating a mouse moving in a maze
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: maze_baker_2003.cpp
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
// The specific maze of interest is circular with landmarks arranged around the 
// periphery. Landmarks are specified by their location in a two dimensional
// space. Each landmark is identified through a group of unspecified "features" which
// can be matched with varying degrees of selectivity through a vector dot product.


#include "maze_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace BAKER_2003;


// --------------------------------------------------------------------
// Landmark class body
// --------------------------------------------------------------------


// Constructors and destructor
Landmark::Landmark()
{
	// Set some initial values just to get started
	_locX = 0;
	_locY = 0;
	_category = 0;
}

Landmark::Landmark(Landmark& lm)
{
	// Copy each field in the object using appropriate accessors
	locX(		lm.locX() );
	locY(		lm.locY() );
	category(	lm.category() );
	features(	lm.features() );
}

Landmark::Landmark(
	Number				x,			// Location x coord
	Number				y,			// Location y coord
	int					cat,		// Category
	int					numFeatEnt,	// Number of feature entries
	Number				features[])	// Feature vector as doubles
{
	int					k;

	_locX = (Number) x;
	_locY = (Number) y;
	_category = cat;
	_features.resize(numFeatEnt);

	for (k=0;k<numFeatEnt;k++) {
		_features[k] = features[k];
	}
}

Landmark::~Landmark() {}


// --------------------------------------------------------------------
// Maze class body
// --------------------------------------------------------------------


// Constructors and destructor
Maze::Maze()
{
	// Set initial values
	_originX = 0;
	_originY = 0;
}

Maze::~Maze() 
{
	int k;

	// Delete held landmarks
	for (k=0;k<_landmarks.size(); k++) {
		delete _landmarks[k];
	}
}

// Add a new landmark using a copy of the one provided
void Maze::add(Landmark* lm)
{
	Landmark* newLM = new Landmark(*lm);

	_landmarks.push_back(newLM);
}

// Return true if the indicated point is inside the maze.
bool Maze::isInsideMaze(Number locX, Number locY)
{
	return boundaryDistance(locX,locY)>=0;
}

// Select points at random from inside the maze, for example, as might
// be used to generate place fields.
void Maze::randomPointInMaze(
	Number&				x,				// X coordinate value (output)
	Number&				y,				// Y coordinate value (output)
	UniformRandom*		unif)			// Source of random numbers
{
	Number s = size();

	if (unif==NULL) {
		unif = UniformRandom::defaultGen();
	}

	do {
		x = originX() + unif->next(-s,s);
		y = originY() + unif->next(-s,s);
	} while ( !isInsideMaze(x,y) );
}

// Return a location in local coordinates. These are distance measures
// in different directions as appropriate for this location in the maze.
// By default distances are returned using the local orientation in
// north, south, east, and west corresponding to vectors as rotated by
// the local orientation.
void Maze::localCoordinates(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				northDist,		// Boundary dist along (0,1)
	Number&				southDist,		// Boundary dist along (0,-1)
	Number&				eastDist,		// Boundary dist along (1,0)
	Number&				westDist)		// Boundary dist along (-1,0)
{
	Number				hx,hy;

	// Get local orientation as vector (hx,hy)
	localOrientation(locX,locY,hx,hy);

	// Get relevant distances rotated by the local orientation
	northDist = boundaryDistance(locX,locY,hy,hx);
	southDist = boundaryDistance(locX,locY,-hy,hx);
	eastDist  = boundaryDistance(locX,locY,hx,hy);
	westDist  = boundaryDistance(locX,locY,-hx,-hy);
}

// Provide a heading defined by local conditions in the maze such as
// the direction of the nearest boundary point. The purpose is to provide
// an orientation for a local coordinate system for this part of the maze.
void Maze::localOrientation(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				hx,				// Unit vector x coord (output)
	Number&				hy)				// Unit vector y coord (output)
{
	// Return a constant default
	hx = 1;
	hy = 0;
}



// --------------------------------------------------------------------
// CircularMaze class body
// --------------------------------------------------------------------



CircularMaze::CircularMaze(Number r, Number origX, Number origY)
{
	radius(r);
	originX(origX);
	originY(origY);
}

CircularMaze::~CircularMaze() {}

// Provide a heading defined by local conditions in the maze.
// Value returned is orientation towards nearest boundary point.
void CircularMaze::localOrientation(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				hx,				// Unit vector x coord (output)
	Number&				hy)				// Unit vector y coord (output)
{
	Number				len = norm(locX,locY);
	Number				epsilon = 1e-8f;

	// Return vector pointing directly to the boundary, which is
	// an easy thing to do in a circle unless you are right at the
	// center, in which case a vector of (1.0) is returned.

	if (len>epsilon) {
		hx = (locX-originX() )/len;
		hy = (locY-originY() )/len;
	}
	else {
		hx = 1;
		hy = 0;
	}
}

// Return a location in local coordinates. These are distance measures
// in different directions as appropriate for this location in the maze.
// By default, distances are returned using the local orientation in
// north, south, east, and west corresponding to vectors as rotated by
// the local orientation.
void CircularMaze::localCoordinates(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				northDist,		// Boundary dist along (0,1)
	Number&				southDist,		// Boundary dist along (0,-1)
	Number&				eastDist,		// Boundary dist along (1,0)
	Number&				westDist)		// Boundary dist along (-1,0)
{
	Number				r = size();
	Number				len = norm(locX-originX(),locY-originY());
	Number				sideLenSq = r*r-len*len;

	// Distances can be computed directly for a circle. Note that the
	// east distance goes negative when (locX,locY) is outside the circle.
	eastDist = r-len;
	westDist = r+len;

	// Make sure (locX,locY) is inside the circle. Otherwise, return a
	// negative value for north and south distances.
	if (sideLenSq>=0) {
		northDist = southDist = sqrt(sideLenSq);
	}
	else {
		// Set distances to a negative value if outside the square.
		// This value is not especially meaningful, but goes to zero as 
		// (locX,locY) approaches the boundary from the outside.
		northDist = southDist = -sqrt(-sideLenSq);
	}
}

// Distance to the maze wall from a given location
Number CircularMaze::boundaryDistance(Number locX, Number locY)
{
	return size()-norm(locX-originX(),locY-originY());
}

// Compute the distance to the maze boundary in a given direction
Number CircularMaze::boundaryDistance(
	Number				locX,		// Current location X coord
	Number				locY,		// Current location Y coord
	Number				hx,			// Unit heading vector x coord
	Number				hy)			// Unit heading vector y coord
{	
	Number				x0 = locX-originX();
	Number				y0 = locY-originY();
	Number				r  = size();

	// Determining the distance is solving for d in 
	// the following simultaneous equations:
	//
	// d>=0
	// x=x0+d*hx 
	// y=y0+d*hy
	// x*x+y*y=r*r 
	//
	// where d is the desired distance. Of course, this is easily done by 
	// solving a quadratic equation resulting from substituting for x,y.
	// We know the solution is real though negative when outside the circle.

	Number	a = hx*hx+hy*hy; // should = 1, but be safe anyway
	Number	b = 2*(x0*hx+y0*hy);
	Number	c = x0*x0+y0*y0-r*r;
	Number	d = (-b+sqrt(b*b-4*a*c))/(2*a);

	return d;
}


// --------------------------------------------------------------------
// RectangularMaze class body
// --------------------------------------------------------------------



RectangularMaze::RectangularMaze(
	Number sizex, Number sizey, 
	Number origx, Number origy)
{
	sizeX(sizex);
	sizeY(sizey);
	originX(origx);
	originY(origy);
}

RectangularMaze::~RectangularMaze() {}

// Return a location in local coordinates. These are distance measures
// in different directions as appropriate for this location in the maze.
// By default, distances are returned using the local orientation in
// north, south, east, and west corresponding to vectors as rotated by
// the local orientation.
void RectangularMaze::localCoordinates(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				northDist,		// Boundary dist along (0,1)
	Number&				southDist,		// Boundary dist along (0,-1)
	Number&				eastDist,		// Boundary dist along (1,0)
	Number&				westDist)		// Boundary dist along (-1,0)
{
	northDist = originY()+sizeY()-locY;
	southDist = locY - originY()+sizeY();

	eastDist  = originX()+sizeX()-locX;
	westDist  = locX - originX()+sizeX();
}

// Provide a heading defined by local conditions in the maze.
// Value returned is orientation towards nearest boundary point.
void RectangularMaze::localOrientation(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number&				hx,				// Unit vector x coord (output)
	Number&				hy)				// Unit vector y coord (output)
{
	Number ndist,sdist,edist,wdist;
	Number xdist,ydist;

	localCoordinates(locX,locY,ndist,sdist,edist,wdist);

	xdist = minval(edist,wdist);
	ydist = minval(ndist,sdist);

	if (xdist<=ydist) {
		hx = edist<=wdist ? 1 : -1;
		hy =0;
	}
	else {
		hx = 0;
		hy = ndist<=sdist ? 1 : -1;
	}
}

// Provide minimum distance to the maze wall from a given location.
// Return a negative value if the current location is outside the maze.
Number RectangularMaze::boundaryDistance(
	Number				locX,			// Current location X coord
	Number				locY)			// Current location Y coord
{
	Number ndist,sdist,edist,wdist;
	Number xdist,ydist;

	localCoordinates(locX,locY,ndist,sdist,edist,wdist);
	xdist = minval(edist,wdist);
	ydist = minval(ndist,sdist);

	return minval(xdist,ydist);
}

// Provide the distance from maze wall along a unit vector.
// Return a negative value if the current location is outside the maze.
Number RectangularMaze::boundaryDistance(
	Number				locX,			// Current location X coord
	Number				locY,			// Current location Y coord
	Number				hx,				// Unit vector x coord
	Number				hy)				// Unit vector y coord
{
	Number ndist,sdist,edist,wdist;
	Number sx,sy;

	localCoordinates(locX,locY,ndist,sdist,edist,wdist);

	// S is the distance along (hx,hy) until the boundary
	// is hit. sx is the x constraint and sy the y constaint.

	if (hx>0)		sx= edist/hx;
	else if (hx<0)	sx=-wdist/hx;
	else			sx=size();

	if (hy>0)		sy= ndist/hy;
	else if (hy<0)	sy=-sdist/hy;
	else			sy=size();

	return minval(sx,sy);
}



// --------------------------------------------------------------------
// RectangularRegion class body
// --------------------------------------------------------------------



// Constructor and destructor
RectangularRegion::RectangularRegion(
	Number				xmin,
	Number				xmax,
	Number				ymin,
	Number				ymax)
{
	_xmin = xmin;
	_xmax = xmax;
	_ymin = ymin;
	_ymax = ymax;
}

RectangularRegion::~RectangularRegion() {}

// Get the bounding box of the region
void RectangularRegion::getBoundingBox(
	Number&				xmin, 
	Number&				xmax,
	Number&				ymin,
	Number&				ymax)
{
	xmin = _xmin;
	xmax = _xmax;
	ymin = _ymin;
	ymax = _ymax;
}

// Set the bounding box of the region
void RectangularRegion::setBoundingBox(
	Number				xmin, 
	Number				xmax,
	Number				ymin,
	Number				ymax)
{
	_xmin = xmin;
	_xmax = xmax;
	_ymin = ymin;
	_ymax = ymax;
}

// Return true if the point supplied is in the region
bool RectangularRegion::containsPoint(
	Number				x,
	Number				y)
{
	return 
		x >= _xmin &&
		x <= _xmax &&
		y >= _ymin &&
		y <= _ymax;
}



// --------------------------------------------------------------------
// CircularRegion class body
// --------------------------------------------------------------------




// Constructor and destructor
CircularRegion::CircularRegion(Number xctr, Number yctr, Number r)
{
	_centerX = xctr;
	_centerY = yctr;
	_radius = r;
}

CircularRegion::~CircularRegion() {}

// Get the bounding box of the region
void CircularRegion::getBoundingBox(
	Number&				xmin, 
	Number&				xmax,
	Number&				ymin,
	Number&				ymax)
{
	xmin = centerX()-radius();
	xmax = centerX()+radius();
	ymin = centerY()-radius();
	ymax = centerY()+radius();
}

// Return true if the point supplied is in the region
bool CircularRegion::containsPoint(
	Number				x,
	Number				y)
{
	Number				dx = x-centerX();
	Number				dy = y-centerY();

	return (dx*dx+dy*dy)<=radius()*radius();
}






