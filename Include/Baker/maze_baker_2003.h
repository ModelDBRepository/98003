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
// File: maze_baker_2003.h
//
// Description:
//
// This header declares classes used to define an environment such a Morris
// Water Maze, though not limited to that specific form of maze. A basic
// definition for associated landmarks, such as might be used for navigation,
// is also provided.
//
// A specific maze of interest is circular with landmarks arranged around the 
// periphery. Landmarks are specified by their location in a two dimensional
// space. Each landmark is identified through a group of unspecified "features" which
// can be matched with varying degrees of selectivity through a vector dot product.
//
// References:
//
// O'Keefe J, Burgess N (1996) Geometric determinants of the place fields
// of hippocampal neurons. Nature 381: 425-428.
//
// Hartley T, Burgess N, Lever C, Cacucci F, O'Keefe J (2000) Modeling
// place fields in terms of the cortical inputs to the hippocampus.
// Hippocampus 10: 369-379.



namespace BAKER_2003 {

	class Landmark;
	class Maze;
		class CircularMaze;
}

// --------------------------------------------------------------------
// Only include the definitions in this header once
// --------------------------------------------------------------------

#ifndef __MAZE_BAKER_2003_H_
#define __MAZE_BAKER_2003_H_


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

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace BAKER_2003 {

	// ----------------------------------------------------------------
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ----------------------------------------------------------------

	class Landmark;

	class Maze;
		class CircularMaze;
		class RectangularMaze;

	class SpatialRegion;
		class RectangularRegion;
		class CircularRegion;


	// ----------------------------------------------------------------
	// Vector and iterator typedefs for classes used here.
	// Note that by convention pointers are stored in
	// vectors so this info is not in the typedef name.
	// ----------------------------------------------------------------

	typedef vector<Landmark*>						LandmarkVector;
	typedef LandmarkVector::iterator				LandmarkVectorIt;

	// ---------------------------------------------------
	// CLASS:	Landmark
	// EXTENDS:	none
	// DESC:	Store information relating to a landmark
	//			including location and visible features.
	// RESP:
	//		1.	Store landmark location (world coordinates)
	//		2.	Store category
	//		3.	Store feature scores.
	//
	// NOTES:	Category is an arbitrary categorization that
	//			can be supplied for the landmark. Similarly,
	//			the content of features is not defined except
	//			by how the values are used outside of this class.
	// ---------------------------------------------------

	class Landmark {
	public:

		// Constructors and destructor
		Landmark();							// Empty constructor
		Landmark(Landmark& lm);				// Copy constructor

		Landmark(							// Alternate constructor
			Number				x,			// Location x coord
			Number				y,			// Location y coord
			int					cat,		// Category
			int					numFeatEnt,	// Number of feature entries
			Number				features[]);	// Feature vector

		virtual ~Landmark();

		// Accessors
		inline  Number			locX() { return _locX; }
		virtual void			locX(Number x) { _locX = x; }

		inline  Number			locY() { return _locY; }
		virtual void			locY(Number y) { _locY = y; }

		inline  int				category() { return _category; }
		virtual void			category(int cat) { _category = cat; }

		inline  NumberArray		features() { return _features; }
		virtual void			features(NumberArray fv) { _features = fv; }

	protected:

		Number					_locX;		// location x coord
		Number					_locY;		// location y coord
		int						_category;	// category designation
		NumberArray				_features;	// feature values
	};
			
	// ---------------------------------------------------
	// CLASS:	Maze
	// EXTENDS:	none
	// DESC:	Abstract class for representing mazes.
	// RESP:
	//		1.	Store maze overall size and origin location.
	//		2.	Store associated distal landmarks.
	//
	// NOTES:	Subclasses must provide functions for locating
	//			walls. Absolute angles are measured in a counterclockwise 
	//			direction using a 0 degree orientation along the
	//			X axis.
	//
	//			Landmark objects are copied when added to the maze.
	//			Hence the contents of a given landmark can only be
	//			changed by accessing the landmarks vector directly.
	//
	//			A local coordinate system is imposed to allow some
	//			generality of maze geometry. Local coordinates do not
	//			necessarily define a unique location in the world
	//			coordinate system, but should correspond with the
	//			type of measurements an animal might make to estimate
	//			position within the maze, at least locally.
	// ---------------------------------------------------

	class Maze {

	public:

		// Constructor and destructor
		Maze();
		virtual ~Maze();

		// Accessors
		inline  Number			originX() { return _originX; }
		virtual void			originX(Number x) { _originX = x; }

		inline  Number			originY() { return _originY; }
		virtual void			originY(Number y) { _originY = y; }

		inline  LandmarkVector&	landmarks() { return _landmarks; }

		// Add a landmark to the existing set copying from the landmark provided
		virtual void			add(Landmark* lm);

		// Utility functions (may be overridden in subclasses) --------

		// Return true if the indicated location is inside the maze.
		// Subclasses may want to override a more efficient implementation.
		virtual bool			isInsideMaze(Number locX, Number locY);

		// Select a point at random from inside the maze, as might 
		// be used to generate place fields. The default method 
		// assumes that the extent of the maze is plus/minus 
		// size() in both X and Y coordinates. Output is placed in 
		// x and y. If the source of random numbers is null, 
		// the default UniformRandom generator is used.
		virtual void			randomPointInMaze(
			Number&				x,				// X coordinate value (output)
			Number&				y,				// Y coordinate value (output)
			UniformRandom*		unif=NULL);		// Source of random numbers

		// Return a location in local coordinates. These are distance measures
		// in different directions as appropriate for this location in the maze.
		// By default, distances are returned using the directios north, south, 
		// east, and west corresponding to vectors as rotated by the local orientation. 
		virtual void			localCoordinates(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number&				northDist,		// Boundary dist along (0,1)
			Number&				southDist,		// Boundary dist along (0,-1)
			Number&				eastDist,		// Boundary dist along (1,0)
			Number&				westDist);		// Boundary dist along (-1,0)

		// Subclass responsibilities ----------------------------------

		// Return the overall size of the maze
		virtual Number			size() = 0;

		// Return the approximate area of the maze
		virtual Number			area() { return 4*size()*size(); } // default only

		// Provide a heading defined by local conditions in the maze such as
		// the direction of the nearest boundary point. The purpose is to provide
		// an orientation for a local coordinate system for this part of the maze.
		// Default is a constant vector (1,0). Subclass should override as needed.
		virtual void			localOrientation(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number&				hx,				// Unit vector x coord (output)
			Number&				hy);			// Unit vector y coord (output)

		// Provide minimum distance to the maze wall from a given location.
		// Return a negative value if the current location is outside the maze.
		virtual Number			boundaryDistance(
			Number				locX,			// Current location X coord
			Number				locY) = 0;		// Current location Y coord

		// Provide the distance from maze wall along a unit vector.
		// Return a negative value if the current location is outside the maze.
		virtual Number			boundaryDistance(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number				hx,				// Unit vector x coord
			Number				hy) = 0;		// Unit vector y coord

	protected:
		Number					_originX;		// Maze origin X coord
		Number					_originY;		// Maze origin Y coord
		LandmarkVector			_landmarks;		// Distal landmarks
	};

	// ---------------------------------------------------
	// CLASS:	CircularMaze
	// EXTENDS:	Maze
	// DESC:	Define a circular maze (ala Water Maze).
	// RESP:
	//		1.	Compute distance to walls in a circle.
	//
	// NOTE:	Origin of maze is the center. If location is
	//			outside the circle, distance to wall is negative.
	// ---------------------------------------------------

	class CircularMaze : public Maze {

	public:

		// Constructor and destructor
		CircularMaze(
			Number radius=25*UOM::cm,			// radius of the circle
			Number originX=0,					// origin X coord
			Number originY=0);					// origin Y coord
		virtual ~CircularMaze();

		// Accessors
		inline  Number			radius() { return _radius; }
		virtual void			radius(Number r) { _radius = r; }
		virtual Number			size() { return _radius; }

		// Return the  area of the maze
		virtual Number			area() { return Pi*radius()*radius(); }

		// Provide a heading defined by local conditions in the maze.
		// Value returned is orientation towards nearest boundary point.
		virtual void			localOrientation(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number&				hx,				// Unit vector x coord (output)
			Number&				hy);			// Unit vector y coord (output)

		// Return a location in local coordinates. These are distance measures
		// in different directions as appropriate for this location in the maze.
		// By default, distances are returned using the local orientation in
		// north, south, east, and west corresponding to vectors as rotated by
		// the local orientation.
		virtual void			localCoordinates(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number&				northDist,		// Boundary dist along (0,1)
			Number&				southDist,		// Boundary dist along (0,-1)
			Number&				eastDist,		// Boundary dist along (1,0)
			Number&				westDist);		// Boundary dist along (-1,0)

		// Provide minimum distance to the maze wall from a given location.
		// Return a negative value if the current location is outside the maze.
		virtual Number			boundaryDistance(
			Number				locX,			// Current location X coord
			Number				locY);			// Current location Y coord

		// Provide the distance from maze wall along a unit vector.
		// Return a negative value if the current location is outside the maze.
		virtual Number			boundaryDistance(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number				hx,				// Unit vector x coord
			Number				hy);			// Unit vector y coord

	protected:
		Number					_radius;
	};

	// ---------------------------------------------------
	// CLASS:	RectangularMaze
	// EXTENDS:	Maze
	// DESC:	Define a rectangular maze.
	// RESP:
	//		1.	Compute distance to walls in a rectangle.
	//
	// ---------------------------------------------------

	class RectangularMaze : public Maze {

	public:

		// Constructor and destructor
		RectangularMaze(
			Number sizeX=25*UOM::cm,			// +/- size in X direction
			Number sizeY=25*UOM::cm,			// +/- size in Y direction
			Number originX=0,					// origin X coord
			Number originY=0);					// origin Y coord

		virtual ~RectangularMaze();

		// Accessors
		inline  Number			sizeX() { return _sizeX; }
		virtual void			sizeX(Number s) { _sizeX = s; }

		inline  Number			sizeY() { return _sizeY; }
		virtual void			sizeY(Number s) { _sizeY = s; }

		// Get the overall size of the maze
		virtual Number			size() { return maxval(sizeX(), sizeY()); }

		// Return the area of the maze
		virtual Number			area() { return 4*sizeX()*sizeY(); }

		// Provide a heading defined by local conditions in the maze.
		// Value returned is orientation towards nearest boundary point.
		virtual void			localOrientation(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number&				hx,				// Unit vector x coord (output)
			Number&				hy);			// Unit vector y coord (output)

		// Return a location in local coordinates. These are distance measures
		// in different directions as appropriate for this location in the maze.
		// By default, distances are returned using the local orientation in
		// north, south, east, and west corresponding to vectors as rotated by
		// the local orientation.
		virtual void			localCoordinates(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number&				northDist,		// Boundary dist along (0,1)
			Number&				southDist,		// Boundary dist along (0,-1)
			Number&				eastDist,		// Boundary dist along (1,0)
			Number&				westDist);		// Boundary dist along (-1,0)

		// Provide minimum distance to the maze wall from a given location.
		// Return a negative value if the current location is outside the maze.
		virtual Number			boundaryDistance(
			Number				locX,			// Current location X coord
			Number				locY);			// Current location Y coord

		// Provide the distance from maze wall along a unit vector.
		// Return a negative value if the current location is outside the maze.
		virtual Number			boundaryDistance(
			Number				locX,			// Current location X coord
			Number				locY,			// Current location Y coord
			Number				hx,				// Unit vector x coord
			Number				hy);			// Unit vector y coord

	protected:
		Number					_sizeX;			// size in x direction
		Number					_sizeY;			// size in y direction
	};

	// ---------------------------------------------------
	// CLASS:	SpatialRegion
	// EXTENDS:	none
	// DESC:	Abstract class to define a region of space
	// RESP:
	//		1.	Given a point determine if the point is
	//			in this region of space (via subclass).
	//		2.	Provide the bounding box, i.e. max and
	//			min X and Y coordinates for the region.
	//
	// NOTE:	A typical use of this hierarchy is to
	//			differentiate regions within a maze.
	// ---------------------------------------------------

	class SpatialRegion  {

	public:

		// Constructor and destructor
		SpatialRegion() {}
		virtual ~SpatialRegion() {}

		// Subclass responsibilities ----------------------------------

		// Get the bounding box of the region
		virtual void			getBoundingBox(
			Number&				xmin, 
			Number&				xmax,
			Number&				ymin,
			Number&				ymax) = 0;

		// Return true if the point supplied is in the region
		virtual bool			containsPoint(
			Number				x,
			Number				y) = 0;
	};

	// ---------------------------------------------------
	// CLASS:	RectangularRegion
	// EXTENDS:	SpatialRegion
	// DESC:	Defines a rectangular region
	// RESP:
	//		1.	Given a point determine if the point is
	//			in this region of space (via subclass).
	//		2.	Provide the bounding box, i.e. max and
	//			min X and Y coordinates for the region.
	// ---------------------------------------------------

	class RectangularRegion : public SpatialRegion {

	public:

		// Constructor and destructor
		RectangularRegion(
			Number				xmin=0,
			Number				xmax=0,
			Number				ymin=0,
			Number				ymax=0);

		virtual ~RectangularRegion();

		// Get the bounding box of the region
		virtual void			getBoundingBox(
			Number&				xmin, 
			Number&				xmax,
			Number&				ymin,
			Number&				ymax);

		// Set the bounding box of the region
		virtual void			setBoundingBox(
			Number				xmin, 
			Number				xmax,
			Number				ymin,
			Number				ymax);

		// Return true if the point supplied is in the region
		virtual bool			containsPoint(
			Number				x,
			Number				y);

	protected:
		Number					_xmin;
		Number					_xmax;
		Number					_ymin;
		Number					_ymax;
	};

	// ---------------------------------------------------
	// CLASS:	CircularRegion
	// EXTENDS:	SpatialRegion
	// DESC:	Defines a circular region
	// RESP:
	//		1.	Given a point determine if the point is
	//			in this region of space (via subclass).
	//		2.	Provide the bounding box, i.e. max and
	//			min X and Y coordinates for the region.
	// ---------------------------------------------------

	class CircularRegion : public SpatialRegion {

	public:

		// Constructor and destructor
		CircularRegion(
			Number				xctr=0,
			Number				yctr=0,
			Number				r=0);

		virtual ~CircularRegion();

		// Parameter accessors

		inline  Number			centerX() { return _centerX; }
		virtual void			centerX(Number x) { _centerX = x; }
		inline  Number			centerY() { return _centerY; }
		virtual void			centerY(Number y) { _centerY = y; }
		inline  Number			radius() { return _radius; }
		virtual void			radius(Number x) { _radius = x; }

		// Get the bounding box of the region
		virtual void			getBoundingBox(
			Number&				xmin, 
			Number&				xmax,
			Number&				ymin,
			Number&				ymax);

		// Return true if the point supplied is in the region
		virtual bool			containsPoint(
			Number				x,
			Number				y);

	protected:
		Number					_centerX;
		Number					_centerY;
		Number					_radius;
	};

};

#endif // #ifndef
