// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: swc_reader.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file contains a utility to read an SWC format morphology file
// and convert it into an easier format to import into a C++ program.
// Only the data table portion is generated here. Appropriate C++
// declarations should be added manually.
//
// X,Y,Z coordinates are estimated along the branch and do not reflect
// the actual track of the associated neurite process. Distance from
// soma does reflect the actual track except where excessive changes
// in coordinates are encountered along the way.
//
// A non-C++ format (.csv or .txt) can also be written for importing 
// the morphology into visualization tools.
//
// Output fields are the following (see BNSF::MorphologyEntry)
//
//		int				type;			-- type of entry
//		int				idnum;			-- identifier of this entry (should = index)
//		int				parent;			-- identifier of parent entry in tree structure
//		int				branch;			-- arbitrary identifier of the branch
//		double			r;				-- radius in microns
//		double			len;			-- length in microns
//		double			dist;			-- distance from soma in microns
//		double			x;				-- X coordinate of this compartment
//		double			y;				-- Y coordinate of this compartment
//		double			z;				-- Z coordinate of this compartment
//		double			dx;				-- Orientation vector X coord
//		double			dy;				-- Orientation vector Y coord
//		double			dz;				-- Orientation vector Z coord
//
// Type values used follow SWC types.
//
//		somaType=1,				axonType=2, 
//		basalDendriteType=3,	apicalDendriteType=4
//		end marker = -1
//
// References:
//
// See the Duke/Southamptom Cell Archive at http://neuron.duke.edu/cells/
// for information about SWC files. The README file associated with
// the program CVAPP has further information about the SWC format.

// Standard C++ Template Library headers
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

// C headers embedded in the std namespace
#include <cstdio>
#include <cmath>

// Incorporate all the names from the std library by reference
using namespace std;

// Internal classes used here

class swcEntry {		// Represents one point read from swc file
public:
	int					type;		// swc type -- see type const below
	float				x,y,z;		// location of point 
	float				r;			// radius at point
	int					parent;		// parent index
	vector<int>			children;	// children (inverse of parent relationship)
	float				dist;		// distance from soma
	bool				isBP;		// is branch point
	int					branch;		// branch containing this point
};

class branchEntry {		// Represents points between end points
public:
	int					type;		// types of points on this branch
	vector<int>			points;		// indexes in points vector
	float				length;		// total length
	float				area;		// total membrane area
	float				radius;		// average radius
	float				distToSoma;	// distance from soma to first point
	int					parent;		// parent branch
	int					nseg;		// number of segments to generate
	int					compNbr;	// number of compartment for 1st segment
};

// Main routine for this utility
void swc_reader(int argc, char* argv[])
{
	cout<<"SWC format conversion"<<endl<<endl;


	// Parameters controlling compartment generation.
	// Lengths are in microns.
	const float		maxCompLen = 50;	// Max compartment size (a few may be larger)
	const float		rChgRelTol = 0.25f;	// Fractional change allowed in radius along branch
	const float		rChgTolLen = 5;		// Minimum branch length when dividing for radius change
	const float		rMaxForDendrite = 5.0;	// Max dendrite size next to soma (larger are merged into soma)	

	const bool		useYAxisForSomaArea = false;	// Assume soma is aligned on Y axis
	const bool		skipAxon = true;				// Skip writing the axon
	const bool		debugBranchRadiusChg = false;	// Notify on each change

	// Input and output file names
	char*			inFileName;
	char*			outFileName;

	// If args provided use them as names.
	// Otherwise use hardcoded names (simplifies change for testing)
	if (argc==3) {
		inFileName = argv[1];
		outFileName = argv[2];
	}
	else if (argc!=1) {
		cerr<<"Usage: swc_reader <infile> <outfile>"<<endl;
	}
	else {

		inFileName = "l56a.swc";
		outFileName = "l56a-50.cpp";
	}
	cout<<"Using the following file names"<<endl
		<<"Input file = "<<inFileName<<endl
		<<"Output file = "<<outFileName<<endl<<endl;

	// Flag indicating whether to generate output in C++ format or not. 
	// Set false if last 4 letters of output file is .csv or .txt. 
	// Otherwise set to true.
	int				ofnl = strlen(outFileName);
	bool			useCPPFormat =
						ofnl>=4 && 
						strcmp(&outFileName[ofnl-4],".csv")!=0 &&
						strcmp(&outFileName[ofnl-4],".txt")!=0;

	// Parameters controlling error correction of the input
	const float		maxXJump = 30;		// Maximum jump in x (microns)
	const float		maxYJump = 30;		// Maximum jump in y (microns)
	const float		maxZJump = 30;		// Maximum jump in z (microns)

	// Constants for swc types
	const int		somaType = 1;
	const int		axonType = 2;
	const int		dendriteType = 3;
	const int		apicalDendriteType = 4;
	const int		otherType = 10;

	// Other constants
	const float		pi = 3.14159f;

	// Counters and statistics
	int				pointsRead = 0;
	int				jumpCorrections = 0;
	int				typeCorrections = 0;
	int				pointsMergedIntoSoma = 0;
	int				branchesRead = 0;
	int				compWritten = 0;
	int				dzCount = 0;
	float			dz1stMoment = 0;

	// Input & output files
	FILE*			in;
	FILE*			out;

	// Working variables
	const int		inputLineLen = 256;
	char			inputLine[inputLineLen];
	char*			pchar;
	swcEntry		point;
	branchEntry		branch, emptyBranch;

	int				fcnt;
	int				next,num;
	int				i,j,k,k0,k1,n,p,child;
	int				errCnt;

	float			diam,len,area,chgRatio;
	float			x,y,z;
	float			dx,dy,dz,dist;
	bool			jumpError;
	bool			stopLoop;
	bool			sphericalSoma;

	// Data describing the morphology
	vector<swcEntry>		points;
	vector<branchEntry>		branches;

	// Size of the soma in terms of a cylinder
	float			somaLength;
	float			somaRadius;
	float			somaArea;

	// Open input and output files
	in=fopen(inFileName,"rb");
	if (in==NULL) {
		perror("Input file open failed");
		exit(1);
	}

	out=fopen(outFileName,"w+");
	if (out==NULL) {
		perror("Output file open failed");
		exit(1);
	}

	// Initialize point values not being set during the first read pass
	point.dist = 0;
	point.branch = -1;
	point.isBP = false;

	// Process the input file (first data pass)
	next=1;
	for(;;) {

		// Get the next input line
		if (NULL==fgets(inputLine,inputLineLen-1,in)) {
			if (ferror(in)) {
				perror("Error reading from input file");
				exit(1);
			}
			break;
		}

		// Copy any comment lines directly to the output.
		// First locate the first non-blank character.
		for (pchar=inputLine;pchar!='\0' && isspace(*pchar); pchar++);
		if (*pchar=='\0' || *pchar=='\n' || *pchar=='\r')
			continue;
		if (*pchar=='#') {
			if (useCPPFormat) {
				// Write out as a C++ comment
				fputs("//",out);
				for (pchar++; *pchar!='\0'; pchar++) {
					if (!iscntrl(*pchar)) {
						fputc(*pchar,out);
					}
				}
				fputc('\n',out);
			}
			continue;
		}

		// Parse the input line
		fcnt=sscanf(inputLine,"%d%d%g%g%g%g%d",
			&num,
			&point.type,
			&point.x,
			&point.y,
			&point.z,
			&point.r,
			&(point.parent) );
		if (fcnt!=7) {
			cerr<<"Error reading swc entry. Input line follows below."<<endl;
			cerr<<inputLine<<endl;
			exit(1);
		}

		if (next!=num) {
			cerr<<"Index number read was not next in sequence"<<endl
				<<"Read: "<<num<<" but expected: "<<next<<endl;
			exit(1);
		}
		next++;				// set for next read
		pointsRead++;		// count the point

		// Adjust the parent number to reflect starting with 0
		// Parent = -1 indicates the root of the tree.
		// This should only occur on the first entry read.
		if (point.parent>0) {
			if (num==1) {
				cerr<<"First entry does not have parent=-1"<<endl;
				exit(1);
			}
			if (num<=point.parent) {
				cerr<<"Parent does not precede child in index order."
					" Child index = "<<num
					<<"Parent index ="<<point.parent<<endl;
				exit(1);
			}
			point.parent--;

			// Make sure that if this is part of the soma, it is
			// not a child of something that is not in the soma.
			if (point.type==somaType &&
				points[point.parent].type!=somaType) {
				cerr<<"Soma point has non-soma parent. Index="<<num;
				
				point.type=points[point.parent].type;
				cerr<<" New type="<<point.type<<endl;
			}
		}
		else {
			if (num!=1) {
				cerr<<"Parent<=0 found on entry other than first"<<endl;
				exit(1);
			}
			if (point.type!=somaType) {
				cerr<<"Tree root in first entry is not at the soma"<<endl;
				exit(1);
			}
			point.parent= -1;
		}

		// Put the line into the vector of points read
		points.push_back(point);
	}

	// Perform a second, in memory, pass of points read.
	// Skip the first entry, which is the tree root.
	for (k=1;k<points.size();k++) {

		// Debug
		point=points[k];

		// Skip other entries
		if (points[k].type==otherType) {
			continue;
		}
		
		// Invert in parent-child relationship.
		p = points[k].parent;
		if (p<0 || p>=points.size()) {
			cerr<<"Item "<<k+1<<" parent is out of range"<<endl;
			exit(1);
		}
		points[p].children.push_back(k);

		// Compute distance from soma root to current point.
		// The maximum change in z is capped because
		// there can be measurement errors in z whenever
		// plane of focus is changed during the reconstruction.
		// Similarly, errors in x or y can arise when joining
		// separate images.
		dx=points[k].x-points[p].x;
		dy=points[k].y-points[p].y;
		dz=points[k].z-points[p].z;

		dzCount++;
		dz1stMoment += fabs(dz);

		if (jumpError= (
			fabs(dx)>maxXJump ||
			fabs(dy)>maxYJump ||
			fabs(dz)>maxZJump )) {
			if (jumpCorrections==0) {
				cerr<<"Warning: coordinate jump corrections at the following indexes:"<<endl;
			}
			cerr<<"\tindex="<<k+1;
		}
		if (fabs(dx)>maxXJump) {
			cerr<<" dx="<<dx;
			jumpCorrections++;
			dx=0;
		}
		if (fabs(dy)>maxYJump) {
			cerr<<" dy="<<dy;
			jumpCorrections++;
			dy=0;
		}
		if (fabs(dz)>maxZJump) {
			cerr<<" dz="<<dz;
			jumpCorrections++;
			dz=0;
		}
		if (jumpError) {
			cerr<<endl;
		}
		points[k].dist=points[p].dist+sqrt(dx*dx+dy*dy+dz*dz);
	}

	// Fix unknown types (-1 occurs now and then).
	// Use either parents or children to find valid type.
	errCnt = 0;
	for (k=1;k<points.size();k++) {

		// Debug
		point=points[k];

		// Skip other entries
		if (points[k].type==otherType) {
			continue;
		}

		// Check for valid type
		p = points[k].parent;
		if (points[k].type != somaType &&
			points[k].type != axonType &&
			points[k].type != dendriteType &&
			points[k].type != apicalDendriteType) {

			// Look to the parent for a type, but only
			// if the parent is not the soma.				
			if (points[p].type!=somaType && points[p].type!=-1) {
				cerr<<"Warning: index "<<k+1<<" unknown type ("<<points[k].type
					<<") changed to that of parent ("<<points[p].type<<")"<<endl;
				points[k].type=points[p].type;
			}
			else {
				// Set to a known invalid value
				points[k].type = -1;
				errCnt++;
			}
		}
		else {
			// See if the parent type was left unresolved
			// and if so, go back and fix it.
			while( p>0 && points[p].type==-1) {
				cerr<<"Warning: index "<<p+1<<" unknown type ("<<points[p].type
					<<") changed to "<<points[k].type<<endl;
				points[p].type = points[k].type;
				errCnt--;
				p = points[p].parent;
			}
		}
	}
	if (errCnt!=0) {
		cerr<<"Unknown types were left unresolved at the following indexes:"<<endl;
		for (k=0;k<points.size();k++) {
			if (points[k].type == -1) {
				cerr<<k+1<<endl;
			}
		}
		exit(1);
	}

	// Check that parents are either the same type as child or else the soma
	for (k=1;k<points.size();k++) {

		// Check for valid type
		p = points[k].parent;
		if (points[p].type!=somaType && 
			points[p].type!=points[k].type &&
			(!skipAxon  || points[k].type!=axonType) ) {
			cerr<<"Warning: index "<<k+1<<" type ("<<points[k].type
				<<") does not match parent type ("<<points[p].type<<")"<<endl;
		}
	}

	// Generate a cylinder of the same area as the soma.
	// In most cases, the soma is modelled as a single equipotential
	// compartment, and the actual geometry is not important.
	// In addition, check that the soma is fully connected.
	// This translates to every parent being in the soma.
	// Also, look for dendrites with radius below rMaxForDendrite.
	// If these are adjacent to the soma, include them in it.

	somaLength = 0;
	somaRadius = 0;
	somaArea = 0;
	
	// Look for soma points. Note that the tree root is skipped.
	// It will be accounted for through its children.
	// The case where the soma is a single point is detected below when 
	// the resulting length is zero.

	// If useYAxisForSomaArea is true, area is computed assuming radius
	// for soma points reflects the radius in a direction orgthogonal to
	// the Y axis. This is an unstated assumption in some reconstructions.

	for (k=1;k<points.size();k++) {

		// Debug
		// point=points[k];

		// Skip other entries
		if (points[k].type==otherType) {
			continue;
		}

		// Skip points that are not in the soma
		p = points[k].parent;
		len = useYAxisForSomaArea
			? fabs(points[k].y-points[p].y)
			: points[k].dist-points[p].dist;
		diam = points[k].r+points[p].r;
		area = pi*diam*len;

		if (points[k].type!=somaType) {
			if (points[k].r<=rMaxForDendrite || points[p].type!=somaType) {
				continue;
			}

			// Record merged point but do not count in soma length
			pointsMergedIntoSoma++;
			cerr<< "Point at index "<<k+1
				<< " merged into soma. radius = "<<points[k].r
				<< " len = "<<len
				<< endl;
			len=0;

		}

		// Check parent -- this should be ok, but be certain
		if (points[points[k].parent].type != somaType) {
			cerr<<"Soma element at index "<<k+1<<" has non-soma parent"<<endl;
			exit(1);
		}

		// Include all found points in the soma. Since parents
		// come before children, this will recursively handle
		// points topologically adjacent to the soma.
		points[k].type=somaType;

		// Add up the length and area that go with the segment
		// from the current point to its parent. With different
		// geometries, the resulting cylinder might not be
		// especially meaningful in terms of length versus radius,
		// but the area will be preserved.
		somaLength += len;
		somaArea += area;
	}

	// Make sure that valid segments for the soma were found.
	// Otherwise, use the first entry as the whole soma.
	if (somaLength!=0) {

		// Get average radius from area and length
		somaRadius = somaArea / ( 2*pi*somaLength);
		sphericalSoma = false;
	}
	else {
		cerr<<"Warning: soma did not contain multiple points."<<endl;
		cerr<<"The first point is taken to define the soma as a ball"<<endl;
		somaRadius = points[0].r;
		somaLength = 2*points[0].r;
		somaArea = 4*pi*points[0].r*points[0].r;
		sphericalSoma = true;
	}

	// Make a prototype empty branch for use in building table below.
	emptyBranch.area = 0;
	emptyBranch.length = 0;
	emptyBranch.distToSoma = 0;

	// Find branch points and build table of points on each branch.
	// For now, the branch is modeled as a single cylinder
	// of equivalent length and radius to preserve area. The
	// maximum variation in radius is limited to rChgTol so that
	// average radius will not be too far off at any point.
	for (k=1;k<points.size();k++) {

		// Debug
		point=points[k];

		// Skip other entries
		if (points[k].type==otherType) {
			continue;
		}
		// Skip somatic branch points and also
		// axon branches if they will not be written.
		if (points[k].type==somaType ||
			(skipAxon && points[k].type == axonType) ) {
			continue;
		}

		// Skip points along the branch unless there is a change
		// in radius of sufficient size to justify treating this
		// as a branch point. This logic looks at the radius of
		// the immediate child and ensures that no earlier point
		// along the branch changes too much from the child's radius. 
		// The branch point is then the last point at which the radius 
		// change tolerance was still satisfied.
		if (points[k].children.size() == 1) {
			child=points[k].children[0];
			k1=k;
			stopLoop=false;
			while (points[k1].type!=somaType && !points[k1].isBP) {
				dist = points[k].dist - points[k1].dist;
				chgRatio = points[child].r / points[k1].r;
				if (dist>rChgTolLen && (chgRatio>1+rChgRelTol || chgRatio<1-rChgRelTol)) {
					if (debugBranchRadiusChg) {
						cerr<<"Branch radius change"
							<<" from "<<points[k1].r<<" at "<<k1+1
							<<" to "<<points[child].r<<" at "<<child+1
							<<endl;
					}
					stopLoop=true;
					break;
				}
				k1=points[k1].parent;
			}
			if (!stopLoop) {
				continue;
			}
		}
		
		// The current point is the last point in the branch
		points[k].isBP = true;

		// Reset the work area for a branch
		branch = emptyBranch;

		// Set the type based on the last point in the branch
		branch.type = points[k].type;

		// Create a new branch and trace its contents
		// by following the parent pointers. Accumulate
		// length and area for the branch along the way.
		k0=k;
		do {
			// Save point and set pointer from point to branch
			branch.points.push_back(k0);
			points[k0].branch=branches.size();

			// Compute distance and area for segment up to parent
			k1 = points[k0].parent;
			dist = points[k0].dist - points[k1].dist;
			if (points[k1].type==somaType) {
				diam = 2*points[k0].r;
				// If soma is a sphere, distance should not count soma radius
				if (sphericalSoma) {
					dist -= somaRadius;
				}
			}
			else {
				// For an n-way branch point parent, the parent radius
				// is not used in computing an average diameter
				if (points[k1].children.size()>1) {
					// Use only the current points radius
					diam = 2*points[k0].r;
				}
				else {
					// Otherwise, get diameter as twice the average radius
					// of adjacent points, i.e. use trapezoid rule.
					diam = points[k0].r +points[k1].r;
				}
			}
			branch.length += dist;
			branch.area += pi*diam*dist;

			// Move one step down the parent chain
			k0 = k1;
		}
		while (points[k0].type!=somaType && !points[k0].isBP);

		// Since the points were inserted backwards, reverse them
		reverse(branch.points.begin(),branch.points.end());

		// Save the parent of this branch
		p = points[k0].branch;
		branch.parent = p;
		if (p!=-1) {
			branch.distToSoma = branches[p].distToSoma+branches[p].length;
		}

		// Get an average radius for this branch
		branch.radius = branch.area/(2*pi*branch.length);

		// Put the current branch in the vector of branches
		branches.push_back(branch);
		branchesRead++;
	}

	// Make a pass through all branches to get number of segments
	n = 1;
	for (k=0;k<branches.size();k++) {

		// Get number of segments to use and compartment number
		// Compartment 0 is reserved for the soma.
		branches[k].nseg = int( branches[k].length/maxCompLen+0.9);
		if (branches[k].nseg == 0) {
			branches[k].nseg = 1;
		}
		branches[k].compNbr = n;
		n += branches[k].nseg;
	}

	// Write out compartments starting with the soma
	if (useCPPFormat) {
		fprintf(out,"{");
	}
	fprintf(out,
		"%2d,%4d,%4d,%5d,%8.3f,%8.3f,%9.3f, "
		"%6.1f,%6.1f,%6.1f, %5.1f,%5.1f,%5.1f",
		somaType,					// type = soma
		0,							// id of the soma
		-1,							// parent of the soma
		-1,							// no branch for the soma
		somaRadius,					// radius
		somaLength,					// length
		0.0f,						// distance to soma
		points[0].x,points[0].y,points[0].z,	// soma center location (x,y,z)
		0.0f, 0.0f, 0.0f);			// orientation vector (dx,dy,dz)

	if (useCPPFormat) {
		fprintf(out,"},");
	}
	fprintf(out,"\n");

	compWritten++;

	// Write the rest of the compartments
	for (k=0; k<branches.size(); k++) {

		// Debug
		branch = branches[k];

		// Put the number of the parent compartment in n
		p = branches[k].parent;
		if (p>=0) {
			n=branches[p].compNbr+branches[p].nseg-1;
		}
		else {
			n = 0;
		}

		// Find starting point of the  branch x,y,z coordinates.
		if (sphericalSoma) {
			// Use the first point in the branch as its origin
			j = branches[k].points[0];
		}
		else {
			// Use the parent of the first point in this branch
			j = branches[k].points.front();
			j = points[j].parent;
		}
		x = points[j].x;
		y = points[j].y;
		z = points[j].z;

		// Get difference per segment between start and end of branch
		i = branches[k].points.back();
		dx = (points[i].x-x) / branches[k].nseg;
		dy = (points[i].y-y) / branches[k].nseg;
		dz = (points[i].z-z) / branches[k].nseg;

		// Adjust location to middle for first compartment
		x += dx/2;
		y += dy/2;
		z += dz/2;

		// Get length of each compartment in the branch
		len = branches[k].length / branches[k].nseg;

		// Generate each compartment in the branch
		for (j=0;j<branches[k].nseg;j++) {

			if (useCPPFormat) {
				fprintf(out,"{");
			}
			fprintf(out,
				"%2d,%4d,%4d,%5d,%8.3f,%8.3f,%9.3f, "
				"%6.1f,%6.1f,%6.1f, %5.1f,%5.1f,%5.1f",
				branches[k].type,					// type of points in branch
				branches[k].compNbr+j,				// index of this entry
				n,									// parent
				i+1,								// swc id of branch
				branches[k].radius,					// radius
				len,								// length
				branches[k].distToSoma+j*len+len/2,	// distance to soma
				x,y,z,								// location of this entry
				dx,dy,dz);							// orientation of compartment
			if (useCPPFormat) {
				fprintf(out,"},");
			}
			fprintf(out,"\n");
			compWritten++;

			// Update compartment location along a line from start to end points
			x += dx;
			y += dy;
			z += dz;

			// Make this compartment the next parent
			n = branches[k].compNbr+j;
		}
	}

	// Write a final end marker (C++ format only)
	if (useCPPFormat) {
		fprintf(out,"\n// end marker\n{");
		fprintf(out,
			"%2d,%4d,%4d,%5d,%8.3f,%8.3f,%9.3f, "
			"%6.1f,%6.1f,%6.1f, %5.1f,%5.1f,%5.1f",
			-1,							// type of points in branch
			-1,							// index of this entry
			-1,							// parent
			-1,							// swc id of branch
			0.0,						// radius
			0.0,						// length
			0.0,						// distance to soma
			0.0, 0.0, 0.0,				// location of this entry
			0.0, 0.0, 0.0);				// location of this entry
		fprintf(out,"}\n");
	}

	// Write some statistics
	cout<<endl;
	cout<<"Number of input points = "<<pointsRead<<endl;
	cout<<"Number of branches read = "<<branchesRead<<endl;

	cout<<"Mean of delta z magnitude (before fixes) = "<<dz1stMoment/dzCount<<endl;	
	cout<<"Number of dx/dy/dz jump fixes = "<<jumpCorrections<<endl;
	cout<<"Number of points merged into soma = "<<pointsMergedIntoSoma<<endl;

	if (useYAxisForSomaArea)	cout<<"Soma is assumed to be Y axis aligned";
	else						cout<<"Soma is not assumed to be Y axis aligned";
	cout<<endl<<"Soma area = "<<somaArea<<" micron^2"<<endl;

	cout<<"Number of compartments written = "<<compWritten<<endl;

	// Wrap up processing
	fclose(in);
	fclose(out);

}
