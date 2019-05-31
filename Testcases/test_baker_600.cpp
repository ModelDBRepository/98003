// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker_600.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Create a network and display connection statistics

#include <iostream>
#include <ctime>

#include "bnsf.h"
#include "mouse_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;

void test_baker_600() 
{
	cout<<"Test case 600 - network creation"<<endl;

	Maze*				mwm = new RectangularMaze;
	Mouse*				mouse = new Mouse(mwm, Mouse::wallFollower);

	// Show statistics
	mouse->printSynapseCounts();

	// Save place cell centers
	cout<<"Saving place field centers"<<endl;
	mouse->ECPlaceCells()->printPlaceFieldCenters("test-baker-ec-pcloc.txt");
	mouse->DGPlaceCells()->printPlaceFieldCenters("test-baker-dg-pcloc.txt");
	mouse->CA3PlaceCells()->printPlaceFieldCenters("test-baker-ca3-pcloc.txt");

	cout<<"Deleting data structures"<<endl;
	delete mouse;
	delete mwm;
	cout<<"Done"<<endl;
}
