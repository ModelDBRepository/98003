// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_utility.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file is a scratchpad utility driver that is changed as needed.
// Individual utilities are in separately compiled files.

// Function prototypes

// void hoc_reader(int argc, char* args[]); // not included here
void swc_reader(int argc, char* args[]);

// Main driver routine
int main(int argc, char* argv[])
{
	swc_reader(argc, argv);
	return 0;
}
