// Rallpack test driver
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: test_rallpack.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file provides driver code for rallpack tests.

// Function prototypes

void test_rallpack_ab();	// Get alpha-beta values for ion channels

void test_rallpack1();		// Test cases -- see assoc Genesis files.
void test_rallpack2();
void test_rallpack3();

// Mainline code
int main(int argc, char* args[])
{
	test_rallpack_ab();
	return 0;
}
