// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: test_baker.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file is a scratchpad test driver that is changed as needed.
// Individual test cases are in separately compiled files.

// Function prototypes
void test_baker_010();			// Individual ion channel cases
void test_baker_020();

void test_baker_100();			// Cell current injection test cases
void test_baker_110();

void test_baker_210();			// Synapse test case

								// Synaptic plasticity test cases (STDP)
void test_baker_300();			// - single paired spikes
void test_baker_310();			// - homosynaptic LTD
void test_baker_320();			// - frequency dependencies

void test_baker_350();			// Synaptic plasticity (presynaptic)

								// Mouse and place cell test cases
void test_baker_600();			// - test cell and network building
void test_baker_610();			// - test afferent spike trains
void test_baker_670();			// - test theta phase precession


int main(int argc, char* args[])
{
	test_baker_670();			// Change as needed to select desired test
	return 0;
}
