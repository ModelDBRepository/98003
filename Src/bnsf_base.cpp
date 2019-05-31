// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_base.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file provides the body of basic BNSF classes and functions.



#include "bnsf_base.h"
#include <map>

using namespace std;
using namespace BNSF;


// ====================================================================
// Global functions
// ====================================================================


// Print a fatal error message and exit.
// This is also a good place to put a
// breakpoint when debugging.
void BNSF::FatalError(char* msg)
{
	cerr << msg << endl;
	exit(1);
}

// Print a fatal error message and exit.
// This is also a good place to put a
// breakpoint when debugging.
void BNSF::FatalError(string msg)
{
	cerr << msg << endl;
	exit(1);
}

// Return a unique identifier for the string supplied.
TokenId BNSF::token(string str)
{
	typedef map<string,unsigned int>	TokenMap;
	typedef TokenMap::iterator			TokenMapIt;
	typedef TokenMap::value_type		TokenEntry;

	// Include the last id and map here as a static to ensure correct order
	// of initialization if tokens are assigned during static initialization.
	// This is needed since the order of such initializations is undefined
	// except within a single file. 
	
	// This is obviously not thread safe as coded.

	static TokenId						lastId = 0;
	static TokenMap						tokens;

	TokenMapIt							it;

	it=tokens.find(str);
	if (it!=tokens.end() ) {
		return it->second;
	}
	else {
		if (++lastId==0) {
			FatalError("(BNSF::token) Too many tokens assigned");
		}
		tokens.insert(TokenEntry(str,lastId));
		return lastId;
	}
}







