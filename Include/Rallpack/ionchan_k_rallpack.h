// K Ion Channel
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: ionchan_k_rallpack.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K channel definition used in Rallpack benchmarks.
 

// Only include this header once
#ifndef __IONCHAN_K_RALLPACK_H_
#define __IONCHAN_K_RALLPACK_H_
#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace RALLPACK {

	// --------------------------------------------
	// K Ion Channel class
	// --------------------------------------------

	class K_channel : public VoltageDepTabChannel {

	public:

		// constructors and destructor
		K_channel(Number gScaled=0) : VoltageDepTabChannel(gScaled) {}
		virtual ~K_channel() {}

		// reversal potential
		inline Number Vrev() { return _Vrev; }

		// alpha beta computations
		virtual Number alphaForTable(Number v);
		virtual Number betaForTable(Number v);

		// state vector label functions
		virtual const char* componentName() {return "K"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "n" }; return sl; }

		// required functions
		virtual Number conductance();

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }
		static const Number			_Vrev;

	private:
		static AlphaBetaEntry*		_abTable;
	};
};


#endif // #ifndef __IONCHAN_K_RALLPACK_H_

