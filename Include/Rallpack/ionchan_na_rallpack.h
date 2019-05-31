// Sodium Ion Channel
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: ionchan_na_rallpack.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the sodium (Na) channel definitions for the Rallpack benchmark.
 

// Only include this header once
#ifndef __IONCHAN_NA_RALLPACK_H_
#define __IONCHAN_NA_RALLPACK_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace RALLPACK {

	// -------------------------------------------------
	// Na Ion Channel Classes
	// -------------------------------------------------

	// m gate

	class Na_m_gate : public VoltageDepTabGate {

	public:

		// Constructors and destructor
		Na_m_gate() {}
		virtual ~Na_m_gate() {}

		// Alpha and beta functions for table loading
		virtual Number alphaForTable(Number v);
		virtual Number betaForTable(Number v);

		// state vector label functions
		virtual const char* componentName() {return "Na"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "m" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;

	};

	// h gate

	class Na_h_gate : public VoltageDepTabGate {

	public:
		// constructors and destructor
		Na_h_gate() {}
		virtual ~Na_h_gate() {}
	
		virtual Number alphaForTable(Number v);
		virtual Number betaForTable(Number v);

		// state vector label functions
		virtual const char* componentName() {return "Na"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "h" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }


	private:
		static AlphaBetaEntry*		_abTable;
	};

	// Na channel

	class Na_channel : public M3HIonChannel {
	
	public:
		// Constructors and destructor
		Na_channel(Number gSpVal=0);
		virtual ~Na_channel();

		// Reversal potential for Na
		inline Number			Vrev() { return _Vrev; }
		
	protected:
		static const Number		_Vrev;
	};
};


#endif // #ifndef __IONCHAN_NA_RALLPACK_H_

