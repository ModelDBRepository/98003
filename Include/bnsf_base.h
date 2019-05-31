// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_base.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the basic classes used as a
// framework for building biological network simulations.
// As is normal for frameworks, many of these classes are
// abstract and subclassed as needed for the simulation.
//
// This file contains common declarations and is intented 
// to be included before other BNSF headers. See below for
// common documentation conventions also.
//
// To force use of double precision in model computations,
// define BNSF_DOUBLE. See declaration of Number below.



// ====================================================================
// Class documentation example.
//
// The approach taken here is derived from: 
// the Class Responsibility Collaborator (CRC) method
// as described in Beck K, & Cunningham W. 1989. 
// A laboratory for teaching object-oriented thinking.
// Proceeding of OOPSLA 1989: 1-6. 
//
// These comments are basically a CRC card without
// explicitly stating collaborators.
// ====================================================================


// ----------------------------------------------------------------
// CLASS:	ExampleClass
// EXTENDS:	Superclass (i.e. class inherited from)
// DESC:	A description of what type of objects make
//			up this class and generally what they are
//			used for. Abstract classes are those which
//			are not instantiated as such but which may
//			provide data or behavior for subclasses.
// RESP:
//		1.	Numbered list of class responsibilities.
//		2.	The order does not imply sequence of events.
//		3.	Abstract classes can have responsiblities.
//
// NOTES:	Optional notes on the class or implementation.
// ----------------------------------------------------------------



// Only include this header once per compilation unit
#ifndef __BSNF_BASE_H_
#define __BSNF_BASE_H_

// ====================================================================
// MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================
#ifdef WIN32

// Disable warning C4786: symbol greater than 255 character,
#pragma warning( disable: 4786)

#endif
// ====================================================================
// END OF MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================



// ====================================================================
// Header files included here by default
// ====================================================================

// Standard C++ Template Library headers
#include <limits>
#include <string>
#include <iostream>
#include <valarray>

// Incorporate all the names from the std library by reference
using namespace std;


// ====================================================================
// Primary namespace for the framework
// ====================================================================


namespace BNSF {


	// ================================================================
	// Common typedefs and globals
	// ================================================================

	// ----------------------------------------------------------------
	// Unit of Measure Constants (Globals)
	// ----------------------------------------------------------------

	namespace UOM {

		// The general scheme is that conversion constants
		// are multiplicative to follow the normal units
		// notation conventions. For example:
		//
		//   y=x*mm;
		//
		// converts from x which is in units of millimeters to y which
		// is in simulation units. Units of measure can be combined by
		// multiplying and dividing as in r=15*ohm/cm_2, 
		// The compiler should optimize away most of the overhead.

		// To avoid having to write UOM:: before everything, the namespace
		// for UOM can be imported where needed via a "using" specification.
		//
		// See a description of the International System of Units (SI)
		// for the derivation of derived units - see NIST.GOV for details.
		//
		// Because temperature does not convert through simple proportions,
		// separate conversion routines are included below.
		//
		// Double precision floating point is used to ensure accuracy
		// of unit conversions. Because base 10 decimals cannot be precisely
		// represented in binary, there is always some degree of numerical
		// error introduced in scaling by factors of 1/10.


		// Basic units ------------------------------------------------

		// Distance units
		static const double 
			meter=		 1,				// 1 meter
			cm=			 1e-2,			// centimeters
			mm=			 1e-3,			// millimeter
			micron=		 1e-6,			// micro-meter (micron)
			angstrom=	 1e-10;			// Angstrom

		// Time units
		static const double
			minute=		60,				// 1 minute
			sec=		 1,				// 1 second
			msec=		 1e-3,			// millisecond
			microsec=	 1e-6,			// microsecond
			nanosec=	 1e-9,			// nanosecond
			picosec=	 1e-12;			// picosecond

		// Mass units
		static const double
			kg=			 1,				// kilogram
			g=			 1e-3,			// 1 gram
			mg=			 1e-6,			// milli-gram
			microgram=	 1e-12;			// micro-gram

		// Electrical current units
		static const double
			ampere=		 1,				// 1 ampere
			mA=			 1e-3,			// milli-ampere
			microA=		 1e-6,			// micro-ampere
			nanoA=		 1e-9,			// nano-ampere
			picoA=		 1e-12;			// pico-ampere

		// Molecular quantity units
		static const double
			mole=		 1,				// 1 Mole
			mMole=		 1e-3,			// milli-Mole
			microMole=	 1e-6,			// micro-Mole
			nanoMole=	 1e-9,			// nano-Mole
			picoMole=	 1e-12;			// pico-Mole
					
		// Derived units ----------------------------------------------

		// Area units
		static const double 
			meter_2=	meter*meter,			// meter^2
			cm_2=		cm*cm,					// centimeters^2
			micron_2=	micron*micron;			// micron^2

		// Volume units
		static const double 
			meter_3=	meter_2*meter,			// meter^3
			cm_3=		cm_2*cm,				// centimeters^3 (cc)
			micron_3=	micron_2*micron,		// micron^3
			liter=		1000*cm_3,				// liter = decimeter^3
 			mLiter=		cm_3;					// milli-liter = cc

		//  Frequency units
		static const double
			Hz=			1/sec,					// Hertz
			kHz=		1e3/sec,				// kilo-hertz
			MHz=		1e6/sec;				// mega-hertz

		//	Force, work, power and charge units (used below)
		static const double
			newton=		meter*kg/sec/sec,		// SI unit of force
			joule=		newton*meter,			// SI unit of work
			watt=		joule/sec,				// SI unit of power
			coulomb=	ampere*sec;				// SI unit of electric charge

		//  Voltage units
		static const double
			volt=		watt/ampere,			// 1 Volt
			mV=			1e-3*volt;				// millivolt

		// Resistance units
		static const double
			ohm=		volt/ampere,			// 1 ohm
			kohm=		1e3*ohm,				// kilo ohm
			megaohm=	1e6*ohm,				// mega ohm
			gigaohm=	1e9*ohm;				// giga ohm				

		// Conductance units
		static const double
			siemen=		ampere/volt,			// 1 Siemen
			mS=			1e-3*siemen,			// milli Siemen
			microS=		1e-6*siemen,			// micro Siemen
			nanoS=		1e-9*siemen,			// nano Siemen
			picoS=		1e-12*siemen;			// pico Siemen

		// Capacitance units
		static const double
			farad=		coulomb/volt,			// 1 Farad
			microF=		1e-6*farad,				// microfarad
			picoF=		1e-12*farad;			// picofarad

		// Units of chemical concentration
		static const double
			molar=		mole/liter,				// 1 molar
			mM =		1e-3*molar,				// milli-molar
			microM=		1e-6*molar,				// micro-molar
			nanoM=		1e-9*molar;				// nano-molar

		// Other unit conversions -------------------------------------

		static const double
			degreesPerRadian = 57.29577951308232,
			radiansPerDegree = 1/degreesPerRadian;
	
		// ============================================================
		// Physical constants adjusted for units of measure
		// Where possible these values are taken from the NIST
		// web site (note that some texts have obsolete values).
		// 
		// Commonly used single character designations include:
		// F = Faraday's constant
		// k = Boltzmann constant
		// Na = Avogadro's number (atoms per mole)
		// R = gas constant (in theory, R=k*Na, though this is not exact)
		//
		// k and R are per degree Kelvin (not shown below).
		//
		// Some authors use k for the gas constant, e.g. 2vF/kT
		// for a term in the GHK current equation for calcium current.
		// The corresponding term here is 2*v*F/(R*T).
		// ============================================================

		static const double Faraday =  9.64853415e4*coulomb/mole;
		static const double Avogadro = 6.02214199e23/mole;			
		static const double TriplePtH2O = 273.15;
		static const double GasConstant = 8.314472*joule/mole;
		static const double Boltzmann = 1.3806503e-23*joule;
	};

	// ----------------------------------------------------------------
	// The following typedef allows values in a model to be switched
	// between float and double. For low precision results, float may 
	// be sufficient. This typedef allows easy changes as needed. 
	// Note that constants specified as number are double by default. 
	// To avoid warnings, cast them to type Number, e.g. Number(3.14). 
	// If accuracy is not critical, this can be specified as 3.14f.
	// Some computations must be done in extended precision and these are
	// explicitly done in double where they occur regardless of Number.
	// ----------------------------------------------------------------

#ifdef BNSF_DOUBLE
	typedef	 double					Number;			// use double precision numbers
#else
	typedef  float					Number;			// use floating point numbers
#endif

	typedef  valarray<Number>		NumberArray;	// stdlib vector type for numbers
	typedef  Number*				NumberPointer;	// pointer to a number

	// ----------------------------------------------------------------
	// The following typedefs and function allow memory allocation to be 
	// done in units of machine words of some common size. This can be 
	// used to provide consistent boundary alignment if multiple structures
	// are placed into a continuous block of allocated memory.
	//
	// Boundary alignment using float works on many processors but may
	// fail on processors (e.g. some RISC) where doubles must absolutely 
	// always be aligned on 8 byte boundaries. In such a case, MachineWord 
	// should be changed to double and/or appropriate #ifdef tests added.
	// ----------------------------------------------------------------

	typedef	 float					MachineWord;	// unit for boundary alignment
	typedef  char					Byte;			// unit for allocate and free

	inline   unsigned int			sizeInBytesRounded(size_t n) {
		const unsigned int wsz = sizeof(MachineWord);
		return wsz*((n+wsz-1)/wsz);
	}

	// ----------------------------------------------------------------
	// Define a typedef and reserved value for string token identifiers.
	// ----------------------------------------------------------------

	typedef	unsigned int			TokenId;
	static const TokenId			NullTokenId = 0;

	// ================================================================
	// Common global function declarations
	// ================================================================

	// ----------------------------------------------------------------
	// Functions for converting temperature values
	// ----------------------------------------------------------------

	inline Number degKfromC(Number x) { return x + UOM::TriplePtH2O; }
	inline Number degCfromK(Number x) { return x - UOM::TriplePtH2O; }
	inline Number degCfromF(Number x) { return (x-32)/1.8; }
	inline Number degFfromC(Number x) { return 32+x*1.8; }

	// ----------------------------------------------------------------
	// Other common functions
	// ----------------------------------------------------------------

	// Print an error message and exit (two alternate forms)
	void		FatalError(string msg);
	void		FatalError(char* msg);

	// Create a token and return the identifier. Strings that compare
	// as equals return the same unique token id. C strings are converted
	// to string objects automatically.
	TokenId		token(string str);

}; // end of namespace

#endif // #ifndef __BSNF_BASE_H_
