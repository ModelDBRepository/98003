// Ih Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_ih_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// a I-h channel modeled on data from Magee (1998).
//
// Properties modeled are for the case in which [Na] is typical of
// physiological values. Proximal parameters are an extrapolation from
// the measurements where [Na]=0. Deactivation time constants are increased
// as per description by Magee for the [Na]>0 case.
//
// References:
// 
// Maccaferri G, Mangoni M, Lazzari A, Difrancesco D (1993)
// Properties of the hyperpolarization-activated current in rat
// CA1 hippocampal pyramidal cells. J Neurophysiol. 69, 2129-2136.
//
// Magee JC, 1998. Dendritic hyperpolarization-activated currents modify
// the integrative properties of hippocampal CA1 pyramidal neurons.
// J. Neuroscience 18(19), 7613-7624.

 
// Only include this header once
#ifndef __IONCHAN_IH_BAKER_2003_H_
#define __IONCHAN_IH_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace BAKER_2003 {

	// Ih channel (proximal dendrites and soma)

	class Proximal_Ih_channel : public Order1EnergyBarrierTabChannel {
	
	public:
		// Constructors and destructor
		Proximal_Ih_channel(Number gSpVal=0); 
		virtual ~Proximal_Ih_channel() {}

		// Parameter values [Na]>0 estimated based on 8mV difference
		// between proximal and distal Vhalf values (Magee Table 1).
		// gamma is rough fit for time constant with [Na]>0
		virtual Number	ratedTempC()	{return 33; }
		virtual Number	Q10()			{return 4.5; }
		virtual Number	Vhalf()			{return	-73*UOM::mV; }
		virtual Number	slope()			{return -7*UOM::mV; }
		
		// Time constants (see Magee 1998 table 1: Na ext case)
		virtual Number	tauMax()		{return 50*UOM::msec; }
		virtual Number	tauMin()		{return 19*UOM::msec; }
		virtual Number	gamma()			{return 0.6f; }

		// Reversal potential (mixed Na,K ions)
		virtual Number	Vrev()			{return _Vrev; }

		// State vector label functions
		virtual const char*	componentName() {return "Ih"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "q" }; return sl; }
		
	protected:
		static const Number			_Vrev;
		virtual AlphaBetaEntry**	pAlphaBetaTable() { return &_abTable; }


	private:
		static AlphaBetaEntry*		_abTable;
	};

	// Ih channel (distal dendrites)

	class Distal_Ih_channel : public Order1EnergyBarrierTabChannel {
	
	public:
		// constructors and destructor
		Distal_Ih_channel(Number gSpVal=0);
		virtual ~Distal_Ih_channel() {}

		// Parameter values derived  from Magee with [Na+]>0 effects.
		// gamma is rough fit for time constant with [Na+]>0
		virtual Number	ratedTempC()	{return 33; }
		virtual Number	Q10()			{return 4.5; }
		virtual Number	Vhalf()			{return	-81*UOM::mV; }
		virtual Number	slope()			{return -7*UOM::mV; }

		// Time constants (see Magee 1998 table 1: [Na+]>0 case)
		virtual Number	tauMax()		{return 50*UOM::msec; }
		virtual Number	tauMin()		{return 19*UOM::msec; }
		virtual Number	gamma()			{return 0.6f; }

		// State vector label functions
		virtual const char*	componentName() {return "Ih"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "q" }; return sl; }

		// Reversal potential (mixed Na+ and K+ ions)
		inline Number Vrev() { return _Vrev; }
		
	protected:
		static const Number			_Vrev;
		virtual AlphaBetaEntry**	pAlphaBetaTable() { return &_abTable; }


	private:
		static AlphaBetaEntry*		_abTable;
	};

	// Blended Ih channel for intermediate locations

	class Blended_Ih_channel : public Order1BlendedIonChannel {
	
	public:
		// Constructors and destructor
		Blended_Ih_channel(Number gSpVal=0, Number bratio=0.5);
		virtual ~Blended_Ih_channel() {}

		// Reversal potential (mixed Na+ and K+ ions)
		virtual Number Vrev() { return _Vrev; }

		// State vector label functions
		virtual const char*	componentName() {return "Ih"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "q" }; return sl; }

	protected:
		static const Number		_Vrev;
	};

};

#endif // #ifndef __IONCHAN_IH_BAKER_2003_H_
