// Calcium Ion Channel Dynamics
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_ca_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the Ca++ channel definitions. Only the most relevant
// channel types are included. Ca-N as modeled here is probably 
// a mixture of medium voltage activated Ca++ channel types. 
//
// Channel models are for the most part taken from Magee & Johnston
// except for Ca-T which is a composite of various results. See
// Metz et al. for Ca++ currents associated with ADP.
//
// Kinetics are drawn from a variety of related studies as noted.
// Per Magee & Johnston, a 15mV adjustment is needed to compensate
// for the use of 20mM Ba++ versus 2mM Ca++ found in vivo.
//
// Order 2 activation gates are used similar to the scheme of Jaffe. 
// Q10 values are very rough estimates, though are consistent with findings
// by Takahashi for Ca-T channels.
//
// References:
//
// Avery RB and Johnston D (1996). Multiple channel types contribute
// to the low-voltage-activated calcium current in hippocampal CA3
// pyramidal neurons. J. Neuroscience 16(18), 5567-5582.
//
// Brown AM, Schwindt PC, and Crill WE (1993). Voltage
// dependence and activation kinetics of pharmacologically
// defined components of the high-threshold calcium currents
// in rat neocortical neurons. J. Neurophysiol. 70(4), 1530-1543.
//
// Fisher RE, Gray R, and Johnston D. (1990). Properties and
// distribution of single voltage-gated calcium channels in
// adult hippocampal neurons. J. Neurophysiology 64, 91-104.
//
// Jaffe DB, Ross WN, Lisman JE, Lasser-Ross N, Miyakawa H, and
// Johnston D. (1994) A model for dendritic Ca++ accumulation in
// hippocampal pryamidal neurons based on fluorescence imaging
// measurements. J. Neurophysiol. 71(3) 1065-1077.
//
// Jung H-y, Staff NP, and Spruston N (2001). Action potential 
// bursting in subicular pyramidal neurons is driven by a calcium
// tail current. J. Neuroscience 21(10), 3312-3321.
//
// Magee JC, Johnston D (1995). Characterization of single
// voltage-gated Na+ and Ca++ channels in apical dendrites
// of rat CA1 pyramidal neurons. J. Physiology 487, 67-90.
//
// Magee JC (1999). Voltage gated ion channels in dendrites, in
// Dendrites, editted by Greg S and Spruston N. 
// New York: Oxford University Press (2001 reprint).
//
// Metz AE, Jarsky T, Martina M, and Spruston N (2005). R-type
// calcium channels contribute to afterdepolarization and
// bursting in hippocampal CA1 pyramidal neurons. 
// J. Neuroscience 25(24), 5763-5773.
//
// Kavalali, ET, Zhuo M, Bito  H, and Tsien RW (1997).
// Dendritic Ca++ channels characterized by recordings from
// isolated hippocampal dendritic segments. Neuron 18, 651-663.
//
// Takahashi K, Ueno S, and Akaike N. (1991). Kinetic properties
// of T-type Ca++ currents in isolated rat hippocampal CA1
// pyramidal neurons. J. Neurophysiology 65, 148-155.



// Only include this header once
#ifndef __IONCHAN_CA_BAKER_2003_H_
#define __IONCHAN_CA_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace BAKER_2003 {

	// ================================================================
	// Ca-T channel classes
	// ================================================================

	// --------------------------------------------
	// Ca-T m gate class
	// --------------------------------------------

	class Ca_T_m_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		Ca_T_m_gate() {}
		virtual ~Ca_T_m_gate() {}

		// Parameter values. Q10 is from Takahashi.
		virtual Number	ratedTempC()	{return  22; }
		virtual Number	Q10()			{return	 1.7f; }

		// Order 2 gate for Vhalf=-30mV, k=6mV (see Avery & Johnston)
		// This is similar to results from Magee & Johnston, but
		// without a -15mV adjustment. See also Kavalali et al.
		virtual Number	Vhalf()			{return	-36.2*UOM::mV;  }
		virtual Number	slope()			{return  7.0*UOM::mV; }

		// Time constants are derived from Takahashi time to peak.
		// tauMin is from tail currents in Kavalali for an order 2 gate.
		virtual Number	tauMax()		{return  20*UOM::msec; }
		virtual Number	tauMin()		{return  3*UOM::msec; }
		virtual Number	gamma()			{return	 0.7f; }

		// State vector label functions
		virtual const char* componentName() {return "CaT"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "m" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }
		static const Number			_Vrev;

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// --------------------------------------------
	// Ca-T h gate class
	// --------------------------------------------

	class Ca_T_h_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		Ca_T_h_gate() {}
		virtual ~Ca_T_h_gate() {}

		// Parameter values. Q10 is from Takahashi.
		virtual Number	ratedTempC()	{return  22; }
		virtual Number	Q10()			{return	 2.5f; }

		// Voltage sensitivity (see Avery & Johnston).
		// This is similar to results from Magee & Johnston
		// with a -15mV adjustment. See also Kavalali et al.
		virtual Number	Vhalf()			{return	-80*UOM::mV; }
		virtual Number	slope()			{return -6.4*UOM::mV; }

		// Time constants (see Takahashi).
		virtual Number	tauMax()		{return  35*UOM::msec; }	
		virtual Number	tauMin()		{return  10*UOM::msec; }
		virtual Number	gamma()			{return  0.8f; }

		// State vector label functions
		virtual const char* componentName() {return "CaT"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "h" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// --------------------------------------------
	// Ca-T Ion Channel class
	// --------------------------------------------

	class Ca_T_channel : public M2HCaIonChannel {

	public:
		// constructors and destructor
		Ca_T_channel(Number gSpVal=0);
		virtual ~Ca_T_channel() {}
	};


	// ================================================================
	// Ca-N channel classes (various MVA Ca++ chan)
	// ================================================================

	// --------------------------------------------
	// Ca-N m gate class
	// --------------------------------------------

	class Ca_N_m_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		Ca_N_m_gate() {}
		virtual ~Ca_N_m_gate() {}

		// Parameter values
		virtual Number	ratedTempC()	{return  22; }
		virtual Number	Q10()			{return	 2; }

		// Order 2 gate for Vhalf=-12mV, k=8mV (-15mV adjustment)
		// from Magee & Johnston 1995.
		virtual Number	Vhalf()			{return	-20.3*UOM::mV; }
		virtual Number	slope()			{return  9.4*UOM::mV; }
		virtual Number	xinfExponent()	{return  2; }

		// Time constants (see Brown et al. 1993).
		// Note that kinetics do not follow independent gating model
		// and that order 2 is for voltage sensitivity only.
		virtual Number	tauMax()		{return  1.2*UOM::msec; }	
		virtual Number	tauMin()		{return  0.4*UOM::msec; }

		// State vector label functions
		virtual const char* componentName() {return "CaN"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "m" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// --------------------------------------------
	// Ca-N h gate class
	// --------------------------------------------

	class Ca_N_h_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		Ca_N_h_gate() {}
		virtual ~Ca_N_h_gate() {}

		// Parameter values
		virtual Number	ratedTempC()	{return  22; }
		virtual Number	Q10()			{return	 2; }

		// Voltages from Magee & Johnston (-15mV adjusted)
		virtual Number	Vhalf()			{return	 -54*UOM::mV; }
		virtual Number	slope()			{return  -9.2*UOM::mV; }

		// Time constants as estimated (see Magee & Johnston)
		virtual Number	tauMax()		{return  100*UOM::msec; }
		virtual Number	tauMin()		{return  100*UOM::msec; }
		
		// State vector label functions
		virtual const char* componentName() {return "CaN"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "h" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// --------------------------------------------
	// Ca-N Ion Channel class
	// --------------------------------------------

	class Ca_N_channel : public M1HCaIonChannel {
	
	public:
		// constructors and destructor
		Ca_N_channel(Number gSpVal=0);
		virtual ~Ca_N_channel() {}

	};

	// ================================================================
	// Ca-L channel class
	// ================================================================

	class Ca_L_channel : public Order1CaEnergyBarrierTabChannel {
	
	public:
		// Constructors and destructor
		Ca_L_channel(Number gSpVal=0);
		virtual ~Ca_L_channel() {}

		// Parameter values (from Magee & Johnston 1995)
		virtual Number	ratedTempC()	{return  22; }
		virtual Number	Q10()			{return	 2; }

		// Order 2 gate for Vhalf=-6mV, k=6mV (-15mV adjustment)
		// from Magee & Johnston 1995.
		virtual Number	Vhalf()			{return	-12.2*UOM::mV; }
		virtual Number	slope()			{return  7.0*UOM::mV; }
		virtual Number	xinfExponent()	{return  2; }

		// Time constants (see Brown et al. 1993).
		// Note that kinetics do not follow independent gating model
		// and that order 2 is for voltage sensitivity only.
		virtual Number	tauMax()		{return  1.2*UOM::msec; }	
		virtual Number	tauMin()		{return  0.4*UOM::msec; }

		// State vector label functions
		virtual const char* componentName()	{return "CaL"; }
		virtual const char** stateLabels()	{ 
			static const char* sl[] = { "s" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};
};


#endif // #ifndef __IONCHAN_CA_BAKER_2003_H_
