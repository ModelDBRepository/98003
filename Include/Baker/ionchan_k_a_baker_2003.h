// K-A Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_a_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K-A channel based on data from the Hoffman 1997 article.
//
// Parameter values are taken from Hoffman. An exponent order of 4
// is used based on data fitting and the observation that currents near
// rest must be small. Bekkers found kinetic fit for order 3
// while other modelers (e.g. Borg-Graham) have used order 4 based on 
// similarity with the squid axon K-Dr and a tetrameric structure.
//
// Values for activation time constants are estimates. Hoffman reports
// an activation time constant of 1 msec for voltages used. Bekkers
// found a much higher tau value in layer 5 cells when more depolarized.
// In younger animals (11-16d rat CA1), Martina finds two transient 
// currents, one partially sensitive to TEA with a larger time constant
// (~6-8 msec) and the other insensitive to TEA with a smaller time constant
// (~1 msec), both of which are largely voltage insensitive. Klee et al.
// found different time constants in CA1 vs. CA3 cells.
//
// The xinfExponent form of multiple gates is used for activation to allow
// the resulting time constant to be independent of the order of the gate
// and to have the same value for activation and deactivation. Martina et al.
// showed rapid deactivation but not sufficiently rapid to justify the
// factor of 4 change in time constants implied by fourth-order gating.
//
// Inactivation time constants are based on apparent minimum values of
// h tau at different voltages as measured by Hoffman et al. and 
// Pan & Colbert. A linoid function is used to provide a smooth curve
// for the inactivation time constant as a function of voltage.
//
// Q10 values here are estimates using typical values for voltage-gated
// channels. Experimental values were derived at near physiological
// temperature so this should have minimal effect. The Q10 value for
// inactivation is derived using the fact that Marina et al. results
// were at room temperature while those of Hoffman et al. were at,
// 35 degC. The accuracy of this estimate is probably limited.
// 
// Note that CA1 pyramidal cells express Kv4.2 but not Kv4.3 channels 
// while interneurons as well as CA3 and DG pyramidal cells express both.
//
// References:
// 
// Bekkers JM, 2000. Properties of voltage-gated potassium currents in
// nucleated patches from large laryer 5 cortical pyramidal neurons
// of the rat. J. Physiology 525.3, 593-609.
//
// Borg-Graham LJ, 1998. Interpretations of Data and Mechanisms for 
// Hippocampal Cell Models, in Cerebral Cortex vol 13. New York: Plenum Press.
// Also available via Surf-Hippo web site.
//
// Hoffman DA, Magee JC, Colbert CM, Johnston D, 1997.
// K+ channel regulation of signal propagation in dendrites
// of hippocampal cells. Nature 387(6636), 869-875.
//
// Klee R, Ficker E, and Heinemann U, 1995. Comparison of voltage-
// dependent postassium currents in rat pyramidal neurons acutely
// isolated from hippocampal regions CA1 and CA3. 
// J. Neurophysiol. 74(5), 1982-1995.
//
// Pan E and Colbert CM, 2001. Subthreshold incactivation of Na+ and
// K+ channels supports activity-dependent enhancement of back-propagating
// action potentials in hippocampal CA1. J. Neurophysiology 85, 1013-1016.
//
// Martina M, Schultz JH, Ehmke H, Monyer H, Jonas P 1998. Functional
// and molecular differences between voltage-gated K+ channels of
// fast-spiking interneurons and pyramidal neurons of rat hippocampus.
// J. Neurosci. 18(20), 8111-8125.

 
// Only include this header once
#ifndef __IONCHAN_K_A_BAKER_2003_H_
#define __IONCHAN_K_A_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace BAKER_2003 {

	// a gate variable (proximal)

	class Proximal_K_A_a_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		Proximal_K_A_a_gate() {}
		virtual ~Proximal_K_A_a_gate() {}

		// Parameter values from Hoffman
		virtual Number	ratedTempC()	{return  35; }
		virtual Number	Q10()			{return	 4; }

		// Order 4 gate for vhalf=11. k=18 (per Hoffman) 
		virtual Number	Vhalf()			{return	-27.1*UOM::mV; }
		virtual Number	slope()			{return  22.9*UOM::mV; }
		virtual Number	xinfExponent()	{return  4; }

		// Time constants
		virtual Number	tauMax()		{return  1.0*UOM::msec; }	
		virtual Number	tauMin()		{return  1.0*UOM::msec; }

		// State vector label functions
		virtual const char*	componentName() {return "K_A"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "a" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// a gate variable (distal)

	class Distal_K_A_a_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		Distal_K_A_a_gate() {}
		virtual ~Distal_K_A_a_gate() {}

		// Parameter values from Hoffman
		virtual Number	ratedTempC()	{return  35; }
		virtual Number	Q10()			{return	 4; }

		// Order 4 gate for vhalf=-1, k=15 (per Hoffman)
		virtual Number	Vhalf()			{return	-32.8*UOM::mV; }
		virtual Number	slope()			{return  19.1*UOM::mV; }
		virtual Number	xinfExponent()	{return  4; }

		// Time constants
		virtual Number	tauMax()		{return  1.0*UOM::msec; }	
		virtual Number	tauMin()		{return  1.0*UOM::msec; }

		// State vector label functions
		virtual const char*	componentName() {return "K_A"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "a" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// a gate variable (blended proximal-distal response)

	class Blended_K_A_a_gate : public BlendedIonGate {
	
	public:
		// Constructors and destructor
		Blended_K_A_a_gate(Number bratio=0.5);
		virtual ~Blended_K_A_a_gate() {}

		// State vector label functions
		virtual const char*	componentName() {return "K_A"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "a" }; return sl; }
	};

	// b gate variable

	class K_A_b_gate : public EnergyBarrierTabGate {

	public:

		// Constructors and destructor
		K_A_b_gate() {}
		virtual ~K_A_b_gate() {}

		// Values derived by Hoffman + Pan&Colbert + Martina et al.
		virtual Number	ratedTempC()	{return  35; }
		virtual Number	Q10()			{return	 1.3f; }

		virtual Number	Vhalf()			{return -56*UOM::mV; }
		virtual Number	slope()			{return -8*UOM::mV; }

		// Special tau computation
		virtual Number	tauForTable(Number v);
		virtual Number	tauMin()		{return  5*UOM::msec; }

		// State vector label functions
		virtual const char*	componentName() {return "K_A"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "b" }; return sl; }

	protected:
		virtual AlphaBetaEntry** pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};

	// K-A channel for channels on proximal dendrite and soma

	class Proximal_K_A_channel : public M1HIonChannel {

	public:
		// Constructors and destructor
		Proximal_K_A_channel(Number gSpVal=0);
		virtual ~Proximal_K_A_channel() {}

		// Reversal potential for K+
		virtual Number Vrev() { return _Vrev; }
		
	protected:
		static const Number		_Vrev;
	};

	// K-A channel for channels on distal dendrite

	class Distal_K_A_channel : public M1HIonChannel {
	
	public:
		// Constructors and destructor
		Distal_K_A_channel(Number gSpVal=0);
		virtual ~Distal_K_A_channel() {}

		// Reversal potential for K+
		virtual Number Vrev() { return _Vrev; }
		
	protected:
		static const Number		_Vrev;
	};

	// Blended K-A channel for intermediate locations

	class Blended_K_A_channel : public M1HIonChannel {
	
	public:
		// Constructors and destructor
		Blended_K_A_channel(Number gSpVal=0, Number bratio=0.5);
		virtual ~Blended_K_A_channel() {}

		// Reversal potential for K+
		virtual Number Vrev() { return _Vrev; }
		
	protected:
		static const Number		_Vrev;
	};
};


#endif // #ifndef __IONCHAN_K_A_BAKER_2003_H_
