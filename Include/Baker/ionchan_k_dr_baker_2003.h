// K-Dr Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_dr_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file contains the classes used to implement the K-DR channel
// based on data from Hoffman et al. Only the activation component is
// modeled because inactivation is relatively slow and no data 
// specifically for dendritic channels is available.
//
// V-half values are not consistent between different experimenters.
// Hoffman values for both K-DR and K-A are used for consistency.
// Note that axons express a different form of K+ channels that 
// may well have different characteristics from those found in the 
// main cell body and dendrites.
//
// Order 1 gating may not be the common model but this improves
// repolarization, especially when K-A and K-C are reduced by ACh.
// Order 4 parameters voltage params are left as a reminder. If
// order 4 gating is used here, K+ conductances must be increased.
//
// Q10 values are an estimate based on other voltage-gated kinetics.
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
// Klee R, Ficker E, Heinemann U, Comparison of voltage-dependent potassium
// currents in rat pyramidal neurons acutely isolated from hippocampal
// regions CA1 and CA3.
//
// Martina M, Schultz JH, Ehmke H, Monyer H, Jonas P 1998. Functional
// and molecular differences between voltage-gated K+ channels of
// fast-spiking interneurons and pyramidal neurons of rat hippocampus.
// J. Neurosci. 18(20), 8111-8125.
//
// Traub R, Jefferys JG, MIles R, Whittington MA, and Toth K (1994).
// A branching dendritic model of a rodent CA3 pyramidal neurone.
// J. Physiol. 481(1), 79-95.

 
// Only include this header once
#ifndef __IONCHAN_K_DR_BAKER_2003_H_
#define __IONCHAN_K_DR_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace BAKER_2003 {

	// K-DR channel (generic delayed rectifier)

	class K_DR_channel : public Order1EnergyBarrierTabChannel {
	
	public:
		// Constructors and destructor
		K_DR_channel(Number gSpVal=0); 
		virtual ~K_DR_channel() {}

		// Temperature sensitivity
		virtual Number	ratedTempC()	{return  22; }
		virtual Number	Q10()			{return  4; }	// assumed Q10

		// Voltage response per Hoffman et al. (order 1 gating)
		virtual Number	Vhalf()			{return	 13*UOM::mV; }
		virtual Number	slope()			{return  11*UOM::mV; }

		// Voltage response per Hoffman et al. (order 4 gating)
		// virtual Number	Vhalf()			{return	-10.3*UOM::mV; }
		// virtual Number	slope()			{return  14.0*UOM::mV; }

		// Estimated time constants (see Martina et al. fig 4E)
		virtual Number	tauMax()		{return  8.0*UOM::msec; }	
		virtual Number	tauMin()		{return  1.0*UOM::msec; }
		virtual Number	gamma()			{return  0.8f; }

		// State vector label functions
		virtual const char*	componentName() {return "K_DR"; }
		virtual const char** stateLabels()	{
			static const char* sl[] = {"n"}; return sl; }

		// Reversal potential for K+
		virtual Number Vrev()			{return _Vrev; }
		
	protected:
		static const Number			_Vrev;
		virtual AlphaBetaEntry**	pAlphaBetaTable() { return &_abTable; }

	private:
		static AlphaBetaEntry*		_abTable;
	};
};

#endif // #ifndef __IONCHAN_K_DR_BAKER_2003_H_
