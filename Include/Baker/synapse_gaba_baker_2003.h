// Common Synapse Dynamics for GABAergic Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: synapse_gaba_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// GABA synapse conductances.
//
// Synaptic responses have temperature dependencies but these are
// not well quantified. For now such dependencies are omitted.
//
// Parameter values adapted from Destexhe, Mainen and Sejnowski are
// based on simplified models in the reference cited and have been
// converted for use in dual exponential formulation.
//
// GABA-b response is based on Otis et al. Paired pulse
// facilitation is used to induce frequency dependency in GABAb
// rather than the more commonly used G-Protein pooling model.
//
// Because both GABA-a, GABA-as, and GABA-b use a form of paired 
// pulse rule for presynaptic plasticity and because of the synaptic
// delay for GABA-b, GABA-b cannot be grouped with other synapses.
// This is not so crazy anyway since GABA-b receptors do not 
// coreside with GABAa synapses in most cases. We now know that the
// K+ component can be resident in synaptic spines (see Huang et al.).
// 
// A strange object model is used for GABA-b responses and plasticity.
// Even so, GABA-b is triggered by GABA synapses and in this model is 
// not directly associated with AMPAR-NMDAR synapses. The relationship 
// is established by GABA-b postsynaptic plasticity whereby GABA-b 
// weights are adjusted based on AMPAR weights, which in turn are 
// derived from NMDAR-dependent plasticity.
// 
// References:
//
// Banks MI, Li T-B, Pearce RA (1998). The synaptic basis of GABA-A,slow.
// J. Neuroscience 18(4), 1305-1317.
//
// Destexhe A, Mainen ZF, Sejnowski TJ (1998). Kinetic models of synaptic
// transmission, in Methods of Neuronal Modeling 2nd edition, ed. Kock C
// and Segev I. MIT Press.
//
// Huang CS, Shi S-H, Ule J, Ruggiu M, Barker LA, Darnell RB, Jan YN,
// Jan YL (2005). Common molecular pathways mediate long-term potentiation
// of synaptic excitation and slow synaptic inhibition. Cell 123, 105-118.
//
// Miles R, Toth K, Gulyas AI, Hajos N, Freund TF (1996). Differences between
// somatic and dendritic inhibition in the hippocampus. Neuron 16, 815-823.
//
// Otis TS, De Koninck Y,  Mody I (1993). Characterization of synaptically
// elicited GABA_b responses using patch-clamp recordings in rat hippocampal
// slices. J. Physiology 463, 391-407.
//
// Scanziani M (2000). GABA Spillover activates postsynaptic GABAb receptors
// to control rhythmic hippocampal activity. Neuron 25, 673-681.


// Only include this header once
#ifndef __SYNAPSE_GABA_BAKER_2003_H_
#define __SYNAPSE_GABA_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation
namespace BAKER_2003 {

	// Class hierarchy
	class GABAa_SynapticResp;
	class GABAb_SynapticResp;
	class GABAab_SynapticResp;

	
	// ================================================================
	// GABA(a) Synaptic Response Class
	// ================================================================

	class GABAa_SynapticResp : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		GABAa_SynapticResp(Number gVal=0.2*UOM::nanoS);
		virtual ~GABAa_SynapticResp();

		// Parameter values from Destexhe, Mainen & Sejnowski
		// tau1 is chosen to give peak at 1 msec after release.
		virtual SimTime		tau1()		{ return 0.3*UOM::msec; }
		virtual SimTime		tau2()		{ return 5.6*UOM::msec; }
		virtual Number		Vrev()		{ return _Vrev; }

		// State vector label functions
		virtual const char* componentName()	{ return "GABAa"; }

	protected:
		static const Number	_Vrev;					// reversal potential
		static DualExpSynapticCondClassCache _CC;	// class cache

		// Return address of the class cache
		DualExpSynapticCondClassCache* pDualExpClassCache() { return &_CC; }
	};
	
	// ================================================================
	// GABA(a-slow) Synaptic Response Class
	// ================================================================

	class GABAas_SynapticResp : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		GABAas_SynapticResp(Number gVal=0.2*UOM::nanoS);
		virtual ~GABAas_SynapticResp();

		// Parameter values from Banks et al. at 35-degC.
		// Vrev is assumed to be the same as Gaba-a (normal Cl-).
		virtual SimTime		tau1()		{ return 5*UOM::msec; }
		virtual SimTime		tau2()		{ return 40*UOM::msec; }
		virtual Number		Vrev()		{ return _Vrev; }

		// State vector label functions
		virtual const char* componentName()	{ return "GABAas"; }

	protected:
		static const Number	_Vrev;					// reversal potential
		static DualExpSynapticCondClassCache _CC;	// class cache

		// Return address of the class cache
		DualExpSynapticCondClassCache* pDualExpClassCache() { return &_CC; }
	};
	
	// ================================================================
	// GABA-b Synaptic Response Class.
	//
	// This is an approximate fit to Otis et al. A triple exponent
	// form allows for efficiency in processing a group of GABAb 
	// synapses as a unit while still retaining a reasonable fit to
	// the data. Because GABA-b receptors bind with extracellular GABA,
	// there may not be a simple one-to-one relationship between
	// afferent GABA-a events and GABA-b stimulation
	// ================================================================

	class GABAb_SynapticResp : public TripleExpSynapticCond {

	public:

		// Constructors and destructor
		GABAb_SynapticResp(Number gVal=0.4*UOM::picoS);
		virtual ~GABAb_SynapticResp();

		// Parameters approximately fit single pulse response in
		// Otis et al. at 34 degC. Additional delay is used to
		// compensate for lack of a fourth order activation process.
		virtual SimTime		tau1()		{ return 70*UOM::msec; }
		virtual SimTime		tau2()		{ return 110*UOM::msec; }
		virtual SimTime		tau3()		{ return 516*UOM::msec; }
		virtual Number		c2()		{ return 0.84f; }		
		virtual Number		Vrev()		{ return _Vrev; }

		// Lag time for G-protein activation and/or GABA diffusion 
		// at 34-36 degC. This delay is only well documented for 
		// single spikes and may not be relevant in all settings.
		virtual SimTime		synapticDelay() { return 20*UOM::msec; }

		// Return the conductance adjusted by any postsynaptic rule
		virtual Number		conductance();

		// State vector label functions
		virtual const char* componentName()	{ return "GABAb"; }

	protected:
		static const Number	_Vrev;							// reversal potential
		static TripleExpSynapticCondClassCache	_CC;		// class cache

		// Return address of the class cache
		TripleExpSynapticCondClassCache* pTripleExpClassCache() { return &_CC; }
	};
};	


#endif // #ifndef

