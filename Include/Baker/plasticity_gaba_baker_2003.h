// Synaptic Plasticity for GABAergic Synapses
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
// GABA synapse plasticity. Two version were considered: one from
// data involving regular spike trains (Scanziani) and the
// other from single pulse pairings (Otis et al.). Of these,
// the spike train data is probably more representative
// though neither should be considered a complete model.
//
// Seeger et al. show that GABA IPSC is modulated by ACh.
// The model here is a fit to that data and the effect is 
// assumed to be presynaptic and thus release probability 
// related. There is no hard data on GABA-b effects, but
// the relationship is assumed to be the same as GABA-a
// because GABA-b receptors are activated by release spillover.
// 
// Results from Scanziani are based on charge transfer and
// as such include any frequency dependencies inherent in
// GABA-b responses.
//
// GABA-b postsynaptic plasticity is derived from Huang et al. but
// considerably simplified. GABA-b synaptic weight is adjusted based
// on the total synaptic weight for AMPAR synapes. To find the related
// synapses, multiple AMPAR response types are allowed (for AC and PP
// synapses in the same compartment). Parameter values are obviously
// speculative and there is no attempt to match the time evolution of
// changes in the GABA-b current. 
//
// While counter-intuitive, the AMPAR adjusted weight for GABA-b is 
// derived from the average AMPAR weight. This is based on an assumption
// that AMPAR and GABA synapses occur in roughly proportional quantities.
// As a matter of implementation, the AMPAR adjusted weight is
// implemented as a change in GABA-b conductance so that it applies to all
// all synapses in the compartment and changes together with AMPAR weights.
//
// References:
//
// Huang CS, Shi S-H, Ule J, Ruggiu M, Barker LA, Darnell RB, Jan YN,
// Jan YL (2005). Common molecular pathways mediate long-term potentiation
// of synaptic excitation and slow synaptic inhibition. Cell 123, 105-118.
//
// Otis TS, De Koninck Y,  Mody I (1993). Characterization of synaptically
// elicited GABA_b responses using patch-clamp recordings in rat hippocampal
// slices. J. Physiology 463, 391-407.
//
// Scanziani M (2000). GABA spillover activates postynaptic GABAb
// receptors to control rhythmic hippocampal activity.
// Neuron 25, 673-681.
//
// Seeger T, Fedorova I, Zheng F, Miyakawa T, Koustova E, Gomeza J,
// Basile AS, Alzheimer C, Wess J (2004). M2 Muscarinic acetylcholine
// receptor knock-out mice show deficits in behavioral flexibility,
// working memory, and hippocampal plasticity. J Neurosci 24, 10117-10127.


// Only include this header once
#ifndef __PLASTICITY_GABA_BAKER_2003_H_
#define __PLASTICITY_GABA_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation
namespace BAKER_2003 {
	
	// ================================================================
	// GABA Presynaptic Rule Class for paired pulse depression
	// using frequency dependent results in Scanziani (see figure 2B2).
	// ================================================================

	class GABA_PresynapticRule : public PlasticityRule {
	public:

		// Constructors and destructor
		GABA_PresynapticRule() {}
		virtual ~GABA_PresynapticRule() {}

		// Include a paired pulse facilitation/depression using:
		//
		//		pp = 1+a0*exp(-isi/tau)
		//
		// pp is to set to the event quantity, which can be interpreted
		// as either a release probability or a continuous quantity.

		// ACh modulation is determined using a Michaelis-Menten formulation
		// from which a multiplier to the quantity is derived.

		// Access the current ACh concentration level
		inline  Number			AChLevel()		{ return _AChLevel; }
		virtual void			AChLevel(Number ach) { _AChLevel = ach; }

		// Apply an ACh neuromodulation rule to adjust F1.
		// id = token("AChModulator"), nv=1, values[0] = ACh concentration
		virtual void			setModParams(TokenId id, int nv, Number* values);

		// Parameter accessors (subclass resp) ------------------------

		virtual SimTime			tau() = 0;	// Time constant for PPF/PPD
		virtual Number			a0() = 0;	// Maximum change (for isi=0)

		virtual Number			AChA()=0;	// Maximum marginal ACh effect
		virtual Number			AChKd()=0;	// Half activation level of ACh 

		virtual bool			useRandomRelease() { return false; }

		// Framework interfaces ---------------------------------------

		// Set event quantity based on rule
		virtual void			finalizeAPEvent(ActionPotentialEvent* apEvent);

		// There is no plasticity state data for this object
		virtual unsigned int	plasticityStateSize() { return 0; }

		// Return an id for this type of plasticity
		virtual TokenId			plasticityTypeId() {return _GABA_PresynapticId; }

	protected:
		Number					_AChLevel;				// ACh concentration

		// Provide a static token id for this type of plasticity
		static const TokenId	_GABA_PresynapticId;	// type id of this rule
	};	

	// ----------------------------------------------------------------
	// GABAa Presynaptic Rule Class for paired pulse depression
	// ----------------------------------------------------------------

	class GABAa_PresynapticRule : public GABA_PresynapticRule {
	public:

		// Constructors and destructor
		GABAa_PresynapticRule() {}
		virtual ~GABAa_PresynapticRule() {}

		// Parameter values (GABAa is effectively PPD)
		virtual SimTime			tau()	{ return 15*UOM::msec; }
		virtual Number			a0()	{ return -0.7f; }

		// Fit ACh modulation to values in Seeger et al.
		virtual Number			AChA()	{ return -0.80f; }
		virtual	Number			AChKd()	{ return 3.3*UOM::microM; }

		// Use binary random release to set final quantity
		virtual bool			useRandomRelease() { return true; }

		// Return a component name for reporting
		virtual const char*		componentName() {return "GABAa_PresynapticRule"; }
	};

	// ----------------------------------------------------------------
	// GABAb Presynaptic Rule Class for paired pulse depression
	// ----------------------------------------------------------------

	class GABAb_PresynapticRule : public GABA_PresynapticRule {
	public:

		// Constructors and destructor
		GABAb_PresynapticRule() {}
		virtual ~GABAb_PresynapticRule() {}

		// Parameter values (GABAb is effectively PPF)
		virtual SimTime			tau()	{ return 33*UOM::msec; }
		virtual Number			a0()	{ return 3.7f; }

		// For lack of better values, use GABAa ACh parameters
		// since GABAb is activated by GABA spillover. In reality the
		// relationship is likely to be more non-linear than assumed here.
		virtual Number			AChA()	{ return -0.80f; }
		virtual	Number			AChKd()	{ return 3.3*UOM::microM; }

		// Return a component name for reporting
		virtual const char*		componentName() {return "GABAb_PresynapticRule"; }
	};

	// ----------------------------------------------------------------
	// GABAb Postsynaptic Rule Class for long-term plasticity
	// ----------------------------------------------------------------

	class GABAb_PostsynapticRule : public PlasticityRule {
	public:

		// Constructors and destructor
		GABAb_PostsynapticRule(
			Number baseWght=1,				// use current gGABAb
			Number amparInc=0,				// no AMPAR weight effect
			SimTime tau=60*UOM::sec);		// 1 min weight time constant

		virtual ~GABAb_PostsynapticRule();

		// Accessors
		inline  Number			baseWeight() { return _baseWeight; }
		virtual void			baseWeight(Number w) { _baseWeight = w; }

		inline  Number			amparIncrement() { return _amparIncrement; }
		virtual void			amparIncrement(Number dw) { _amparIncrement=dw; }

		inline  SimTime			tauW() { return _tauW; }
		virtual void			tauW(SimTime tau) { _tauW = tau; }

		inline  Number			weight() { return _weight; }

		// Return the size to allocate for plasticity state data.
		virtual unsigned int	plasticityStateSize() { return 0; };

		// Take action at the end of the time step (before AP purge)
		virtual void			applyEndOfStep(ActionPotentialEventQueue& apQueue);

		// Return a component name for reporting
		virtual const char*		componentName() {return "GABAb_PostsynapticRule"; }

	protected:

		// Associated AMPAR with NMDAR dependent plasticity
		SynapticResponseVector	_ampar;

		// Current adjusted weight for use in computing condutance
		Number					_weight;

		// Synapse weight when there are no AMPAR associated increments
		Number					_baseWeight;

		// GABA-b weight increment per increment of average AMPAR weight.
		// Average AMPAR weight is used to scale with different comparment
		// sizes on the assumptions that GABA and AMPA synapses are roughly
		// in the same ratio regardless of compartment size.
		Number					_amparIncrement;

		// Flag indicating that _ampar has been initialized (even if empty)
		bool					_amparInitialized;

		// Time constant by which actual weight and target weight converge.
		// This is assumed to be large relative to the time step size.
		SimTime					_tauW;

		// Locate associated AMPAR responses and leave in _ampar
		virtual void			locateAMPAR();
	};

};	


#endif // #ifndef

