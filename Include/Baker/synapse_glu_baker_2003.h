// Synapse Dynamics for Glutamate Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: synapse_glu_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// glumate synaptic conductances within a CA3 model.
//
// Synaptic responses have temperature dependencies but these are
// not well quantified. For now such dependencies are omitted except
// for NMDAR channels, for which an estimated Q10 is available.
//
// Parameter values adapted from Destexhe, Mainen and Sejnowski are
// based on simplified models in the reference cited and have been
// converted for use in dual exponential formulation.
//
// Total synaptic weights for AMPA receptors is provided to allow
// GABAb responses (via GABAb triggered GIRKs coresident in synapses)
// to be affected by synaptic weight changes.
// 
// References:
//
// Allen C and Stevens CF (1994). An evaluation of causes for unreliability
// of synaptic transmission. PNAS USA 91, 10380-10383.
//
// Andrasfalvy BK and Magee JC (2001). Distance-dependent increase in AMPA 
// receptor number in the dendrites of adult hippocampal CA1 pyramidal
// neurons. J. Neuroscience 21, 9151-9159. 
//
// Benke TA Luthi A, Plamer MJ, Wikstrom MA, Anderson WW, Isaac  JTR, 
// Collingridge GL (2001). Mathematical modelling  of non-stationary 
// fluctuation analysis for studying channel properties of synaptic
// AMPA receptors. J. Physiology 537, 407-429.
//
// Debanne D, Guerineau NC, Gahwiler BH, Thompson SM (1995). Physiology
// and pharmacology of unitary synaptic connections between pairs of
// cells in areas CA3 and CA1 of rat hippocampal slice cultures.
// J. Neurophysiology 73, 1282-1294.
//
// Destexhe A, Mainen ZF, Sejnowski TJ (1998). Kinetic models of synaptic
// transmission, in Methods of Neuronal Modeling 2nd edition, ed. Kock C
// and Segev I. MIT Press.
//
// Flint AC, Maisch US, Weishaupt JH, Kriegstein AB, Moyner H (1997).
// NR2A subunit expression shortens NMDA receptor synaptic currents in
// developing neocortex. J. Neuroscience 17(7), 2469-2476.
//
// Grishin AA, Gee CE, Gerber U, Benquet P (2004). Differential
// calcium-dependent modulation of NMDA currents in CA1 and CA3
// hipppocampal pyramidal cells. J. Neuroscience 24(2), 350-355.
//
// Grishin AA, Benquet P, Gerber U (2005). Muscarinic receptor 
// stimulation reduces NMDA responses in CA3 hippocampal pyramidal 
// cells via Ca++ - dependent activation of tyrosine phosphatase. 
// Neuropharmacology (in press).
//
// Jahr, CW and Stevens, CF (1990). Voltage dependence of NMDA-activated
// macroscopic conductances predicted by single-channel kinetics.
// J. Neuroscience 10, 3178-3182.
//
// Kampa BM, Clements J, Jonas P, Stuart GJ (2004). Kinetics of Mg++ unblock
// of NMDA receptors: implications for spike-timing dependent synaptic
// plasticity. J. Physiology (London) 556, 337-345.
//
// Svoboda K, Tank DW, Denk W (1996). Direct meansurement of coupling between
// dendritic spines and shafts. Science 272, 716-719. 
//
// Williams SH and Johnston D (1991). Kinetic properties of two anatomically
// distinct excitatory synapses in hippocampal CA3 pyramidal neurons.
// J. Neurophysiol. 66(3), 1010-1020.


// Only include this header once
#ifndef __SYNAPSE_GLU_BAKER_2003_H_
#define __SYNAPSE_GLU_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation
namespace BAKER_2003 {

	// Class hierarchy
	class AMPA_SynapticResp;

	class NMDA_SynapticResp;
		class NMDA_NR2A_SynapticResp;
		class NMDA_NR2B_SynapticResp;

	class AC_Glu_SynapticResp;
	class PP_Glu_SynpaticResp;
	class MF_Glu_SynapticResp;


	// ================================================================
	// AMPA Synaptic Conductance Class (mostly abstract class)
	// ================================================================

	class AMPA_SynapticResp : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		AMPA_SynapticResp (Number gVal=0);
		virtual ~AMPA_SynapticResp ();

		// Parameters from Andrasfalvy and Magee, 2001.
		// See also Benke et al. for similar values.
		virtual SimTime		tau1() { return 0.6*UOM::msec; }
		virtual SimTime		tau2() { return 2.8*UOM::msec; }
		virtual Number		Vrev() { return _Vrev; }

		// Provide a nominal spine neck resistance (Svoboda et al.)
		virtual Number		spineNeckResistance() { return 100*UOM::megaohm; }

		// Provide a default release probability for AMPAR synapses.
		virtual Number		releaseProbability() { return 0.24f; }

		// State vector label functions (subclass should override)
		virtual const char* componentName() {return "AMPAR"; }

	protected:
		static const Number	_Vrev;					// reversal potential
		static DualExpSynapticCondClassCache _CC;	// class cache

		// Return address of the class cache
		DualExpSynapticCondClassCache* pDualExpClassCache() { return &_CC; }
	};

	// ----------------------------------------------------------------
	// NMDAVmTableEntry
	//
	// This public class stores precomputed values for the gate 
	// variable associated with the NMDAR channel Mg++ gate.
	//
	// Factors for computing Ca++ currents are not currently 
	// implemented. Note that this is a bit more complex than
	// GHK equations for calcium channels because Ca++ is not
	// the only ion with non-zero permissivity for the channel.
	// ----------------------------------------------------------------

	class NMDAVmTableEntry {
	public:
		Number			MgGate;			// Value of Mg gate at a voltage
	};

	// ================================================================
	// (Abstract) NMDA Synaptic Conductance Class
	// The following special properties apply:
	//
	// RESP:	
	//		1.	Provide relationship between voltage and Mg++ plug.
	//		2.	Adjust other parameter values for average voltage.
	//		3.	Set conductance based on ACh level.
	//
	// NOTES:	The state of the Mg++ plug affects ligand binding kinetics.
	//			The unbinding time constant tau2 is function of Vm 
	//			as shown in Kampa et al. figure 3. Otherwise the Mg++ 
	//			plug itself is treated as instantaneous.
	//
	//			Current for Ca++ ions is not currently implemented.
	//			Applying GHK equations involves the permability of
	//			different ions besides just Ca++.
	//
	//			Variable tau2 is also not supported owing to performance
	//			implications. A nominal value derived from mean membrane
	//			potential is used instead. Q10 factors are applied to 
	//			the derived tau2 value.
	//
	//			External Mg++ concentration is difficult to determine.
	//			Physiological values are 1-2 mM. 2mM is typical for slice 
	//			work but this is largely to prevent excess spiking. 
	//			Rats have somewhat lower Mg++ levels than humans.
	// ================================================================

	class NMDA_SynapticResp : public DualExpSynapticCond {

	public:

		// Constructors and destructor for abstract class
		NMDA_SynapticResp(Number gVal=0);
		virtual ~NMDA_SynapticResp();

		// Factors for Mg++ block behavior (Jahr & Stevens) -----------

		// These factors are used to load the VmTable and a separate table
		// must be provided in subclasses that provide different parameters.

		// Voltage sensitivity (analogous to zFRT in energy barrier models)
		virtual Number		MgVmMult()	{ return 0.062/UOM::mV; }

		// Factors associated with external Mg++ concentration.
		// [Mg++] at the receptor is assumed to be the same as in CSF. 
		// For rats this is in the range of 1 to 2 mM, though for some 
		// reason an accurate value does not seem to be readily available.
		// From human and dog studies, the value 1.3 mM is representative of 
		// CSF Mg++. Note that blood plasma values for [Mg++] are different.
		virtual Number		MgExtConc()	{ return 1.3*UOM::mM; } // concentration
		virtual Number		MgKdAtV0()	{ return 3.57*UOM::mM; }	// Mg++ Kd

		// Access the Mg++ gate value via table lookup
		inline  Number		MgGateValue()
		{	NMDAVmTableEntry*	ent = *pNMDAVmTable()+container()->VmIndex();

			return VTableInterp(container()->VmRem()/VStepForIndex,
				(ent-1)->MgGate,ent->MgGate, (ent+1)->MgGate,(ent+2)->MgGate);
		}

		// Factors for ACh modulation adjustments ---------------------

		inline  Number		AChLevel() { return _AChLevel; }
		virtual void		AChLevel(Number ach);

		// Return the current ACh conductance modulation
		inline  Number		AChMod() { return _AChMod; }

		// ACh parameters. Conductance is adjusted based on up to two
		// pathways interacting. The formula for conductance modulation is:
		// 
		// (1+a1*ach/(ach+kd1)) * (1+a2*ach/(ach+kd2))
		//
		// Values supplied here are defaults only and result in no modulation.

		virtual Number		AChA1()		{ return 0; }
		virtual Number		AChKd1()	{ return 1*UOM::microM; }
		virtual Number		AChA2()		{ return 0; }
		virtual Number		AChKd2()	{ return 1*UOM::microM; }

		// Apply an ACh neuromodulation rule to adjust F1.
		// id = token("AChModulator"), nv=1, values[0] = ACh concentration
		virtual void		setModParams(TokenId id, int nv, Number* values);

		// Factors and functions for electrical behavior --------------

		// Get conductance using Mg++ plug and membrane potential
		virtual Number		conductance();

		// Current using Mg++ plug and Ohm's law behavior
		virtual Number		Iion();

		// Total current reversal potential.
		virtual Number		Vrev() { return _Vrev; }

		// Factors and functions for variable tau2 --------------------

		// Voltage dependency between tau2 and Vm is based on a
		// linear fit for Kampa et al. fig 3 normalized
		// half-duration (nhd) values for [Mg++] = 1mM as in:
		//
		//	tau2 = tau2Max*max(nhdMin, nhdAtV0+ndhSlope*Vm)
		//
		// To simplify processing a nominal mean Vm is assumed
		// and that value is used to derive tau2.

		virtual Number		nhdAtV0()		{ return 0.74f; }
		virtual Number		nhdSlope()		{ return 6.4e-3/UOM::mV; }
		virtual Number		nhdMin()		{ return 0.2f; }

		// Return the membrane potential corresponding to tau2Max
		virtual Number		ratedVm() = 0;		// subclass responsibility

		// Return the nominal mean membrane potential for use in
		// computing the constant tau2 value used in dual exp form. 
		virtual Number		nominalVm() = 0;	// subclass responsibility

		// Return the tau2 value corresponding to the voltage
		// supplied in ratedVm(). The above formula ia used
		// to derive a tau2Max and from that tau2.
		virtual SimTime		ratedTau2() = 0;		// subclass responsibility

		// Return a time constant dependent on mean membrane potential.
		// (Consider caching this for a performance improvement).
		virtual	SimTime		tau2();

		// ODE Solver interface ---------------------------------------

		virtual void		simulationStarted();
		virtual void		simulationEnded();

	protected:

		// ACh modulation values
		Number				_AChLevel;				// current concentration
		Number				_AChMod;				// current modulation

		static const Number	_Vrev;					// reversal potential
		static NMDAVmTableEntry*	_NMDAVmTable;	// Vm dependent values table

		// Locate the NMDA Vm table (subclass may override if needed)
		virtual	NMDAVmTableEntry**	pNMDAVmTable() { return &_NMDAVmTable; }

		// Load the NMDA Vm table
		virtual void		loadNMDAVmTable();

		// Return the Mg++ plug state for a given voltage
		virtual Number		MgGateValueForTable(Number vm); 
	};

	// ================================================================
	// NMDA NR2A Synaptic Conductance Class
	// ================================================================

	class NR2A_SynapticResp : public NMDA_SynapticResp {

	public:

		// Constructors and destructor
		NR2A_SynapticResp(Number gVal=0.2*UOM::nanoS) : NMDA_SynapticResp(gVal) {}
		virtual ~NR2A_SynapticResp() {}

		// Parameter values from Flint et al. Tau1 is derived from Kampa et al.
		virtual Number		ratedTempC()	{ return 25; }
		virtual	Number		Q10()			{ return 3.0f; }
		virtual Number		ratedVm()		{ return -30*UOM::mV; }
		virtual Number		nominalVm()		{ return -60*UOM::mV; }
		virtual SimTime		tau1()			{ return 2*UOM::msec; }
		virtual SimTime		ratedTau2()		{ return 116*UOM::msec; }

		// ACh modulation parameters. See Grishin et al for
		// differences between CA3 and CA1. These values roughly
		// represent the finding in Grishin but are hardly the
		// only interpretation. Internal calcium concentrations are
		// not simulated but have been shown (Grishin) to be critical
		// in selecting up versus down regulation of NMDAR currents.
		virtual Number		AChA1()		{ return +0.39f; }
		virtual Number		AChKd1()	{ return 1*UOM::microM; }

		// Provide a nominal spine neck resistance (Svoboda et al.)
		virtual Number		spineNeckResistance() { return 100*UOM::megaohm; }

		// State vector label functions
		virtual const char* componentName()		{return "NR2A"; }

	protected:
		static DualExpSynapticCondClassCache _CC;

		// Return address of the class cache
		DualExpSynapticCondClassCache* pDualExpClassCache() { return &_CC; }
	};

	// ================================================================
	// NMDA NR2B Synaptic Conductance Class
	// ================================================================

	class NR2B_SynapticResp : public NMDA_SynapticResp {

	public:

		// Constructors and destructor
		NR2B_SynapticResp(Number gVal=0.2*UOM::nanoS) : NMDA_SynapticResp(gVal) {}
		virtual ~NR2B_SynapticResp() {}

		// Parameter values from Flint et al. Tau1 is derived from Kampa et al.
		virtual Number		ratedTempC()	{ return 25; }
		virtual	Number		Q10()			{ return 3.0f; }
		virtual Number		ratedVm()		{ return -30*UOM::mV; }
		virtual Number		nominalVm()		{ return -60*UOM::mV; }
		virtual SimTime		tau1()			{ return 2*UOM::msec; }
		virtual SimTime		ratedTau2()		{ return 256*UOM::msec; }

		// ACh modulation parameters. See Grishin et al for
		// differences between CA3 and CA1. These values roughly
		// represent the finding in Grishin but are hardly the
		// only interpretation. Internal calcium concentrations are
		// not simulated but have been shown (Grishin) to be critical
		// in selecting up versus down regulation of NMDAR currents.
		virtual Number		AChA1()		{ return +0.39f; }
		virtual Number		AChKd1()	{ return 1*UOM::microM; }

		// Provide a nominal spine neck resistance (Svoboda et al.)
		virtual Number		spineNeckResistance() { return 100*UOM::megaohm; }

		// State vector label functions
		virtual const char* componentName()		{return "NR2B"; }

	protected:
		static DualExpSynapticCondClassCache _CC;

		// Return address of the class cache
		DualExpSynapticCondClassCache* pDualExpClassCache() { return &_CC; }
	};

	// ================================================================
	// CA3 Associational Collateral Glutamate Synaptic Conductance
	// ================================================================

	class AC_Glu_SynapticResp : public SynapticGroupResponse {

	public:

		// Constructors and destructor
		AC_Glu_SynapticResp (
			Number gAMPAR=1.0*UOM::nanoS,		// AMPA conductance 
			Number gNMDAR=0.2*UOM::nanoS);		// NMDA conductance

		virtual ~AC_Glu_SynapticResp () {}

		// Accessors
		virtual AMPA_SynapticResp*		ampa() { return _ampa; }
		virtual NMDA_SynapticResp*		nmda() { return _nmda; }

		// State vector label functions
		virtual const char* componentName()		{return "AC_GluR"; }

	protected:
		AMPA_SynapticResp*		_ampa;
		NMDA_SynapticResp*		_nmda;
	};

	// ================================================================
	// CA3 Perforant Path Glutamate Synaptic Conductance
	// ================================================================

	class PP_Glu_SynapticResp : public SynapticGroupResponse {

	public:

		// Constructors and destructor
		PP_Glu_SynapticResp (
			Number gAMPAR=1.0*UOM::nanoS,		// AMPA conductance 
			Number gNMDAR=0.2*UOM::nanoS);		// NMDA conductance

		virtual ~PP_Glu_SynapticResp () {}

		// Accessors
		virtual AMPA_SynapticResp*		ampa() { return _ampa; }
		virtual NMDA_SynapticResp*		nmda() { return _nmda; }

		// State vector label functions
		virtual const char* componentName()		{return "PP_GluR"; }

	protected:
		AMPA_SynapticResp*		_ampa;
		NMDA_SynapticResp*		_nmda;
	};

	// ================================================================
	// CA3 Mossyy Fiber AMPA Synaptic Conductance Class.
	// ================================================================

	class MF_AMPA_SynapticResp : public DualExpSynapticCond {

	public:

		// Constructors and destructor
		MF_AMPA_SynapticResp (Number gVal = 1.0*UOM::nanoS);
		virtual ~MF_AMPA_SynapticResp ();

		// Parameters adapted from Williams and Johnston.
		// MF have similar unitary parameters though these
		// values may omit slow currents mentioned in Henze et al.
		virtual SimTime		tau1() { return 3.0*UOM::msec; }
		virtual SimTime		tau2() { return 3.0*UOM::msec; }
		virtual Number		Vrev() { return _Vrev; }

		// State vector label functions
		virtual const char* componentName() {return "MF_AMPAR"; }

	protected:
		static const Number	_Vrev;					// reversal potential
		static DualExpSynapticCondClassCache _CC;	// class cache

		// Return address of the class cache
		DualExpSynapticCondClassCache* pDualExpClassCache() { return &_CC; }
	};

	// ================================================================
	// CA3 Mossy Fiber Glutamate Synaptic Conductance
	// ================================================================

	class MF_Glu_SynapticResp : public SynapticGroupResponse {

	public:

		// Constructors and destructor
		MF_Glu_SynapticResp (
			Number gAMPAR=2.5*UOM::nanoS,		// AMPA conductance 
			Number gNMDAR=0.2*UOM::nanoS);		// NMDA conductance

		virtual ~MF_Glu_SynapticResp () {}

		// Accessors
		virtual MF_AMPA_SynapticResp*	ampa() { return _ampa; }
		virtual NMDA_SynapticResp*		nmda() { return _nmda; }

		// State vector label functions
		virtual const char* componentName()		{return "MF_GluR"; }

	protected:
		MF_AMPA_SynapticResp*	_ampa;
		NMDA_SynapticResp*		_nmda;
	};
};


#endif // #ifndef

