// Mammalian CA3 cell model based on reconstructed cell morphology
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// File: neuron_baker_2003.cpp
//
// Description:
//
// This header file contains the declarations used to implement
// a mammalian CA3 cell model based on reconstructed cell morphology.
//
// The morphology table must have one entry for the entire soma.
// All other entries must be for dendrites. An initial segment
// and axon compartments are added based on hard-coded parameters
// during the initialization process.
//
// Only cell L56a is used in the current simulations. Other cell
// morphologies are available for future use (not included here).
//
// References (see individual channels for further references):
//
// Alroy G. Hailing S, Yaari Y (1999) Protein kinase C mediates
// muscarinic block of intrinsic bursting in rat hippocampal neurons.
// J Physiol. 518, 71-79. 
//
// Cantrell AR, Ma JY, Scheuer T, Catterall WA (1996). Muscarinic
// modulation of sodium current by activation of protein kinase C
// in rat hippocampal neurons. Neuron 16, 1019-1026.
//
// Colbert CM and Pan E (2002). Ion channel properties underlying axonal 
// action potential initiation in pyramidal neurons.
// Nature Neuroscience 5(6), 533-538.
//
// Fishan A, Yamada M, Duttaroy A, Gan J-W, Deng C-X, McBain CJ,
// Wess J (2002). Muscarinic Induction of hippocampal gamma
// oscillations requires coupling of the M1 receptor to two
// mixed cation currents. Neuron 33, 615-624.
//
// Fisher RE, Gray R, and Johnston D. (1990). Properties and
// distribution of single voltage-gated calcium channels in
// adult hippocampal neurons. J. Neurophysiology 64, 91-104.
//
// Fischer RE and Johnston D (1990). Differential modulation of single
// voltage-gated calcium channels by cholinergic and adrenergic
// agonists in adult hippocampal neurons. J Neurophysiol. 64, 1291-1302.
//
// Frick A, Magee J, Koester HJ, Migliore M, and Johnston D. (2003).
// Normalization of Ca++ signals by small oblique dendrites of CA1
// pyramidal neurons. J. Neuroscience 23(8), 3243-3250.
//
// Hoffman DA, Magee JC, Colbert CM, Johnston D, 1997.
// K+ channel regulation of signal propagation in dendrites
// of hippocampal cells. Nature 387(6636), 869-875.
//
// Jaffe DB, Ross WN, Lisman JE, Lasser-Ross N, Miyakawa H, and Johnston D.
// 1994. A model for dendritic Ca++ accumulation in hippocampal pyramidal
// neurons based on fourescence imaging measurements. 
// J. Neurophysiology 71, 1065-1077.
//
// Kavalali, ET, Zhuo M, Bito  H, and Tsien RW (1997).
// Dendritic Ca++ channels characterized by recordings from
// isolated hippocampal dendritic segments. Neuron 18, 651-663.
//
// Lazarewica MT, Migliore M, Ascoli GA, 2002. A new bursting model of CA3
// pyramidal cell physiology suggests multiple locations for spike initiation.
// BioSystems 67, 129-137.
//
// Magee JC, 1998. Dendritic hyperpolarization-activated currents modify
// the integrative properties of hippocampal CA1 pyramidal neurons.
// J. Neuroscience 18(19), 7613-7624.
//
// Major G, Larkman AU, Jonas P, Sakmann B, Jack JJB 1994. Detailed
// passive cable models of whole-cell recorded CA3 pyramidal neurons in
// rat hippocampal slices. J. Neurosci. 14(8), 4613-4638.
//
// Menschik ED 1999. Cholinergic neuromodulation in hippocampal
// function and disease. University of Pennsylvania PhD dissertation.
//
// Migliore M, Cook EP, Jaffe DB, Turner DA, and Johnston D, 1995.
// Computer simulations of morphologically reconstructed CA3
// hippocampal neurons. J. Neurophysiol. 73(3), 1157-1168.
//
// Poolos NP and Johnston D, 1999. Calcium-activated potassium
// conductances contribute to action potential repolarization at the
// soma but not the dendrites of hippocampal CA1 pyramidal neurons.
// J. Neurosci. 19(13), 5205-5212.
//
// Scanziani M (2000). GABA Spillover activates postsynaptic GABAb receptors
// to control rhythmic hippocampal activity. Neuron 25, 673-681.
//
// Tsubokawa H and Ross WN (1997). Muscarinic modulation of spike back-
// propagation in the apical dendrites of hippocampal CA1 pyramidal neurons.
// J Neurosci. 17, 5782-5791. 
 

// Only include this header once
#ifndef __NEURON_BAKER_2003_H_
#define __NEURON_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace BAKER_2003 {

	// ----------------------------------------------------------------
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ----------------------------------------------------------------

	class PyramidalCell;
		class CA3PyramidalCell;

	// -------------------------------------------
	// Pryamidal cell neuron
	// -------------------------------------------

	class PyramidalCell : public MorphologicalNeuron {

	public:

		// Constructors and destructor (see initialize below)
		PyramidalCell(Model* m=NULL);
		virtual ~PyramidalCell();

		// Accessors for controlling ACh modulation
		inline  Number			AChLevel() { return _AChLevel; }
		virtual void			AChLevel(Number AChExtConc);

		// Accessors for controlling Na s-gate status.
		// S-gate disable is set by the most recent call to either
		// NaSGateDisabled(b) or as a side-effect of AChLevel(x).
		inline  bool			NaSGateDisabled() { return _NaSGateDisabled; }
		virtual void			NaSGateDisabled(bool disableSGate);

		inline  Number			NaSGateDisabledACh() { return _NaSGateDisabledACh; }
		inline  Number			NaSGateDisabledValue() { return _NaSGateDisabledValue; }

		// Knockout NMDA receptor effects selectively
		virtual void			knockoutNMDAR(
			bool					disableCurrent = true,
			bool					disablePlasticity = true,
			bool					disableAChMod = false,
			bool					disableCaDepSupp = false);
		
		// Special compartment accessors (others are inherited)
		virtual Compartment*	initialSegment() { return _compartments[1]; }

		// Accessors for selected parameters (those used in other objects)
		inline  Number			orientationX() { return _orientationX; }
		inline  Number			orientationY() { return _orientationY; }
		inline  Number			orientationZ() { return _orientationZ; }

		// Other framework interfaces ---------------------------------

		// Voltage and deriv to use for determining whether firing occurred
		virtual Number			firingVoltage() { return axonComp(_numAxonComp)->Vm(); }
		virtual Number			firingDeriv() { return axonComp(_numAxonComp)->VmDot(); }

		// Default ODESolver to use when neuron owns the model.
		virtual ODESolver*		defaultSolver();

	protected:

		// Variables for defining current (static) cell state ---------

		// Nominal external [ACh]. This will typically only have one of
		// two values, 100 microM or 0 microM reflecting ACh present or not.
		Number					_AChLevel;

		// Flag indicating whether Na channel s-gate is to be used.
		// Note that Na s-gate is automatically disabled when ACh is present.
		bool					_NaSGateDisabled;

		// Flag indicating whether ACh modulation of Inmdar is applied or not
		bool					_AChInmdarModDisabled;

		// Electrophysiology parameters -------------------------------

		// Factor to increase compartment membrane area to allow for spines
		// and other irregularities of structure with respect to a perfect
		// cylinder. This adjustment is applicable to dendrites only.
		Number					_areaAdjustment;

		// Electrophysiology parameters across the cell
		Number					_Rm;				// Effective membrane resistance
		Number					_Ri;				// Axial resistance
		Number					_Cm;				// Membrane capacitance

		Number					_RmAxon;			// Rm for axon only
		Number					_RiAxon;			// Ri for axon only
		Number					_CmAxon	;			// Cm for axon only
		Number					_VleakAxon;			// Leak reversal for axon only

		// Parameters for a leak current attributable to K+ etc. (non-ACh case).
		Number					_RmLeak;			// Leak resistance
		Number					_Vleak;				// Leak reversal

		// Initial value of Vm at soma (+IS+Axon) and at dendrites.
		Number					_Vinit;					// Vm initial value at soma
		Number					_VinitDendriteSlope;	// dendrite distance dependency
		Number					_VinitAxonSlope;		// axon distance dependency

		// Parameters covering IS and axon geometry.
		// Superclass variable _numAxonComp holds the number of axon compartments
		Number					_initSegLen;		// IS compartment length
		Number					_axonCompLen;		// Axon compartment length

		Number					_initSegRadius;			// IS radius
		Number					_axonProximalRadius;	// Axon radius at IS
		Number					_axonDistalRadius;		// Axon radius away from cell

		// Calcium buffer parameters ----------------------------------

		Number					_CaXrest;			// Resting [Ca++]
		Number					_CaXinit;			// Initial value for [Ca++]
		Number					_CaUbr;				// Unbound Ca++ ion ratio
		Number					_CaUbrOblique;		// UBR for thin oblique dendrites
		Number					_CaKd;				// MM half activation for pump
		Number					_CaVmax;			// MM peak pump rate
		Number					_mdShellDepth;		// micro-domain subshell depth
		Number					_somaShellDepth;	// allowance for nucleus in soma

		// Ion channel conductances -----------------------------------
		// These are under physiological conditions (e.g. 37-38 deg C)

		Number					_gNaTSoma;			// Na density for soma	
		Number					_gNaTProximal;		// Na density for proximal dendrites	
		Number					_gNaTDistal;		// Na density for distal dendrites	
		Number					_gNaPSoma;			// Persistent Na for soma
		Number					_gNaPProximal;		// Persistent Na for proximal dendrites
		Number					_gNaPDistal;		// Persistent Na for distal dendrites
		Number					_gIh;				// Ih conductance
		Number					_gKdr;				// Kdr conductance
		Number					_gKa;				// K-A transient cond
		Number					_gKc;				// Ca++ activated K+ cond
		Number					_gKm;				// K-M conductance
		Number					_gKahp;				// K-AHP conductance

		// Ca++ conductances under physiological conditions (e.g. 37-38 deg C)
		Number					_gCaTSoma;			// Ca-T for soma only
		Number					_gCaTDendrite;		// Ca-T for dendrites
		Number					_gCaN;				// Medium voltage act Ca++ chan
		Number					_gCaL;				// High voltage act Ca++ chan

		// Synaptic conductances. These are point conductance per synapse.
		Number					_gAMPA_PP;			// Peak AMPAR cond for PP
		Number					_gAMPA_AC;			// Peak AMPAR cond for AC
		Number					_gAMPA_MF;			// Peak AMPAR cond for MF
		Number					_gNR2A_PP;			// Peak NMDA-NR2A cond for PP
		Number					_gNR2A_AC;			// Peak NMDA-NR2A cond for AC
		Number					_gNR2A_MF;			// Peak NMDA-NR2A cond for MF
		Number					_gGABAa;			// Peak GABA-a conductance
		Number					_gGABAas;			// Peak GABA-a-slow conductance
		Number					_gGABAb;			// Peak GABA-b conductance

		// Special conductances for the axon
		Number					_gNaAxon;			// Axon Na+ channel density
		Number					_gKdrAxon;			// Axon Kdr channel density

		// Distance dependent ion channel conductance parameters ------

		Number					_maxDensityDist;	// max dist for setting cond
		Number					_maxBlendDist;		// max dist for blended chan
		Number					_maxCaDepKDist;		// max dist for Ca++ dep K+ chan

		Number					_gLeakSlope;		// dg/dx for gLeak=1/Rm
		Number					_gKdrSlope;			// dg/dx for Kdr
		Number					_gKaSlope;			// dg/dx for K-A
		Number					_gIhSlope;			// dg/dx for Ih

		// Provide a maximum radius of oblique (terminal) dendrites.
		// This allows adjustment of key conductances in small dendrites.
		// The power is applied as a multiplier on conductance as in:
		// g *= pow(obliqueRadius/radius,power).
		Number					_obliqueRadius;		// max oblique radius
		Number					_obliqueKaPower;	// power for scaling gKa

		// Parameters controlling ACh neuromodulation -----------------

		// Parameters for a mixed cation current modulated by ACh.
		Number					_gLeakAChModMC;		// mixed cation conductance
		Number					_VleakAChModMC;		// reversal potential
		Number					_VinitAChAdjust;	// change in Vinit from ACh

		// Values controlling ACh disable of the Na s-gate.
		// The s-gate is disabled when ACh>threshold and enabled otherwise
		Number					_NaSGateDisabledACh;	// threshold value
		Number					_NaSGateDisabledValue;	// value when disabled

		// Generally ACh neuromodulation is expressed as an adjustment to
		// peak conductance (even though many other effects are also likely).
		// These generally follow a Michaelis-Menten formula (well sort of):
		//
		// g_as_modulated = g_normal * (1+a*ACh/(ACh+Kd))
		//
		// Only a subset of channel types have ACh modulation defined.
		
		Number					_NaT_AChA;		// NaT params
		Number					_NaT_AChKd;		
		Number					_NaP_AChA;		// NaP params
		Number					_NaP_AChKd;		
		Number					_Ih_AChA;		// Ih params
		Number					_Ih_AChKd;
		Number					_Ka_AChA;		// K-A params
		Number					_Ka_AChKd;
		Number					_Km_AChA;		// K-M params
		Number					_Km_AChKd;
		Number					_Kahp_AChA;		// K-AHP params
		Number					_Kahp_AChKd;
		Number					_CaN_AChA;		// Ca-N params
		Number					_CaN_AChKd;
		Number					_CaL_AChA;		// Ca-L params
		Number					_CaL_AChKd;
		Number					_CaT_AChA;		// Ca-T params
		Number					_CaT_AChKd;

		// Structure parameters --------------------------------------

		// Provide orientation of cell with respect to laminar structure.
		// These coordinates define a unit normal vector for the
		// various strata of layers associated with the cell.
		Number					_orientationX;
		Number					_orientationY;
		Number					_orientationZ;

		// Functions for controlling neuromodulation ------------------

		// Set synaptic modulation based on ACh level
		virtual void			setGABAMod();
		virtual void			setGluMod();

		// Functions for building the cell ----------------------------

		// Build the different components of the cell
		// Subclass should explicitly invoke createCell to start the build.
		virtual void			createCell();			// creation driver
		virtual void			createSoma();			// build the soma compartment
		virtual void			createAxonAndIS();		// build axon and IS comps
		virtual void			createDendrites();		// build dendrite comps

		// Subclass responsibilities ----------------------------------			
		virtual void			setMorphology() = 0;	// Locate morphology table
		virtual void			setParameters() = 0;	// Initialize cell params
	};

	// -------------------------------------------
	// CA3 pyramidal cell neuron
	// -------------------------------------------

	class CA3PyramidalCell : public PyramidalCell {

	public:

		// Constructors and destructor
		CA3PyramidalCell(Model* m = NULL, bool doInit = true);
		virtual ~CA3PyramidalCell();

	protected:
		virtual void			setParameters();
	};

	// -------------------------------------------
	// L56a pyramidal cell neuron
	// -------------------------------------------

	class L56aPyramidalCell : public CA3PyramidalCell {

	public:

		// Constructors and destructor
		L56aPyramidalCell(Model* m = NULL, bool doInit = true);
		virtual ~L56aPyramidalCell();

	protected:
		virtual void			setMorphology();
	};

	// -----------------------------------------------
	// Function declarations for cell morphologies
	// which can potentially be used for neurons.
	// -----------------------------------------------

	// L56a CA3b cell from Duke/Southampton archive
	MorphologyEntry* cell_l56a_5_micron();
	MorphologyEntry* cell_l56a_25_micron();
	MorphologyEntry* cell_l56a_50_micron();

};

#endif // #ifndef
