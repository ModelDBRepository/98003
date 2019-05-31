// Synaptic Plasticity for Glutamate Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: plasticity_glu_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// synaptic plasticity for glumate synapses.
//
// AMPAR presynaptic plasticity is based on the model in
// Dittman et al (2000). Parameters for associational collaterals
// are based on studies of Schaffer collaterals. Parameters for
// mossy fiber synapses are based on PPF studies (Salin et al.)
// but have not been fully characterized experimentally.
//
// ACh modulation of glutamate release was measured by Hasselmo et al.
// but they did not address the effect on multiple pulses (PPF).
// Modulation included here is hypothetical based on the assumption
// that the net of effect of ACh is to reduce initial release
// probability and decrease presynaptic calcium influx in roughly
// the same measure. Results seem plausible but no experimental
// measurements of the effect of ACh on PPF have been done so far.
//
// Even more speculative is a phasic modulation of glutamate release
// within the theta cycle (Molyneaux & Hasselmo). Parameters are
// supplied to implement such a modulation, but experimental data
// does not address realistic spike trains in terms of theta phase
// effects (Ohliger-Frerking et al.). Possibly a more reasonable
// assumption is that interneuron firing frequency can be related 
// to a modulation of glutamate release. This is not implemented 
// here because interneuron firing rates are fixed. nAChR might
// also be a candidate mechanism among other possibilities.
//
// Synaptic conductance value releaseProbability is conceptually
// similar to the F1 parameter in the PPF model by Dittman et al.
// The F1 value is used here exclusively (see RandomReleaseRule 
// for a case where the synaptic conductance value is used instead).
//
// NMDAR dependent plasticity is based on experimental results
// from Debanne et al. (1998). This is specific to CA3 but uses
// cells derived from rather young animals aged in culture. 
// See Bi & Po for results using a different preparation. 
// See Sjostrom et al. for a similar type of plasticity in 
// L5 neocortical cells. More than a little extrapolation
// is needed to get to a composite model. 
//
// Froemke, Poo & Dan show Ca++ dependent NMDAR inactivation 
// following postsynaptic spikes leading to longer LTD times. 
// Ca++ dependent inactivation is included as a feature of the 
// model but parameters are difficult to set because the
// underlying mechanism is unclear (though see Grishin)
// 
// Note that in Grishin et al. (2004) CA1 and CA3 respond 
// differently in terms of ACh modulation of NMDAR currents. 
// When using sharp electrodes so as not to affect internal
// concentrations, there is a potentiation of Inmdar in both 
// CA1 and CA3 when muscarine is applied but the effect is 
// much smaller in CA3. Grishin et al. (2005) shows the
// pathway involved but getting the necessary [Ca++] and
// associated rate constants for simulation is another matter.
//
// A general pattern is that abstract classes implement the
// plasticity algorithm and subclasses supply parameter values.
// Some parameter values are stored as instance variables for 
// efficiency of access. These parameters can be changed during 
// the simulation, but this is not the objective.
// 
// References:
//
// Bi G-g and Poo M-m (1998). Synaptic Modifications in cultured
// hippocampal neurons: dependence onspike timing, synaptic strength,
// and postsynaptic cell type. J. Neurosci. 18, 10464-10472.
//
// Debanne D, Gahwiler BH, Thompson SM (1998). Long-term synaptic
// plsticity between pairs of individual CA3 pyramidal cells in rat
// hippocmapal slice cultures. J. Physiology 507, 237-247.
//
// Dittman JS, Kreitzer AC, and Regehr WG (2000). Interplay between
// facilitation, depression, and residual calcium at three presynaptic
// terminals. J. Neurosci. 20, 1374-1385.
//
// Flint AC, Maisch US, Weishaupt JH, Kriegstein AB, Moyner H (1997).
// NR2A subunit expression shortens NMDA receptor synaptic currents in
// developing neocortex. J. Neuroscience 17(7), 2469-2476.
//
// Froemke RC, Poo M-m, Dan Y (2005). Spike-timing-dependent synaptic
// plasticity depends on dendritic location. Nature 434: 221-225.
//
// Grishin AA, Gee CE, Gerber U, Benquet P (2004). Differential
// calcium-dependent modulation of NMDA currents in CA1 and CA3
// hipppocampal pyramidal cells. J. Neuroscience 24(2), 350-355.
//
// Hasselmo ME, Schnell E, Barkai E (1995). Dynamics of learning and recall
// at excitatory recurrent synapses and cholinergic modulation in rat
// hippocampal region CA3. J. Neurosci. 15, 5249-5262.
//
// Jahr, CW and Stevens, CF (1990). Voltage dependence of NMDA-activated
// macroscopic conductances predicted by single-channel kinetics.
// J. Neuroscience 10, 3178-3182.
//
// Kampa BM, Clements J, Jonas P, Stuart GJ (2004). Kinetics of Mg++ unblock
// of NMDA receptors: implications for spike-timing dependent synaptic
// plasticity. J. Physiology (London) 556, 337-345.
//
// Marino JM, Rouse ST, Levey AL, Potter LT, Conn PJ (1998). Activation
// of the genetically defined m1 muscarinic receptor potentiates
// N-methyl-D-aspartate (NMDA) receptor currents in hippocampal cells.
// PNAS 95, 11465-11470.
//
// Molyneaux BJ, Hasselmo ME (2002) GABA-b presynaptic inhibition has an
// in vivo time constant sufficiently rapid to allow modulation at
// theta frequency. J Neurophysiology 87, 1196-1205. 
//
// Ohliger-Frerking P, Wiebe SP, Staubli U, Frering M (2003). GABA-b
// receptor-mediated presynaptic inhibition has history-dependent
// effects on synaptic transmission during physiologically relevant
// spike trains. J Neurosci 23, 4809-4814.
//
// Salin PA, Scanziani M, Malenka RC, Nicoll RA (1996). Distinct short-term
// plasticity at two excitatory synapses in the hippocampus. 
// PNAS USA 93, 13304-13309.
//
// Sjostrom PJ, Turrigiano GG, Nelson SB (2001). Rate, timing, and
// cooperativity jointly determine cortical synaptic plasticity.
// Neuron 32, 1149-1164.


// Only include this header once
#ifndef __PLASTICITY_GLU_BAKER_2003_H_
#define __PLASTICITY_GLU_BAKER_2003_H_

#include "bnsf.h"
#include "synapse_glu_baker_2003.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation
namespace BAKER_2003 {

	// Class hierarchy
	class AMPARPresynapticState;
	class NMDARDepPlasticityState;

	class AMPA_PresynapticStateClassCache;
	class NMDARDepPlasticityRuleClassCache;

	class AMPARPresynapticRule;
		class CA3_AC_PairedPulseRule;
		class CA3_PP_PairedPulseRule;
		class CA3_MF_PairedPulseRule;

	class NMDARDepPlasticityRule;
		class CA3_STDPRule;



	// ----------------------------------------------------------------
	// Framework classes for Glu plasticity (mostly abstract)
	// ----------------------------------------------------------------

	// ================================================================
	// AMPAR Presynaptic Plasticity State
	// ================================================================

	class AMPA_PresynapticState {
	public:

		// State values from Dittman, Kreitzer, Wegher model.
		// See table 1 of the paper for the notation used there.
		Number					CaXF;		// Concentration of calcium-bound site XF
		Number					CaXD;		// Concentration of calcium-bound site XD
		Number					F;			// Facilitation variable
		Number					D;			// Depression variable

		// Flag to keep track of which synapses have received a spike
		// in the current time step. This is initialized to false.
		bool					spikeThisStep;	// true if spike just occurred
	};

	// ================================================================
	// AMPAR Presynaptic Plasticity Rule. By default, this implements
	// the model of Dittman et al. for Shaffer collaterals. Subclasses
	// may override to customize for other synapse types.
	// ================================================================

	class AMPA_PresynapticRule : public PlasticityRule {
	public:

		// Constructor and destructor
		AMPA_PresynapticRule();
		virtual ~AMPA_PresynapticRule();

		// Parameter value accessors 
		inline  Number			rho()				{ return _rho; }
		inline  Number			F1base()			{ return _F1base; }
		inline  SimTime			tauF()				{ return _tauF; }
		inline  SimTime			tauD()				{ return _tauD; }
		inline  Number			k0()				{ return _k0; }
		inline  Number			kmax()				{ return _kmax; }
		inline  Number			KD()				{ return _KD; }
		inline  Number			NT()				{ return _NT; }
		inline  Number			AChA()				{ return _AChA; }
		inline  Number			AChKd()				{ return _AChKd; }
		inline	Number			TPModA()			{ return _TPModA; }
		inline  Number			TPModFreq()			{ return _TPModFreq; }
		inline  Number			TPModPhase()		{ return _TPModPhase; }

		// Derived value accessors 
		inline  Number			KF()				{ return _KF; }
		inline  Number			F1()				{ return _F1; }
		inline  Number			deltaF()			{ return _deltaF; }
		inline  Number			deltaD()			{ return _deltaD; }
		inline  Number			AChMod()			{ return _AChMod; }

		// Accessors for setting values
		virtual void			rho(Number x);
		virtual void			F1base(Number x);
		virtual void			tauF(SimTime tau);
		virtual void			tauD(SimTime tau);
		virtual void			k0(Number x)			{ _k0=x; }
		virtual void			kmax(Number x)			{ _kmax=x; }
		virtual void			KD(Number x)			{ _KD=x; }
		virtual void			NT(Number x)			{ _NT=x; }
		virtual void			AChA(Number a);
		virtual void			AChKd(Number kd);
		virtual void			TPModA(Number a)		{ _TPModA=a; }
		virtual void			TPModFreq(Number hz)	{ _TPModFreq=hz; }	
		virtual void			TPModPhase(Number tp)	{ _TPModPhase=tp; }

		// Access the current ACh concentration level
		inline  Number			AChLevel()			{ return _AChLevel; }
		virtual void			AChLevel(Number ach);

		// Apply an ACh neuromodulation rule to adjust F1.
		// id = token("AChModulator"), nv=1, values[0] = ACh concentration
		virtual void			setModParams(TokenId id, int nv, Number* values);

		// Access a control value that specifies if release quantity
		// is to be either 0 or 1 based on a interpreting the
		// final model quantity as a probability. The alternative is
		// for a finalized event quantity to be a continuous value
		// as described in the original Dittman model.
		virtual bool			useRandomRelease() { return false; }

		// Access synapse release probability (for debug)
		virtual Number			releaseProbability(Synapse* syn);

		// Framework interfaces ---------------------------------------

		// Return the size to allocate for plasticity state data.
		virtual unsigned int	plasticityStateSize();

		// Create a plasticity state object in a synapse buffer.
		virtual void			createPlasticityState(Synapse* syn, Number wght);

		// Take action at the end of the time step (before AP purge)
		virtual void			applyEndOfStep(ActionPotentialEventQueue& apQueue);

		// Update state for a single AP event.
		virtual void			applyRule(ActionPotentialEvent* apEvent);

		// Set event quantity based on continuous quantity rule (default).
		virtual void			finalizeAPEvent(ActionPotentialEvent* apEvent);

		// Return a component name for reporting (should be overridden)
		virtual const char*		componentName() {return "AMPARPresynapticRule"; }

	protected:

		// Parameter values
		Number					_rho;		// Faciliatation ratio for closely spaced AP
		Number					_F1base;	// Initial release probability at [ACh]=0
		SimTime					_tauF;		// Decay time constant of CaXF after AP
		SimTime					_tauD;		// Decay time cosntant of CaXD after AP
		Number					_k0;		// Baseline recovery rate from refractory state
		Number					_kmax;		// Maximal recovery rate from refractory state
		Number					_KD;		// Affinity of CaXD for release site
		Number					_NT;		// Number of release sites (result scale factor)
		Number					_AChA;		// Maximum ACh effect (MM formulation)
		Number					_AChKd;		// ACh half activation concentration
		Number					_AChLevel;	// Current ACh concentration

		// The following parameters implement a modulation of release
		// probability that is modulated in a phasic manner within the
		// theta cycle (e.g. via GABA-b receptor activation ala Hasselmo).
		// Release probability multiplier is : 
		//
		// 1+a*(1+cos(2*pi*freq*t-phase))/2 where t is time.
		//
		// The multiplier is bounded below to be >=0. Note that phase is in radians.
		// By default, no theta phase modulation is provided, i.e. a=0.

		Number					_TPModA;	// Theta modulation amplitude (a)
		Number					_TPModFreq;	// Theta modulation frequency (freq)
		Number					_TPModPhase; // Theta modulation peak phase (phase)

		// Derived values
		Number					_KF;		// Affinity of CaXF for release site
		Number					_F1;		// Initial release probability as modulated
		Number					_deltaF;	// F-pathway calcium impulse as modulated
		Number					_deltaD;	// D-pathway calcium impulse as modulated
		Number					_AChMod;	// Derived ACh modulation effect

		// Cached values based on current time step size
		double					_ExpHTauF;	// exp(-h/tauF)
		double					_ExpHTauD;	// exp(-h/tauD)

		// Update derived values
		virtual void			updateAChMod();
		virtual void			updateKF();
		virtual void			updateCacheForStepSize();

		// Set the AP event quantity based on ppf or ppd value
		virtual void			updateQuantity(
			ActionPotentialEvent*	apEvent,
			Number					ppfd);
	};
	
	// ================================================================
	// NMDAR Dependent Postsynaptic Plasticity State Class
	// ================================================================

	class NMDARDepPlasticityState {
	public:

		// Time of previous presynaptic action potential
		Number					preSpikeTime;

		// State values representing fraction of hypothetical LTP molecules.
		// Such molecules are inspired by CamKII but do not attempt to model it.
		//
		// State model is:
		//
		//             >--(LTD Event)-->   
		//     +---> P                   D >--+
		//     |       <--(LTP Event)--<      |
		//     |                              | Timed decay
		//     |                              |
		//     +----(LTP Event)--< N <--------+
		//
		// States P and D lead to AMPAR expression. N=1-(P+D) is the fraction
		// that do not lead to AMPAR expression. Thus, transitions from N to P 
		// lead to LTP. Transitions from P to D do not immediately change 
		// anything, but as D decays to N, LTD occurs.

		Number					potentiationFraction;	// (P) potentiation state fraction
		Number					depressionFraction;		// (D) pre-depression state fractopm

		// Estimated fraction of NMDAR currently bound with glutamate.
		// The glycine receptor is assumed to be saturated in all cases.
		Number					gluBoundFraction;
	};

	// ================================================================
	// NMDAR Dependent Postsynaptic Plasticity Rule Class (Abstract)
	// ================================================================

	class NMDARDepPlasticityRule : public PlasticityRule {
	public:

		// Constructors and destructor
		NMDARDepPlasticityRule(NMDA_SynapticResp* nmdar);
		virtual ~NMDARDepPlasticityRule();

		// Accessor for the associated NMDAR response in this synapse
		inline  NMDA_SynapticResp* nmdar() { return _nmdar; }

		// Parameter value accessors ----------------------------------

		inline  Number			Wmax()					{return _Wmax; }
		inline	SimTime			tauW()					{return _tauW; }
		inline  SimTime			postSpikeMinISI()		{return _postSpikeMinISI; }
		inline  Number			ltpPostVm()				{return _ltpPostVm; }
		inline  Number			ltdPostVm()				{return _ltdPostVm; }
		inline  Number			Nnmdar()				{return _Nnmdar; }
		inline  Number			Popen()					{return _Popen; }
		inline  Number			PopenQ10()				{return _PopenQ10; }
		inline  Number			ratedTempC()			{return _ratedTempC; }
		inline  Number			Pbind()					{return _Pbind; }
		inline  SimTime			Topen()					{return _Topen; }
		inline  Number			nLTP()					{return _nLTP; }
		inline  Number			postSpikeNMDARSupp()	{return _postSpikeNMDARSupp; }
		inline  SimTime			tauCaDepSupp()			{return _tauCaDepSupp; }
		inline  SimTime			tauLTD()				{return _tauLTD; }
		inline	Number			PltdMin()				{return _PltdMin; }
		inline	Number			PltdMax()				{return _PltdMax; }
		inline  Number			PltdTau()				{return _PltdTau; }
		inline  Number			kLTP()					{return _kLTP; }
		inline  Number			kLTD()					{return _kLTD; }
		inline  Number			kFDP()					{return _kFDP; }

		// Accessors for setting parameter values ---------------------

		virtual void			Wmax(Number wm)					{_Wmax = wm; }
		virtual void			postSpikeMinISI(SimTime isi)	{_postSpikeMinISI = isi; }
		virtual void			ltpPostVm(Number vm)			{_ltpPostVm = vm; }
		virtual void			ltdPostVm(Number vm)			{_ltdPostVm = vm; }
		virtual void			Pbind(Number prob)				{_Pbind = prob; }
		virtual void			Topen(SimTime t)				{_Topen = t; }
		virtual void			nLTP(Number n)					{_nLTP = n; }
		virtual void			postSpikeNMDARSupp(Number a)	{_postSpikeNMDARSupp = a; }
		virtual void			PltdMin(Number prob)			{_PltdMin = prob; }
		virtual void			PltdMax(Number prob)			{_PltdMax = prob; }
		virtual void			PltdTau(Number tau)				{_PltdTau = tau; }
		virtual void			kLTP(Number k)					{_kLTP=k; }
		virtual void			kLTD(Number k)					{_kLTD=k; }
		virtual void			kFDP(Number k)					{_kFDP=k; }

		virtual void			Nnmdar(Number num);
		virtual void			Popen(Number pmax);
		virtual void			PopenQ10(Number q10);
		virtual void			ratedTempC(Number temp);
		virtual void			tauW(Number tau);
		virtual void			tauLTD(SimTime tau);
		virtual void			tauCaDepSupp(SimTime tau);

		// Accessors for param controlling selective knockout
		inline  bool			disableAChMod()			{ return _disableAChMod; }
		virtual void			disableAChMod(bool x);

		inline  bool			disableCaDepSupp()		{ return _disableCaDepSupp; }
		virtual void			disableCaDepSupp(bool x);

		// Accessor to get current target weight for a synapse (for debug)
		virtual Number			targetWeight(Synapse* syn);

		// Framework interfaces ---------------------------------------

		// Return the size to allocate for plasticity state data.
		virtual unsigned int	plasticityStateSize();

		// Create a plasticity state object in a synapse buffer.
		virtual void			createPlasticityState(Synapse* syn, Number wght);

		// Take action at the end of the time step (before AP purge)
		virtual void			applyEndOfStep(ActionPotentialEventQueue& apQueue);

		// Update state for a single AP event.
		virtual void			applyRule(ActionPotentialEvent* apEvent);

		// Return a component name for reporting (should be overridden)
		virtual const char*		componentName() {return "NMDARDepPlasticityRule"; }

		// Return an id for this type of plasticity
		virtual TokenId			plasticityTypeId() {return _NMDARDepPlasticityId; }

	protected:

		// Associated NMDAR to be referenced by this rule
		NMDA_SynapticResp*		_nmdar;

		// Parameter values (see related state variables for a discussion of states)
		// Subclass must set after base class construction.
		Number					_Wmax;		// Maximum weight value (when P+D=1)
		SimTime					_tauW;		// Time constant by which weights change

		SimTime					_postSpikeMinISI;	// Minimum ISI for LTP/LTD events
		Number					_ltpPostVm;			// Local Vm needed for LTP events
		Number					_ltdPostVm;			// Local Vm needed for LTD events

		Number					_Nnmdar;	// Average number of NMDAR per synapse
		Number					_Popen;		// Peak channel open prob after paired spikes
		Number					_PopenQ10;	// Q10 for adjusting Popen for temperature
		Number					_ratedTempC; // Temperature at which Popen is measured
		Number					_Pbind;		// Prob of binding an unbound NMDAR with Glu
		SimTime					_Topen;		// Time of first opening after Glu release.				
		Number					_nLTP;		// Threshold number of NMDAR open to induce LTP

		Number					_postSpikeNMDARSupp;	// NMDAR suppression after spike
		SimTime					_tauCaDepSupp;			// time const for Ca++ dependent supp

		SimTime					_tauUnbind;	// Time constant for unbinding Glu from NMDAR
		SimTime					_tauLTD;	// Time constant for moving from state D to N

		Number					_PltdMin;	// Min probability of LTD (eg when unpaired)
		Number					_PltdMax;	// Max probability of LTD (when delta t<0)
		Number					_PltdTau;	// Time constant for LTD probability

		Number					_kLTP;		// Fraction N->P for LTP inducing event
		Number					_kLTD;		// Fraction P->D for LTD inducing event
		Number					_kFDP;		// Fraction D->P for LTP inducing event

		// Parameters controlling selective knockout
		bool					_disableAChMod;		// Disable ACh effects on plasticity
		bool					_disableCaDepSupp;	// Disable Ca++ dep. suppression

		// Derived values
		Number					_Eopen;		// expected open NMDAR (Popen*Nnmdar)
		double					_ExpHTauW;			// exp(-h/tauW)
		double					_ExpHTauUnbind;		// exp(-h/tauUnbind)
		double					_ExpHTauLTD;		// exp(-h/tauLTD)
		double					_ExpHTauCaDepSupp;	// exp(-h/tauCaDepSupp)

		// Local state (other than synapse state)
		SimTime					_ltpPostEventTime;	// Time of previous LTP threshold crossing
		SimTime					_ltdPostEventTime;	// Time of previous LTD threshold crossing
		Number					_prevVm;			// Vm value for previous time step
		Number					_prevVmDot;			// VmDot at end of previous time step
		Number					_NMDARCaDepSupp;	// NMDAR Ca dep activity suppression

		// Provide a static token id for this type of plasticity
		static const TokenId	_NMDARDepPlasticityId;

		// Check for changes in time step size and adjust cached values
		virtual void			updateCacheForStepSize();

		// Update precomputed values
		virtual void			updateEopenValue();

		// Return the estimated peak Vm of the postsynaptic compartment.
		// Current Vm is saved for the next step as a side-effect so this
		// should only be invoked once per time step.
		virtual Number			peakPostsynapticVm();

		// Update presynaptic plasticity state based on pre-before-post
		// spike pairings and the passage of time.
		virtual void			updatePresynapticState(Synapse* syn, bool spikeNow);

		// Get probability of LTP inducing event given expected number of NMDAR open
		virtual Number			probLTP(Number nopen);

		// Get probability of LTD inducing event given postSpikeTime-PreSpikeTime
		virtual Number			probLTD(SimTime dt);
	};


	// ----------------------------------------------------------------
	// Concrete instance classes for Glu plasticity
	// ----------------------------------------------------------------

	// ================================================================
	// CA3 Mossy Fiber Paired-pulse Rule
	// ================================================================

	class CA3_MF_PairedPulseRule : public AMPA_PresynapticRule {
	public:

		// Constructor and destructor
		CA3_MF_PairedPulseRule()
		{
			using namespace UOM;

			// Update parameter values. These values 
			// are found by optimal fitting to PPF
			// values found in Salin et al. for MF
			// synapses onto CA3 pyramidal cells. The 
			// longer term presynaptic plasticity found 
			// in Salin et al. is not provided here.

			// This is a low probability release synapse. 
			// The fitting process does not supply an F1 
			// value (optimization of the fit decreases F1 
			// to zero which is impossible). Assuming
			// an NT value of 30 gives the F1 value used.
			// The optimization process is used to select 
			// rho and tauF only.

			F1base(	 0.033f );
			rho( 3.14f );
			tauF( 170*msec);
			NT( 30 );
		}

		virtual ~CA3_MF_PairedPulseRule() {}

		// Return a component name for reporting
		virtual const char*	componentName() {return "CA3_MF_PPRule"; }
	};

	// ================================================================
	// CA3 Associational Collateral Paired-pulse Rule
	// ================================================================

	class CA3_AC_PairedPulseRule : public AMPA_PresynapticRule {
	public:

		// Constructor and destructor
		CA3_AC_PairedPulseRule()
		{
			// Update parameter values. The value of rho
			// is changed from that in the SC model to better
			// match PPF values found in Salin et al.
			rho( 2.0f );

			// Provide ACh parameters based on a fit to Hasselmo et al.
			// This modulation is applied to AC synapses only.
			// SC synapses are similar but have different sensitivities.
			AChA( -0.8f );
			AChKd( 45*UOM::microM );
		}

		virtual ~CA3_AC_PairedPulseRule() {}

		// Use binary random release to set final quantity
		virtual bool			useRandomRelease() { return true; }

		// Return a component name for reporting
		virtual const char*		componentName() {return "CA3_AC_PPRule"; }
	};

	
	// ================================================================
	// CA3 spike-time dependent plasticity rule for assoc. collaterals
	// ================================================================

	class CA3_AC_STDPRule : public NMDARDepPlasticityRule {
	public:

		// Constructors and destructor
		CA3_AC_STDPRule(NMDA_SynapticResp* nmdar) : NMDARDepPlasticityRule(nmdar)
		{
			using namespace UOM;

			// Set initial parameter values.
			// These are chosen to match results for CA3 synapses,
			// primarily results from Debanne et al assuming 60 (or 100)
			// pairings for the spike pair STDP curve shown there.

			// Parameters are adjusted for current temperature using an
			// NMDAR Q10 value of 3. Nopen Q10 is set to match STDP data
			// in Debanne et al. by temperature adjusting Bi & Poo results.

			Wmax			( 3.2f );		// Max weight (Debanne max LTP)
			tauW			( 10*sec );		// See Shouval et al BC 87:383-391

			postSpikeMinISI	( 3*msec );		// Minimum time between spikes
			ltpPostVm		( -30*mV );		// Postsynaptic threshold for LTP
			ltdPostVm		( -40*mV );		// Postsynaptic threshold for LTD
			
			ratedTempC		( 25 );			// Temp at which Popen was measured
			Popen			( 0.02f );		// From Kampa et al. for brief AP
			PopenQ10		( 1.8f );		// Q10 for nopen probLTP computation

			Nnmdar			( 400 );		// Cottrell et al.
			nLTP			( 5 );			// Fit based on Bi & Poo + Debanne at +15ms

			Pbind			( 0.6f );		// Mainen et al. estimate
			Topen			( 10*msec );	// Interval over which LTD dominates

			tauLTD			( 70*msec );	// See Wang et al. 2005 + Sjostrom et al.
			PltdMin			( 0.06f );		// Debanne homosynaptic LTD fit
			PltdMax			( 0.45f );		// Debanne LTD fit extrapolated to dt=0;
			PltdTau			( 90*msec );	// Debanne LTD fit at -15ms and -200ms

			kLTP			( 0.012f );		// Fit to Debanne LTP at dt=+15ms (60 pairings)
			kLTD			( 0.012f );		// Fit to Debanne LTD at dt=-15ms (100 pairings)
			kFDP			( 1.0f );		// LTP cancels immediately preceeding LTD

			tauCaDepSupp	( 90*msec );	// See Froemke Nature 434:221-225 fig 2
			postSpikeNMDARSupp ( 0.3f );	// See Froemke. See also Grishin et al.
		}

		virtual ~CA3_AC_STDPRule() {}

		// Return a component name as needed for reporting
		virtual const char*	componentName() {return "CA3_AC_STDP"; }
	};

	// ================================================================
	// CA3 spike-time dependent plasticity rule for perforant path
	// ================================================================

	class CA3_PP_STDPRule : public NMDARDepPlasticityRule {
	public:

		// Constructors and destructor
		CA3_PP_STDPRule(NMDA_SynapticResp* nmdar) : NMDARDepPlasticityRule(nmdar)
		{
			using namespace UOM;

			// Set initial parameter values. See comments for AC values.
			// The PP rule exists separately to permit experimentation
			// with different parameter settings.

			Wmax			( 3.2f );		// Max weight (Debanne max LTP)
			tauW			( 10*sec );		// See Shouval et al BC 87:383-391

			postSpikeMinISI	( 3*msec );		// Minimum time between spikes
			ltpPostVm		( -30*mV );		// Postsynaptic threshold for LTP
			ltdPostVm		( -40*mV );		// Postsynaptic threshold for LTD
			
			ratedTempC		( 25 );			// Temp at which Popen was measured
			Popen			( 0.02f );		// From Kampa et al. for brief AP
			PopenQ10		( 1.8f );		// Q10 for nopen probLTP computation

			Nnmdar			( 400 );		// Cottrell et al.
			nLTP			( 5 );			// Fit based on Bi & Poo + Debanne at +15ms

			Pbind			( 0.6f );		// Mainen et al. estimate
			Topen			( 10*msec );	// Interval over which LTD dominates

			tauLTD			( 70*msec );	// See Wang et al. 2005 + Sjostrom et al.
			PltdMin			( 0.06f );		// Debanne homosynaptic LTD fit
			PltdMax			( 0.45f );		// Debanne LTD fit extrapolated to dt=0;
			PltdTau			( 90*msec );	// Debanne LTD fit at -15ms and -200ms

			kLTP			( 0.012f );		// Fit to Debanne LTP at dt=+15ms (60 pairings)
			kLTD			( 0.012f );		// Fit to Debanne LTD at dt=-15ms (100 pairings)
			kFDP			( 1.0f );		// LTP cancels immediately preceeding LTD

			tauCaDepSupp	( 90*msec );	// See Froemke Nature 434:221-225 fig 2
			postSpikeNMDARSupp ( 0.3f );	// See Froemke. See also Grishin et al.
		}


		virtual ~CA3_PP_STDPRule() {}

		// Return a component name as needed for reporting
		virtual const char*	componentName() {return "CA3_PP_STDP"; }
	};
};


#endif // #ifndef

