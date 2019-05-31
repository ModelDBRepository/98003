// Synaptic Plasticity for Glutamate Synapses
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: plasticity_glu_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// glutamate synapse synaptic plasticity.
//
// See header file for references.

#include "plasticity_glu_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace BAKER_2003;


// =================================================
// AMPA_PresynapticRule Class Body
// =================================================


// Constructors and destructor
AMPA_PresynapticRule::AMPA_PresynapticRule()
{
	using namespace UOM;

	// Set default values from Dittman et al. 
	// for Schaffer collateral synapses.
	// Subclass may reset as needed.
	_rho		= 2.2f;
	_F1base		= 0.24f;
	_tauF		= 100*msec;
	_tauD		= 50*msec;
	_k0			= 2/sec;
	_kmax		= 30/sec;
	_KD			= 2;
	_NT			= 1;

	// Set ACh modulation to have no effect
	_AChA		= 0;		// no effect
	_AChKd		= 1*microM;	// arbitrary non-zero value
	_AChLevel	= 0;		// assume no ACh at start

	// Set theta modulation to have no effect
	_TPModA		= 0;		// amplitude
	_TPModFreq	= 8*Hz;		// theta frequency
	_TPModPhase = 0;		// default only

	// Set values based on ACh level even
	// though this is basically trivial here.
	// Subclasses would use normal accessors
	// that include automatic updates of AChMod.
	updateAChMod();

	// Initialize KF to force evaluation on first use.
	// This allows F1base and rho to be set in subclasses
	// without triggering an error when the derived KF<0.
	_KF			= 0;
}

AMPA_PresynapticRule::~AMPA_PresynapticRule() {}

// Set rho value
void AMPA_PresynapticRule::rho(Number x)
{
	_rho = x;
	_KF  = 0;
}

// Set F1 base value
void AMPA_PresynapticRule::F1base(Number x)
{
	// Save value and update any derived values
	_F1base = x;
	_KF  = 0;			// Defer update of KF until later
	updateAChMod();		// Get new F1 value
}

// Set tauF value
void AMPA_PresynapticRule::tauF(SimTime tau)
{
	_tauF = tau;
	updateCacheForStepSize();
}

// Set tauD value
void AMPA_PresynapticRule::tauD(SimTime tau)
{
	_tauD = tau;
	updateCacheForStepSize();
}

// Set ACh maximum effect
void AMPA_PresynapticRule::AChA(Number a)
{
	_AChA = a;
	updateAChMod();
}

// Set ACh half activation value
void AMPA_PresynapticRule::AChKd(Number kd)
{
	_AChKd = kd;
	updateAChMod();
}

// Set ACh concentration
void AMPA_PresynapticRule::AChLevel(Number ach)
{
	// Save the value and update derived values
	_AChLevel = ach;
	updateAChMod();
}

// Compute the net modulation effect of ACh
void AMPA_PresynapticRule::updateAChMod()
{
	// Compute the derived effect of ach
	_AChMod = 1 + AChA()*AChLevel()/( AChLevel()+AChKd() );

	// Update derived values affected by ACh
	_F1 = F1base() * _AChMod;
	_deltaF = 1 * _AChMod;		// unit impulse at ACh=0
	_deltaD = 1 * _AChMod;		// unit impulse at ACh=0
}

// Update derived values based on current settings
// This is done when KF is first used to avoid the
// problem of invalid values when f1 and rho must
// be changed together.
void AMPA_PresynapticRule::updateKF()
{
	// Get values for f1 and rho prior to neuromodulation.
	// The assumption here is that neuromodulation does not
	// affect the KF value derived here.
	Number f1 = F1base();

	// See Dittman et al. KF can be inferred from rho and F1.
	_KF = (1-f1)/(rho()*f1/(1-f1)-f1) - 1;

	// Make sure that a valid value results
	if ( !(_KF>0) ) {
		FatalError("(AMPA_PresynapticRule::updateDerivedValues) "
			"rho and F1 together imply a negative KF value");
	}
}

// Update cached values based on time step
void AMPA_PresynapticRule::updateCacheForStepSize()
{
	// Make sure there is a model from which to
	// get the current step size. If not stop now.
	if (model()==NULL)
		return;

	SimTime h = timeStepSize(); // Size of previous step or else 0

	_ExpHTauF = exp( -h/tauF() );
	_ExpHTauD = exp( -h/tauD() );
}

// Return the size to allocate for plasticity state data.
unsigned int AMPA_PresynapticRule::plasticityStateSize()
{
	return sizeInBytesRounded(sizeof(AMPA_PresynapticState));
}

// Create a plasticity state object in a synapse buffer.
void AMPA_PresynapticRule::createPlasticityState(Synapse* syn, Number wght)
{
	AMPA_PresynapticState* psd = 
		(AMPA_PresynapticState*) plasticityData(syn);

	psd->CaXF	= 0;
	psd->CaXD	= 0;
	psd->F		= F1();
	psd->D		= 1;

	psd->spikeThisStep = false;
}

// Take action at the end of the time step (before AP purge)
void AMPA_PresynapticRule::applyEndOfStep(ActionPotentialEventQueue& apQueue)
{
	SimTime h = timeStepSize(); // Size of previous step or else 0

	// Check for out-of-date derived values
	if (_KF==0) {
		updateKF();
		updateAChMod();
	}
	if (model()->stepSizeChanged() ) {
		updateCacheForStepSize();
	}

	// Let superclass loop through current AP events.
	// This just tags synapses for which there was an AP
	// with an event within the time interval of this step.
	PlasticityRule::applyEndOfStep(apQueue);

	// Go through synapses one at a time and update the state
	Synapse* syn = synapticCond()->synapses();
	while(syn!=NULL) {

		AMPA_PresynapticState* psd = 
			(AMPA_PresynapticState*) plasticityData(syn);

		// Solve the following ODEs explicitly: 
		//
		// dCaXF/dt = -CaXF/tauF + deltaF*ispk(t)
		// dCaXD/dt = -CaXD/tauD + deltaD*ispk(t)
		// dD/dt = (1-D)*krecov - D*F*ispk(t)
		// 
		// where: ispk(t) = unit impulse if spike at time t,
		// F(t) = F1+(1-F1)*CaXF(t)/(CaXF(t)+KF), and
		// krecov(t) = k0+(kmax-k0)*CaXD(t)/(CaXD(t)+KD)
		//
		// Note that we have to solve these equations step-by-step
		// because of non-linearities. When a  spike arrives, F 
		// changes discontinuously. For purposes  of solving the 
		// ODEs  we treat the spike as arriving at the end of the 
		// time step so that the discontinuity is not reflected in 
		// the new value of D at that time. The solution here is
		// still an approximation based on the assumption that
		// time steps are much shorter than relevant time constants
		// but is more than close enough for the present purposes.

		Number krecov = k0() + (kmax()-k0())*psd->CaXD/(psd->CaXD+KD());

		psd->CaXF = psd->CaXF * _ExpHTauF;
		psd->CaXD = psd->CaXD * _ExpHTauD;

		psd->F = F1() + (1-F1())*psd->CaXF/(psd->CaXF+KF());
		psd->D = 1-(1-psd->D)*qdexp(-h*krecov);

		if (psd->spikeThisStep) {

			// Adjust for the unit impulse is(t).
			psd->D -= psd->D*psd->F;
			psd->CaXF += deltaF();
			psd->CaXD += deltaD();

			// Clear the spike flag for next time
			psd->spikeThisStep = false;
		}

		// Move on to the next synapse
		syn = syn->nextPostsynaptic();
	}
}

// Update state for a single AP event.
void AMPA_PresynapticRule::applyRule(ActionPotentialEvent* apEvent)
{
	// Tag the synapse as having received an AP this time step
	AMPA_PresynapticState* psd = 
		(AMPA_PresynapticState*) plasticityData(apEvent->synapse());

	psd->spikeThisStep = true;
}

// Access synapse release probability for debug
Number AMPA_PresynapticRule::releaseProbability(Synapse* syn)
{
	AMPA_PresynapticState* psd = 
		(AMPA_PresynapticState*) plasticityData(syn);

	return NT() * psd->F * psd->D;
}

// Set event quantity based on rule
void AMPA_PresynapticRule::finalizeAPEvent(ActionPotentialEvent* apEvent)
{
	// Set quantity to average response and tag event as final.

	// This uses the plasticity state as of the start of the time
	// step in which the AP event was processed. There is a 
	// possibility that an adaptive solver may decide to use a
	// smaller time step, in which case the AP event can be
	// finalized with a quantity that would be changed slightly
	// in a later time step. This should be a very minor effect 
	// given the size of tauF and tauD in terms of time steps.

	// Note that NT is just a scale factor for this process.
	// In many cases, NT will be approximately 1/F1 to give
	// a unit response to the first spike of a series.

	AMPA_PresynapticState* psd = 
		(AMPA_PresynapticState*) plasticityData(apEvent->synapse());

	Number ppfd = NT() * psd->F * psd->D;

	// Apply theta phase modulation, if any.
	if (TPModA()!=0) {

		SimTime		t = apEvent->eventTime();
		Number		a = TPModA();
		Number		f = TPModFreq();
		Number		p = TPModPhase();

		ppfd *= plusval(1+a*(1+cos(2*Pi*f*t-p))/2);
	}	

	updateQuantity(apEvent, ppfd);
	apEvent->isFinal( true );
}

// Update the AP event based on a ppf/ppd value
void AMPA_PresynapticRule::updateQuantity(
	ActionPotentialEvent* apEvent, 
	Number ppfd)
{
	if ( useRandomRelease() ) {

		// Use the current postsynaptic release probability to set
		// quantity to 1 (release) or 0 (failure). A Bernoulli-trial
		// release is all that is supported here.
		apEvent->quantity( synapticCond()->runif()<ppfd ? 1 : 0 );
	}
	else {

		// Otherwise treat release as continuous 
		// This is as in the original model by Dittman et al.
		apEvent->quantity( ppfd );
	}
}

// Apply an ACh neuromodulation rule
void AMPA_PresynapticRule::setModParams(TokenId id, int nv, Number* values)
{
	using namespace UOM;

	static const TokenId AChMod = token("AChModulator");

	// Skip any other forms of modulation and check num of params
	if (id!=AChMod) return;
	if (nv<1) {
		FatalError("(AMPA_PresynapticRule::setModParams::setModParams) "
			"Too few params");
	}

	// Set the new concentration value
	AChLevel( values[0] );
}



// =================================================
// NMDARDepPlasticityRule Class Body
// =================================================



// Static values
const TokenId NMDARDepPlasticityRule::_NMDARDepPlasticityId 
	= token("NMDARDepPlasticityRule");

// Constructors and destructor
NMDARDepPlasticityRule::NMDARDepPlasticityRule(NMDA_SynapticResp* nmdar)
{
	using namespace UOM;

	// Save associated NMDAR reference
	_nmdar = nmdar;

	// Set initial values for local state
	_ltpPostEventTime = InfinitePast;
	_ltdPostEventTime = InfinitePast;
	_prevVm = VMinForIndex;
	_prevVmDot = 0;
	_NMDARCaDepSupp = 0;
	_Eopen = 0;

	// Fill in parameters with initial values.
	// Values here show typical units of measure.
	// These are only placeholders and should be
	// updated during subclass construction.
	_ratedTempC			= 25;
	_Wmax				= 1;
	_tauW				= 1*sec;
	_postSpikeMinISI	= 1*msec;
	_ltpPostVm			= -30*mV;
	_ltdPostVm			= -40*mV;
	_Nnmdar				= 0;
	_Popen				= 0;
	_PopenQ10			= 1;
	_Pbind				= 0;
	_Topen				= 1*msec;
	_nLTP				= 0;
	_postSpikeNMDARSupp	= 0;
	_tauLTD				= 100*msec;
	_tauCaDepSupp		= 100*msec;
	_PltdMin			= 0;
	_PltdMax			= 0;
	_PltdTau			= 100*msec;
	_kLTP				= 0;
	_kLTD				= 0;
	_kFDP				= 0;

	// Set selective knockouts to have no effect
	_disableAChMod		= false;
	_disableCaDepSupp	= false;
}

NMDARDepPlasticityRule::~NMDARDepPlasticityRule() {}

// Set the temperature corresponding to Popen
void NMDARDepPlasticityRule::ratedTempC(Number temp)
{
	_ratedTempC = temp;
	updateEopenValue();
}

// Set number of NMDAR and update cached value
void NMDARDepPlasticityRule::Nnmdar(Number num)
{
	_Nnmdar = num;
	updateEopenValue();
}


// Set peak open probability of an NMDAR and reset cache
void NMDARDepPlasticityRule::Popen(Number pmax)
{
	_Popen = pmax;
	updateEopenValue();
}

// Set peak open probability of an NMDAR and reset cache
void NMDARDepPlasticityRule::PopenQ10(Number q10)
{
	_PopenQ10 = q10;
	updateEopenValue();
}

// Set tau value and update derived values
void NMDARDepPlasticityRule::tauW(Number tau)
{
	_tauW = tau;
	updateCacheForStepSize();
}

// Set tau value and update derived values
void NMDARDepPlasticityRule::tauCaDepSupp(SimTime tau)
{
	_tauCaDepSupp = tau;
	updateCacheForStepSize();
}

// Set tau value and update derived values
void NMDARDepPlasticityRule::tauLTD(SimTime tau)
{
	_tauLTD = tau;
	updateCacheForStepSize();
}

// Set the ACh Modulation status
void NMDARDepPlasticityRule::disableAChMod(bool x)
{
	_disableAChMod = x;
	updateEopenValue();
}

// Set the Ca++ dependent suppression status
void NMDARDepPlasticityRule::disableCaDepSupp(bool x)
{
	_disableCaDepSupp = x;
}

// Return target weight for the synapse provided
Number NMDARDepPlasticityRule::targetWeight(Synapse* syn)
{
	// Locate the postsynaptic state data
	NMDARDepPlasticityState* psd = 
		(NMDARDepPlasticityState*) plasticityData(syn);

	// Compute the current target weight and return it
	return Wmax() * (psd->potentiationFraction + psd->depressionFraction);
}

// Update derived Eopen value
void NMDARDepPlasticityRule::updateEopenValue()
{
	_Eopen = _Nnmdar * _Popen;
	_Eopen *= pow(PopenQ10(), (nmdar()->currentTempC()-ratedTempC())/10);
}

// Adjust cached values to reflect current time step size
void NMDARDepPlasticityRule::updateCacheForStepSize()
{
	// Make sure there is a model from which to get step size.
	// Otherwise, stop now.
	if (model()==NULL)
		return;

	// In this, we assume that nmdar tau2 does not change
	// after the initial start up setting done here.

	SimTime h = timeStepSize();
	SimTime tauUnbind = nmdar()->tau2()/nmdar()->Q10Factor();

	_ExpHTauW			= exp(-h/tauW() );
	_ExpHTauUnbind		= exp(-h/tauUnbind );
	_ExpHTauLTD			= exp(-h/tauLTD() );
	_ExpHTauCaDepSupp	= exp(-h/tauCaDepSupp() );
}

// Return the size to allocate for plasticity state data.
unsigned int NMDARDepPlasticityRule::plasticityStateSize()
{
	return sizeInBytesRounded(sizeof(NMDARDepPlasticityState));
}

// Create a plasticity state object in a synapse buffer.
void NMDARDepPlasticityRule::createPlasticityState(Synapse* syn, Number wght)
{
	NMDARDepPlasticityState* psd = 
		(NMDARDepPlasticityState*) plasticityData(syn);

	psd->preSpikeTime = InfinitePast;
	psd->potentiationFraction = wght/Wmax();
	psd->depressionFraction = 0;
	psd->gluBoundFraction = 0;
}

// Take action at the end of the time step (before AP purge)
void NMDARDepPlasticityRule::applyEndOfStep(ActionPotentialEventQueue& apQueue)
{
	SimTime		Tnow = currentTime();
	bool		postSpikeNow = false;

	// Check for change in time step size and adjust derived values.
	if (model()->stepSizeChanged() ) {
		updateCacheForStepSize();
	}
	
	// Use superclass logic to invoke applyRule for each AP event.
	// This handles the post-before-pre pairings leading to LTD.
	// Note that this uses postsynaptic spike times as of the
	// start of the time step interval. Postsynaptic spikes
	// during the time step are handled as a special case of
	// updatePresynapticState.
	PlasticityRule::applyEndOfStep(apQueue);

	// Advance shared postsynaptic state to end of step.
	// See note at the end of this logic wrt time lag
	// in calcium dependent suppression.
	_NMDARCaDepSupp *= _ExpHTauCaDepSupp;

	// Get the estimated peak voltage during the time step.
	// Note that this function must only be invoked once
	// per time step since it saves the previous voltage
	// as a side-effect of the making the estimate.
	Number		peakVm = peakPostsynapticVm();

	// See if the postsynaptic potential indicates that
	// a postsynaptic spike occurred by exceeding the LTP threshold.
	// MinISI rules for post synaptic spike detection filters out
	// small oscillations that would otherwise be multiple spikes.
	// LTD threshold is not necessarily a spike though exceeding the
	// LTD threshold triggers LTD processing of presynaptic spikes.
	if (peakVm >= ltpPostVm() && 
		Tnow >= _ltpPostEventTime+postSpikeMinISI() ) {
		_ltpPostEventTime = Tnow;
		postSpikeNow = true;
	}
	if (peakVm>=ltdPostVm() ) {
		_ltdPostEventTime = Tnow;
	}

	// Go through all synapses one at a time to
	// update state to the end of the time step.
	// Any postsynaptic spikes are treated as if they
	// occurred at the end of the time step.
	Synapse*	syn = synapticCond()->synapses();
	while (syn!=NULL) {

		// Update presynaptic plasticity state
		updatePresynapticState(syn,postSpikeNow);

		// Move on to the next synapse
		syn = syn->nextPostsynaptic();
	}

	// If a postsynaptic spike has just occurred, set a new value
	// for Ca++ dependent NMDAR suppression for use in the next step.
	// This is done after the previous loop under the assumption that
	// there is a non-zero time lag before NMDAR suppression, or,
	// put another way, Ca++ dependent suppression is included as an
	// integral part of the probability of LTP events.
	if (postSpikeNow) {

		// Do not set depression amount if currently disabled.
		// Otherwise set to value immediately after the spike.
		if (disableCaDepSupp() ) {
			_NMDARCaDepSupp = 0;
		}
		else {
			// Set depression amount such that even with the
			// current ACh modulation applied multipliciatively,
			// the result is initially the depression specified
			// by postSpikeNMDARSupp(). This is based on the
			// observation that post spike suppression, Ca++
			// inactivation and ACh modulation with raised [Ca++]i
			// yield about the same suppression of Inmdar.
			_NMDARCaDepSupp = 1-(1-postSpikeNMDARSupp())/nmdar()->AChMod();
		}
	}
}

// Update presynaptic plasticity state based on postsynaptic activity
void NMDARDepPlasticityRule::updatePresynapticState(
	Synapse* syn,			// Synapse with state to be updated
	bool spikeNow )			// Flag indicating an LTP inducing spike
{
	// Locate the postsynaptic state data and synapse weight
	NMDARDepPlasticityState* psd = 
		(NMDARDepPlasticityState*) plasticityData(syn);

	SimTime		Tnow = currentTime();
	Number		weight = synapticWeight(syn);
	Number		targetWeight,newWeight;

	// Update synapse states to reflect the passage of time
	psd->gluBoundFraction *= _ExpHTauUnbind;
	psd->depressionFraction *= _ExpHTauLTD;

	// Adjust state based on any pre-before-post pairing (mostly LTP).
	// LTP events both induce new LTP (N->P) and cancel pending LTD (D->P).
	if (spikeNow) {
		if (Tnow >= psd->preSpikeTime) {

			// Handle LTP Event

			Number		nopen;			// expected number of open receptors
			Number		Nfraction;		// N(eutral) state fraction where N+P+D=1
			Number		pltp;			// probability of LTP event
			Number		cltd;			// cancelled LTD fraction

			// Based on the number of NMDAR open, determine the probability
			// of LTP induction and adjust state values accordingly.

			// For the interval between the presynaptic spike arrival up to
			// a time Topen later, NMDAR are assumed to open at a linear
			// rate so that initially the number open is zero and after Topen
			// it is whatever the decay time constant indicates. This precludes
			// LTP induction in the very early interval of a simultaneous
			// pre- and postsynatic spike pairing.

			nopen = _Eopen * psd->gluBoundFraction *  (1-_NMDARCaDepSupp);
			if (Tnow - psd->preSpikeTime < Topen() ) {
				nopen *= (Tnow - psd->preSpikeTime)/Topen();
			}
			if (!disableAChMod() ) {
				nopen *= nmdar()->AChMod();
			}
			pltp = probLTP(nopen);

			cltd = kFDP() * pltp * psd->depressionFraction;
			Nfraction = 1-(psd->potentiationFraction+psd->depressionFraction);
			psd->potentiationFraction += cltd + kLTP()*pltp*Nfraction;	
			psd->depressionFraction -= cltd;

			// The following logic attempts to handle collisions between pre- and 
			// postsynaptic activity which should result in net LTD.

			if (Tnow - psd->preSpikeTime < Topen()) {

				// If a backpropagating spike is detected before NMDAR would be open,
				// then LTD is assumed in an amount proportional to the fraction not
				// open fraction found for the NMDAR. This rule is primarily important
				// when pre- and postsynaptic activity occur exactly together, for 
				// example, because of local dendritic spikes. 

				Number pltd = probLTD(Tnow - psd->preSpikeTime);
				Number dltd = kLTD()*pltd*psd->potentiationFraction;
				psd->potentiationFraction -= dltd;
				psd->depressionFraction += dltd;
			}
		}
	}

	// Get the new target weight and update the synapse weight in that direction.
	// This is a pretty simple version of what happens in the synapse following
	// LTP induction, but at least it decouples LTP events and an immediate
	// change in synapse weights. Minimum synaptic weight is set at 0 to avoid
	// small roundoff errors from forcing the weight slightly negative.
	targetWeight = Wmax()*(psd->potentiationFraction + psd->depressionFraction);
	newWeight = targetWeight+(weight-targetWeight) * _ExpHTauW;
	synapticWeight(syn,plusval(newWeight));
}

// Update state for a single AP event.
void NMDARDepPlasticityRule::applyRule(ActionPotentialEvent* apEvent)
{
	// Quantile amount of 0 indicates that no release occurred.
	// In that case, skip the rest of event processing.
	if (apEvent->quantity() == 0)
		return;

	Synapse* syn = apEvent->synapse();
	NMDARDepPlasticityState* psd = 
		(NMDARDepPlasticityState*) plasticityData(syn);

	SimTime					dt;
	Number					dltd;
	Number					pb;

	// Get the new number of receptors bound based on binding probability.
	// Estimate probability of binding as proportional to AP quantity
	// (usually quantity=1) with an upper bound of pb=1.
	pb = minval(1, apEvent->quantity()*Pbind() );
	psd->gluBoundFraction += pb*(1-psd->gluBoundFraction);

	// Adjust state for post-before-pre case (LTD).
	dt = _ltdPostEventTime - apEvent->eventTime();
	dltd = kLTD() * probLTD(dt) * psd->potentiationFraction;
	psd->potentiationFraction -= dltd;
	psd->depressionFraction += dltd;

	// Save the time of the presynaptic spike
	psd->preSpikeTime = apEvent->eventTime();
}

// Return an estimated peak Vm of the postsynaptic compartment.
Number NMDARDepPlasticityRule::peakPostsynapticVm()
{
	// Compute an estimate based on the Vm as of the end of the
	// current step, derivative, and the previous Vm and derivative
	// as of the end of the previous time step.

	Number		vm		= synapticCond()->Vm();
	Number		vmDot	= synapticCond()->VmDot();
	SimTime		h		= model()->timeStepSize();
	Number		peakVm;

	// Start with the maximum of current and previous Vm
	peakVm = maxval(vm,_prevVm);

	// If there is a maximum voltage within the time step
	// make an estimate of it using derivatives
	if (_prevVmDot>0 && vmDot<=0) { 
		
		// Estimate the peak by finding the point where
		// linear approximations to vm(t) match at beginning
		// and ending of the time step, i.e. solve
		//
		// v(t0+dt) = prevVm + dt*prevVmDot = vm + (dt-h)*vmDot
		//
		// where t0 is the step start time,
		// dt is a time offset into the step, and
		// h is the time step size (thus 0<=dt<=h).
		//
		// If there is no solution, then leave the
		// earlier estimate in place.

		SimTime	dt = (vm-h*vmDot-_prevVm)/(_prevVmDot-vmDot);

		if (0<dt && dt<h) {
			peakVm = _prevVm + dt*_prevVmDot;
		}
	}

	// Save the current vm and vmdot for the next step
	_prevVm = vm;
	_prevVmDot = vmDot;

	// Return the result
	return peakVm;
}

/**********************NOT CURRENTLY USED****************************** 

// Get probability of LTP inducing event given expected number of NMDAR open
Number NMDARDepPlasticityRule::probLTP(Number nopen)
{
	// NMDA receptor openings are assumed to be independent and
	// longer than the width of the backpropagating action potential.
	// A Poisson distribution with a mean of nopen is used to determine
	// the probability that the number open at the BPAP exceeds the
	// threshold for inducing LTP.

	int		k;
	double	pltp = exp(-nopen);			// poisson probability starting with 0 open
	double	qltp = 1 - pltp;			// poisson residual quantile (start at >0 open)
	
	for (k=1;k<nLTP();k++) {			// sum probabilities up to threshold
		pltp *= nopen/k;				// probability of k open
		qltp -= pltp;					// residual quantile for >k open
	}
	return (Number) qltp;				// all done - round to Number precision
}

**********************************************************************/

// Get probability of LTP inducing event given expected number of NMDAR open
Number NMDARDepPlasticityRule::probLTP(Number nopen)
{
	// Assign probLTP based on fourth order kinetic formula

	// When nopen is exactly zero (however unlikely) return 0.
	if (nopen==0) {
		return 0;
	}

	// Otherwise, use nLTP as the half activation value and
	// apply a fourth order kinetic formulation.
	Number	x = nLTP() / nopen;

	return 1/(1+x*x*x*x);
}


// Get probability of LTD inducing event given postSpikeTime-PreSpikeTime
Number NMDARDepPlasticityRule::probLTD(SimTime dt)
{
	// Get probability of LTD based on a rough fit to the
	// LTD part of the STDP curve from Debanne (similar to Feldman).
	// PltdMin handles the case of unpaired homogenous LTD, that is,
	// LTD caused by repeated presynaptic stimulation without
	// any paired postsynaptic activity. PltdMax results from
	// simultaneous pre- and postsynaptic pairings.

	// Handle special cases both for efficiency and to prevent mishandling
	// of exp(-infinity), which functions differently in release versus
	// debug mode when using MS compilers. 
	
	// dt>=0 can happen only when both pre and post synaptic action potentials
	// fall into the same time step or when a simultaneous pairing occurs.
	// In this case, a previous homosynaptic LTD event is assumed for these
	// types of pairings and is offset against the resulting paired LTD rate.
	// The probability of LTD is prorated over the Topen interval so that it
	// is maximum at dt=0 and zero at dt=Topen;

	if (dt>=0) {
		if (dt<Topen() ) { 
			return (PltdMax()-PltdMin()) * (1-dt/Topen());
		}
		else {
			return 0;
		}
	}

	if (dt==InfinitePast)
		return PltdMin();

	// Handle the remaining case of dt>0 but less than infinity.
	// For dt near 0, this probably overstates LTD but there is
	// little data from which to make a better model.
	return PltdMin()+(PltdMax()-PltdMin())*exp(dt/PltdTau());
}
