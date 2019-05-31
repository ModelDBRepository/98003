// Mammalian CA3 cell model based on reconstructed cell morphology
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: neuron_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file contains the class bodies used to implement
// a realistic CA3 Pyramidal cell model with channel properties
// taken from recent and hopefully updated measurements.
//
// References (see header file)

#include "ionchan_na_baker_2003.h"
#include "ionchan_ca_baker_2003.h"
#include "ionchan_k_dr_baker_2003.h"
#include "ionchan_k_c_baker_2003.h"
#include "ionchan_k_ahp_baker_2003.h"
#include "ionchan_k_a_baker_2003.h"
#include "ionchan_k_m_baker_2003.h"
#include "ionchan_ih_baker_2003.h"
#include "synapse_glu_baker_2003.h"
#include "plasticity_glu_baker_2003.h"
#include "synapse_gaba_baker_2003.h"
#include "neuron_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;

// Specify default namespace(s)
using namespace BAKER_2003;


// -----------------------------------------------
// PyramidalCellNeuron class body
// -----------------------------------------------

// Constructors and a destructor
PyramidalCell::PyramidalCell(Model* m) : MorphologicalNeuron(m) 
{
	// Initialize state values
	_AChLevel = 0;
	_NaSGateDisabled = false;
	_AChInmdarModDisabled = false;

	// Set a default orientation vector
	_orientationX = 0;
	_orientationY = 1;
	_orientationZ = 0;
}

PyramidalCell::~PyramidalCell() {}

// Create a default ODE solver for this neuron.
ODESolver* PyramidalCell::defaultSolver()
{
	NeuronSolver* mySolver = new NeuronSolver;

	mySolver->timeStep(1.0*msec);
	mySolver->errTol(0.001f);
	return mySolver;
}

// Initialize the cell by setting parameters and adding compartments
void PyramidalCell::createCell()
{
	// Set the morphology table that controls cell construction
	setMorphology();

	// Set all common parameters (could be morphology dependent)
	setParameters();

	// Build the various compartments of the cell
	createSoma();

	// Build the axon and IS compartments
	createAxonAndIS();

	// Build the dendrite compartments
	createDendrites();

	// Ensure (re)setting values based on AChLevel
	// if is set to other than the default value.
	if (AChLevel()>0) {
		AChLevel( AChLevel() );
	}
}

// Build the soma compartment
void PyramidalCell::createSoma()
{
	// Working variables
	Compartment*		soma;
	SimpleCalciumPool*	pool;
	SimpleCalciumPool*	mdp;
	IonChannel*			CaTChan;
	IonChannel*			CaNChan;
	IonChannel*			CaLChan;
	Number				area,radius;

	// Explicitly create a compartment for the soma
	// The soma is treated as an equipotential ball with a
	// radius chosen to preserve membrane area.
	area=2*Pi*_morphology[0].r*_morphology[0].len*micron_2;
	radius=sqrt(area/(4*Pi));

	soma = new SphericalCompartment(radius,_Cm,_Rm);
	soma->componentName("Soma");
	soma->numericIdentifier(0);

	// Make calcium pools for the soma.
	// pool is a pool for all calcium while mdp is a 
	// micro-domain pool used for Ca++ activated channels.
	pool = new SimpleCalciumPool(
				_CaXrest,_CaXinit,
				_CaUbr,_CaVmax,_CaKd,
				_somaShellDepth);

	mdp  = new SimpleCalciumPool(
				_CaXrest,_CaXinit,
				_CaUbr,_CaVmax,_CaKd,
				_mdShellDepth);
	mdp->componentName("CaMDP");

	// Set effective leakage reversal potential
	soma->Vleak(_Vleak);
	soma->Vinit(_Vinit);
	
	// Populate the soma with channels
	soma->add(new Soma_Na_channel(_gNaTSoma,NaSGateDisabled(),NaSGateDisabledValue() ));
	soma->add(new Persistent_Na_channel(_gNaPSoma));
	soma->add(new Proximal_Ih_channel(_gIh));
	soma->add(new K_DR_channel(_gKdr));
	soma->add(new Proximal_K_A_channel(_gKa));
	soma->add(new K_M_channel(_gKm));
	soma->add(new K_C_channel(_gKc,mdp));
	soma->add(new K_AHP_channel(_gKahp,mdp));
	soma->add(CaTChan = new Ca_T_channel(_gCaTSoma));
	soma->add(CaNChan = new Ca_N_channel(_gCaN));
	soma->add(CaLChan = new Ca_L_channel(_gCaL));

	// Selectively associate calcium sources with the micro-domain pool
	mdp->addSourceChannel(CaNChan);
	mdp->addSourceChannel(CaLChan);

	// Add synaptic conductances
	soma->add(new GABAa_SynapticResp(_gGABAa));

	// Add the calcium pool after all ion channels
	soma->add(pool);
	soma->add(mdp);

	// Make the soma is part of the cell model
	add(soma);
	_numSomaComp = 1;

	// Optionally, debug the calcium pool values
	if (0) {
		cerr<<"Soma calcium pool beta = "<<pool->beta()
			<<" phi="<<pool->phi()<<endl;
		cerr<<"Soma microdomain pool beta = "<<mdp->beta()
			<<" phi="<<mdp->phi()<<endl;
	}
}

// Build the IS and axon compartments
void PyramidalCell::createAxonAndIS()
{
	// Working variables
	Compartment*		pcomp;
	Compartment*		initialSegment;

	int					k;
	char				nameBuf[32];
	Number				radius;
	Number				taperRate;
	Number				dist;

	// Add the initial segment and populate with channels
	_ISCompOffset = _compartments.size();
	_numISComp = 1;

	initialSegment = new Compartment(
		_initSegRadius,
		_initSegLen,
		_Cm,
		_Ri,
		_Rm,
		_areaAdjustment);
	initialSegment->componentName("IS");
	initialSegment->numericIdentifier(10000);
	initialSegment->distFromSoma(_initSegLen/2);

	initialSegment->add(new Soma_Na_channel(
		_gNaTSoma,
		NaSGateDisabled(),
		NaSGateDisabledValue() ));

	initialSegment->add(new K_DR_channel(_gKdr));
	initialSegment->add(new K_M_channel(_gKm));
	initialSegment->add(new GABAa_SynapticResp(_gGABAa));
	initialSegment->Vleak(_Vleak);
	initialSegment->Vinit(_Vinit);
	add(initialSegment);

	// Add axon segments and populate with channels
	_axonCompOffset = _compartments.size();
	radius = _axonProximalRadius;
	taperRate = pow(
		double(_axonDistalRadius/_axonProximalRadius),
		1.0/(_numAxonComp-1));
	dist = _initSegLen-_axonCompLen/2;

	for (k=0; k<_numAxonComp; k++,radius*=taperRate) {

		// Get distance to center of compartment
		dist += _axonCompLen;

		// Create the compartment and set its params
		pcomp = new Compartment(
			radius,
			_axonCompLen,
			_CmAxon,
			_RiAxon,
			_RmAxon);
		sprintf(nameBuf,"Axon%d",k+1);
		pcomp->componentName(nameBuf);
		pcomp->numericIdentifier(10000+k);
		pcomp->distFromSoma(dist);

		pcomp->add( new Axon_Na_channel(_gNaAxon));
		pcomp->add( new K_DR_channel(_gKdrAxon));

		pcomp->Vinit(_Vinit + _VinitAxonSlope*dist);
		pcomp->Vleak(_VleakAxon);
		add(pcomp);

		// Connect axon segments with each other.
		// Note that axons are numbered starting with 1.
		if (k>0) {
			new ElectricalCoupling( axonComp(k),axonComp(k+1) );
		}
	}
	
	// Wire together soma, initial segment and axon segments
	new ElectricalCoupling( initialSegment, axonComp(1) );
	new ElectricalCoupling( somaComp(), initialSegment );
}

// Build the dendrite compartments
void PyramidalCell::createDendrites()
{
	// Declare a dummy channel object just to get the number of
	// state variables. This will be used later to fill in a 
	// place holder in compartments without such channels.
	K_C_channel			dummyKC;
	const int			numKCVar = dummyKC.numStateVar();

	// Working variables
	Compartment*		pcomp;
	SimpleCalciumPool*	pool;
	SimpleCalciumPool*	mdp;
	IonChannel*			CaTChan;
	IonChannel*			CaNChan;
	IonChannel*			CaLChan;

	int					k;
	char				nameBuf[32];
	Number				dist,densityDist;
	Number				compProxDist,compDistalDist,limitDist;
	Number				len,radius;
	Number				blendRatio;
	Number				gCalc;
	Number				xCaDepK;

	// Set starting offset and number of comp.
	// All entries in morphology table are dendrites except first.
	_dendriteCompOffset = _compartments.size();
	_numDendriteComp = _morphologySize-1;

	// Populate dendrites with channels and calcium pools
	for (k=1; k<_morphologySize; k++) {

		// Determine the effective distance used for computing
		// channel densities based on distance.
		dist = _morphology[k].dist*micron;			// distance to comp middle point
		densityDist = minval(dist,_maxDensityDist);	// distance for calculations
		
		// Create a new dendrite compartment with adjustments to Cm and Rm for spines
		pcomp = new Compartment(
				_morphology[k].r * micron,			// radius
				_morphology[k].len * micron,		// length
				_Cm,								// membrane capacitance
				_Ri,								// axial resistance
				_Rm/(1+_gLeakSlope*densityDist),	// membrane resistance
				_areaAdjustment );					// adjustment for spines etc.
		pcomp->distFromSoma(dist);					// set distance from soma

		// Compute effective blend ratio as the average ratio
		// thoughout the compartment. Blending is assumed to
		// be a linear process starting with a ratio value of 0
		// at the soma and increasing linearly up to a value of 1
		// at maxBlendDist from the soma.
		len = _morphology[k].len*micron;	// compartment length
		radius = _morphology[k].r*micron;	// compartment radius
		compDistalDist = dist+len/2;		// distance to distal end
		compProxDist = plusval(dist-len/2);	// distance to proximal end

		// Test if compartment lies beyond blend region
		if (_maxBlendDist<=compProxDist) {
			blendRatio = 1.0;
		}
		// Otherwise get location when blend region ends in compartment
		else {
			limitDist = minval(compDistalDist,_maxBlendDist);

			// Prorate the blending over the part that lies in the
			// blend region with any part extending past it.
			blendRatio =
				(limitDist+compProxDist)/2/_maxBlendDist	// blended part
					*(limitDist-compProxDist)/len		// weighted by length
				+ (compDistalDist-limitDist)/len; 		// distal part (ratio=1)
		}

		// Get the distance dependent conductance ratio for Ca++
		// dependent K+ channels. These are treated as conductances 
		// decreasing to 0 at maxCaDepKDist from soma. Because this may 
		// lie within a compartment, it is necessary to prorate the
		// conductance over the compartment by considering different cases.
		if (compProxDist>_maxCaDepKDist) {
			xCaDepK = 0;
		}
		else {
			if (compDistalDist<_maxCaDepKDist) {
				// Entire compartment lies within maxKcDist limit.
				// Use midpoint distance to get average in compartment.
				xCaDepK = (1-dist/_maxCaDepKDist);
			}
			else {
				// xCaDepK>0 is over length maxCaDepKDist-compProxDist
				// with maximum at proximal end and a minimum of 0 at distal, 
				// Prorate xCaDepK as half the maximum value over the fraction of
				// the compartment where the conductance is non-zero.
				xCaDepK = (1-compProxDist/_maxCaDepKDist)/2
						* (_maxCaDepKDist-compProxDist)/len;			
			}
		}

		// Assign a name and id to the compartment for reporting
		sprintf(nameBuf,"D%04d",k);
		pcomp->componentName(nameBuf);
		pcomp->numericIdentifier(k);

		// Set the leakage reversal potential and initial potential
		pcomp->Vleak(_Vleak);
		pcomp->Vinit(_Vinit+dist*_VinitDendriteSlope);

		// Make calcium pools for the dendrite.
		// pool is a pool for all calcium while mdp is a 
		// micro-domain pool used for Ca++ activated channels.
		pool = new SimpleCalciumPool(
					_CaXrest,_CaXinit,
					_CaUbr,_CaVmax,_CaKd);

		mdp  = new SimpleCalciumPool(
					_CaXrest,_CaXinit,
					_CaUbr,_CaVmax,_CaKd,
					_mdShellDepth);

		pool->componentName("CaPool");
		mdp->componentName("CaMDP");

		// Populate the compartment with ion channels.

		// Add transient Na channels. Conductance follows blend ratio
		// in terms of proximal vs distal channels. Where possible,
		// blended channels are not used for performance reasons.
		gCalc = (1-blendRatio)*_gNaTProximal + blendRatio*_gNaTDistal;
		if (blendRatio==0) {
			pcomp->add(new Proximal_Na_channel(gCalc,
							NaSGateDisabled(), 
							NaSGateDisabledValue() ));
		}
		else if (blendRatio<1) {
			pcomp->add(new Blended_Na_channel(gCalc,blendRatio,
							NaSGateDisabled(),
							NaSGateDisabledValue() ));
		}
		else {
			pcomp->add(new Distal_Na_channel(gCalc,
							NaSGateDisabled(), 
							NaSGateDisabledValue() ));
		}

		// Add persistent Na channels
		gCalc = (1-blendRatio)*_gNaPProximal + blendRatio*_gNaPDistal;
		pcomp->add(new Persistent_Na_channel(gCalc));
		
		// Add Ih channels
		gCalc=_gIh*(1+_gIhSlope*densityDist);
		if (blendRatio==0) 
			pcomp->add(new Proximal_Ih_channel(gCalc) );
		else if (blendRatio>=1)
			pcomp->add(new Distal_Ih_channel(gCalc) );
		else 
			pcomp->add(new Blended_Ih_channel(gCalc,blendRatio) );
		
		// Add Kdr channels
		pcomp->add(new K_DR_channel(_gKdr*(1+_gKdrSlope*densityDist)));
		
		// Add K-A channels
		gCalc=_gKa*(1+_gKaSlope*densityDist);
		if (radius<_obliqueRadius) {
			gCalc *= pow(radius/_obliqueRadius,_obliqueKaPower);
		}
		if (blendRatio==0)
			pcomp->add(new Proximal_K_A_channel(gCalc) );
		else if (blendRatio>=1)
			pcomp->add(new Distal_K_A_channel(gCalc) );
		else
			pcomp->add(new Blended_K_A_channel(gCalc,blendRatio) );

		// Add K-M channels
		pcomp->add(new K_M_channel(_gKm));

		// Add K-C channels. These are Ca++ dependent channels.
		// If conductance is zero reserve state vector entry only.
		if (xCaDepK>0) {
			pcomp->add(new K_C_channel(gCalc*xCaDepK,mdp));
		}
		else {
			pcomp->add(new DummyIonChannel(numKCVar,mdp));
		}

		// Add other channels
		pcomp->add(new K_AHP_channel(_gKahp,mdp));
		pcomp->add(CaTChan = new Ca_T_channel(_gCaTDendrite));
		pcomp->add(CaNChan = new Ca_N_channel(_gCaN));
		pcomp->add(CaLChan = new Ca_L_channel(_gCaL));
		
		// Selectively associate calcium sources with microdomain pool
		mdp->addSourceChannel(CaNChan);
		mdp->addSourceChannel(CaLChan);

		// Add synaptic conductances. All possible synapse types are
		// included even if this particular compartment might not have
		// any synapses of the associated type.
		pcomp->add(new AC_Glu_SynapticResp(_gAMPA_AC,_gNR2A_AC));
		pcomp->add(new PP_Glu_SynapticResp(_gAMPA_PP,_gNR2A_PP));
		pcomp->add(new MF_Glu_SynapticResp(_gAMPA_MF,_gNR2A_MF));
		pcomp->add(new GABAa_SynapticResp(_gGABAa));
		pcomp->add(new GABAas_SynapticResp(_gGABAas));
		pcomp->add(new GABAb_SynapticResp(_gGABAb));

		// Add the calcium pools after all ion channels
		pcomp->add(pool);
		pcomp->add(mdp);

		// Add the compartment to the neuron's collection
		add(pcomp);

		// Optionally, debug the calcium pool values
		if (0) {
			cerr<<"Dendrite comp  "<<k<<" calcium pool beta = "<<pool->beta()
				<<" phi="<<pool->phi()<<endl;
			cerr<<"Dendrite comp  "<<k<<" microdomain pool beta = "<<mdp->beta()
				<<" phi="<<pool->phi()<<endl;
		}

	}

	// Electrically connect compartments with each other and the soma.
	connectDendrites();
}

// Accessor for controlling ACh modulation.
void PyramidalCell::AChLevel(Number AChExt)
{
	// All this obviously does not belong in the Neuron
	// hierarchy, but neuromodulation as defined here is 
	// too specific (and simplistic) to include as a general
	// IonChannel behavior. A policy object is probably a better 
	// design, but for now a short-cut is used.

	static const TokenId	AChMod	= token("AChModulator");

	static const TokenId	Ih		= token("Ih");
	static const TokenId	NaT		= token("Na");
	static const TokenId	NaP		= token("NaP");
	static const TokenId	Ka		= token("K_A");
	static const TokenId	Km		= token("K_M");
	static const TokenId	Kahp	= token("K_AHP");
	static const TokenId	CaN		= token("CaN");
	static const TokenId	CaL		= token("CaL");
	static const TokenId	CaT		= token("CaT");

	Number					vinit,vleak,gleak;
	Compartment*			comp;
	int						k;

	// Save the value supplied
	_AChLevel = AChExt;

	// If createCell is not done yet, skip the rest
	if (somaComp()==NULL)
		return;

	// Get the effective membrane resistance and reversal potential
	// by combining K leak with a mixed cation leak. Note that as a
	// simplification of experimental results, the net membrane 
	// resistance (Rm) is left unaffected by ACh modulation. Only
	// the leakage reversal potential is thus affected. Results from
	// Seeger & Alzheimer suggest that the effect is perisomatic in
	// CA1 pyramidal cells, but others (McBain etc.) have suggested
	// that cation currents are involved in dendritic oscillations.

	if (AChLevel() > 0) {
		gleak = 1/_Rm - _gLeakAChModMC;
		vleak = _Rm*(_Vleak*gleak + _VleakAChModMC*_gLeakAChModMC);
	}
	else {
		vleak = _Vleak;
	}

	// Get the new soma Vinit value
	vinit = _Vinit + (AChLevel()>0 ? _VinitAChAdjust : 0);

	// Set leak parameters..
	somaComp()->Vinit(vinit);
	somaComp()->Vleak(vleak);
	initialSegment()->Vinit(vinit);
	initialSegment()->Vleak(vleak);

	// Apply leak parameters to each dendrite compartment
	for (k=1;k<=numDendriteComp();k++) {
		comp = dendriteComp(k);
		comp->Vinit(vinit + _VinitDendriteSlope * comp->distFromSoma() );
		comp->Vleak(vleak);
	}

	// Apply leak parameters to each axon compartment.
	// vleak is not affected for axons.
	for (k=1;k<=numAxonComp();k++) {
		comp = axonComp(k);
		comp->Vinit(vinit + _VinitAxonSlope * comp->distFromSoma() );
	}

	// Apply neuromodulation factors for each channel
	setMichaelisMentenMod(NaT,	AChLevel(),	_NaT_AChA,		_NaT_AChKd);
	setMichaelisMentenMod(NaP,	AChLevel(),	_NaP_AChA,		_NaP_AChKd);
	setMichaelisMentenMod(Ih,	AChLevel(),	_Ih_AChA,		_Ih_AChKd);
	setMichaelisMentenMod(Ka,	AChLevel(),	_Ka_AChA,		_Ka_AChKd);
	setMichaelisMentenMod(Km,	AChLevel(),	_Km_AChA,		_Km_AChKd);
	setMichaelisMentenMod(Kahp,	AChLevel(),	_Kahp_AChA,		_Kahp_AChKd);
	setMichaelisMentenMod(CaN,	AChLevel(),	_CaN_AChA,		_CaN_AChKd);
	setMichaelisMentenMod(CaL,	AChLevel(),	_CaL_AChA,		_CaL_AChKd);
	setMichaelisMentenMod(CaT,	AChLevel(),	_CaT_AChA,		_CaT_AChKd);

	// Set the s-gate status taking the new AChLevel into account.
	NaSGateDisabled( AChLevel()>=NaSGateDisabledACh() );

	// Set synapse modulation
	setGABAMod();
	setGluMod();
}

// Set ACh Level for GABA synapses
void PyramidalCell::setGABAMod()
{

	static const TokenId	AChMod	= token("AChModulator");

	static const TokenId	GABAa	= token("GABAa");
	static const TokenId	GABAas	= token("GABAas");
	static const TokenId	GABAb	= token("GABAb");

	static const int		numGABAIds = 3;		
	static const TokenId	GABASynapseIds[numGABAIds] = { GABAa, GABAas, GABAb };

	int						i,k;
	Compartment*			comp;

	// Apply ACh neuromodulation to soma, initial segment, and dendrites.
	for (k=-1; k<=numDendriteComp();k++) {

		if		(k ==-1)	comp = ISComp();
		else if (k == 0)	comp = somaComp();
		else				comp = dendriteComp(k);

		for (i=0;i<numGABAIds;i++) {
			SynapticResponse* resp;

			resp = comp->findSynapticResponse(GABASynapseIds[i],false);
			if (resp!=NULL) {
				resp->setModParam(AChMod,AChLevel());
			}
		}
	}

}

// Set ACh level and modulation for Glu synapses
void PyramidalCell::setGluMod()
{

	static const TokenId	AChMod	= token("AChModulator");

	static const TokenId	ACGlu	= token("AC_GluR");
	static const TokenId	PPGlu	= token("PP_GluR");
	static const TokenId	MFGlu	= token("MF_GluR");

	static const int		numGluIds = 3;		
	static const TokenId	GluSynapseIds[numGluIds] = { MFGlu, ACGlu, PPGlu };

	int						i,k;
	Compartment*			comp;

	// Apply ACh neuromodulation to synapses affected.
	for (k=1; k<=numDendriteComp();k++) {

		comp = dendriteComp(k);

		// Pass ACh level to each synaptic response type
		for (i=0;i<numGluIds; i++) {	
			SynapticResponse* resp;

			resp = comp->findSynapticResponse(GluSynapseIds[i],false);
			if (resp!=NULL) {
				resp->setModParam(AChMod,AChLevel());
			}
		}
	}
}

// Knockout NMDA receptor effects selectively
void PyramidalCell::knockoutNMDAR(
	bool	disableCurrent,
	bool	disablePlasticity,
	bool	disableAChMod,
	bool	disableCaDepSupp)
{
	static const TokenId	ACGlu	= token("AC_GluR");
	static const TokenId	PPGlu	= token("PP_GluR");
	static const TokenId	MFGlu	= token("MF_GluR");
	static const TokenId	NR2A	= token("NR2A");
	static const TokenId	AMPAR	= token("AMPAR");
	static const TokenId	STDP	= token("NMDARDepPlasticityRule");

	static const int		numNR2AIds = 3;		
	static const TokenId	NR2ASynapseIds[numNR2AIds] = { MFGlu, ACGlu, PPGlu };
	const Number			NR2AGMax[numNR2AIds] = { _gNR2A_MF, _gNR2A_AC, _gNR2A_PP };

	int i,k;


	// Set NMDAR conductance depending on AChMod (if changed)
	if (_AChInmdarModDisabled != disableAChMod) {
		setGluMod();
	}
	_AChInmdarModDisabled = disableAChMod;

	// Set NMDAR properties based on knockout effects selected
	for (k=1; k<=numDendriteComp();k++) {

		Compartment* comp = dendriteComp(k);

		// Apply to each synaptic response type containing NR2A
		for (i=0;i<numNR2AIds; i++) {
			
			SynapticResponse* resp;
			SynapticConductance* nr2achan;
			SynapticConductance* ampar;

			resp = comp->findSynapticResponse(NR2ASynapseIds[i],false);
			if (resp==NULL)
				continue;
				
			// Set the NMDAR channel conductance
			nr2achan = (SynapticConductance*) resp->findIonChannel(NR2A);
			nr2achan->gMax(disableCurrent ? 0 : NR2AGMax[i]);

			// Set AMPAR NMDAR-dependent plasticity, but not for MFGlu
			// Also set other parameters specific to NMDAR plasticity
			if (i>0) {
				ampar = (SynapticConductance*) resp->findIonChannel(AMPAR);

				// Make sure the right type of plasticity is provided.
				// At present there is only one valid option.
				if (ampar->postsynapticRule()->plasticityTypeId()==STDP) {

					NMDARDepPlasticityRule* nmdarRule
						= (NMDARDepPlasticityRule*) (ampar->postsynapticRule());

					nmdarRule->isDisabled( disablePlasticity );
					nmdarRule->disableAChMod( disableAChMod );
					nmdarRule->disableCaDepSupp( disableCaDepSupp );
				}
				else {
					FatalError("PyramidalCell::knockoutNMDAR) "
						"Unknown plasticity type found");
				}
			}
		}
	}
}

// Accessor for controlling Na s-gate status
void PyramidalCell::NaSGateDisabled(bool disableSGate)
{	
	static TokenId	Na = token("Na");
	int				k;

	Abstract_Na_channel* na;

	// Save the value supplied
	_NaSGateDisabled = disableSGate;

	// If createCell is not done yet, skip the rest
	if (somaComp()==NULL)
		return;

	// Apply to Soma Na channel
	na=(Abstract_Na_channel*) (somaComp()->findIonChannel(Na));
	na->sGateDisabled( disableSGate );
	na->sGateDisabledValue( NaSGateDisabledValue() );

	// Apply to all dendrites
	for (k=1;k<=numDendriteComp();k++) {
		na = (Abstract_Na_channel*) (dendriteComp(k)->findIonChannel(Na));
		na->sGateDisabled( disableSGate );
		na->sGateDisabledValue( NaSGateDisabledValue() );
	}
}



// -----------------------------------------------
// CA3 PyramidalCellNeuron class body
// -----------------------------------------------



// Constructors and a destructor
CA3PyramidalCell::CA3PyramidalCell(Model* m, bool doInit) 
: PyramidalCell(m) 
{
	// Build all components of the cell unless this
	// is to be deferred to the subclass constructor.
	if (doInit) {
		createCell();
	}
}

CA3PyramidalCell::~CA3PyramidalCell() {}

// Set the parameters controlling cell construction
void CA3PyramidalCell::setParameters()
{
	// Passive properties ---------------------------------------------

	_Cm					= 1.0*microF/cm_2;	// membrance capacitance
	_Ri					= 100*ohm*cm;		// axial resistivity
	_Rm					= 200*kohm*cm_2;	// membrane resistivity (leak)
	_Vleak				= -90*mV;			// Leak reversal potential (K+)
	_VleakAChModMC		= 0*mV;				// non-selective mixed-cation reversal
	_gLeakAChModMC		= 1.5*microS/cm_2;	// mixed-cation channels modulated by ACh

	_CmAxon				= 1.0*microF/cm_2;	// capacitance for axon only
	_RiAxon				= 100*ohm*cm;		// axial resistivity for axon
	_RmAxon				= 750*ohm*cm_2;		// membrane resistivity for axon
	_VleakAxon			= -74*mV;			// leak reversal for axon (near local rest)

	// Factor to increase compartment membrane area to allow for spines
	// and other irregularities of structure with respect to a perfect
	// cylinder. The value use here is an estimate only. Perhaps a better
	// solution would be to associate some amount of membrane area with 
	// spines once synapses are added. This is left for the future.
	_areaAdjustment		= 1.25f;

	// Parameters covering IS and axon geometry. This is roughly the
	// axon found in cell L56a in the Duke-Southampton archive.
	// Length values are per compartment.Radius of the axon tapers so
	// the ratio of radius from one compartment to the next is constant.
	_numAxonComp		= 10;
	_initSegLen			= 30*micron;
	_axonCompLen		= 30*micron;

	_initSegRadius		= 0.7*micron;
	_axonProximalRadius	= 0.5*micron;
	_axonDistalRadius	= 0.2*micron;

	// Calcium buffer parameters --------------------------------------

	// See Jaffe et al. for parameter background.
	// Pool dynamics are derived from Sabatini et al. table 1.
	// and also Maravall et al. (referenced in Sabatini).
	_CaXinit			= 50*nanoM;		// initial [Ca++]-in
	_CaXrest			= 50*nanoM;		// resting [Ca++]-in for pump
	_CaUbr				= 0.013f;		// Unbound Ca++ ion ratio (Ke=75)
	_CaKd				= 2*microM;		// half activation for Ca++ pump
	_CaVmax				= 1.5e-12*mMole/msec/cm_2; // peak Ca++ pump rate
	_mdShellDepth		= 0.1*micron;	// microdomain subshell depth
	_somaShellDepth		= 1*micron;		// soma pool shell depth

	// Ion channel parameters -----------------------------------------

	// Set axon conductances based on whole cell response.
	// See Colbert and Pan for ratio between axon and soma Na currents
	// as measured in neocortical pyramidal cells (2-3 times soma).
	// In axons there is less active repolarization and no K-A. gKdr 
	// is lower than otherwise expected while membrane conductance 
	// has been increased to compensate. Kdr in axons may have somewhat
	// different properties than modeled here. Axon Rm and gKdr are chosen
	// to permit repolarization of the axon during a burst response.
	_gNaAxon			= 60*mS/cm_2;
	_gKdrAxon			= 10*mS/cm_2;

	// Soma and dendrite ion channel conductances. These are set to
	// represent physiological conditions (e.g. 37-38 deg C).
	_gNaTSoma			= 30.0*mS/cm_2;	// Density for soma	
	_gNaTProximal		= 40.0*mS/cm_2;	// Density for proximal dendrites
	_gNaTDistal			= 10.0*mS/cm_2;	// Density for distal dendrites

	// ACh disables Na slow inactivation (Tsubokawa & Ross, 1997). 
	// This parameter provides the fixed value for NaT channel 
	// s-gate when the gate is disabled. See also Cantrell et al. 
	// for ACh effects on CA1 Na+ currents. A decrease in Na+ current 
	// is expected but the distribution  along the somato-dendritic axis 
	// is unknown. Dose-response is derived from Cantrell et al., but
	// variability from cell to cell suggests a range of values.
	_NaSGateDisabledACh	= 0.1*microM;	// arbitrary threshold
	_NaSGateDisabledValue = 0.6f;		// approx value at rest
	_NaT_AChA			= -0.35f;		// -29.8% for [Cb]=50 microM
	_NaT_AChKd			= 95*microM;	// assuming 10:1 for Cb:ACh.

	// NaP affects effective resting membrane resistance, ADP, and
	// repetitive spiking. ACh down regulation is from Alroy et al.
	// See also Mittmann & Alzheimer (1998). NaP is suppressed to 
	// inhibit bursting wnen ACh>0. AChKd is otherwise arbitrary.
	// French et al. concluded that NaP was primarily near the soma 
	// since severing dendrites had minimal effect at the soma.
	// Cantrell et al. suggest that the suppression might not be
	// total, but this is necessary to prevent bursting when [ACh]>0.
	_gNaPSoma			= 0.2*mS/cm_2;	// 1-2% peak Ina
	_gNaPProximal		= 0.2*mS/cm_2;	// 1-2% peak Ina
	_gNaPDistal			= 0.0*mS/cm_2;	// see note above
	_NaP_AChA			= -1.0f;		// total suppression
	_NaP_AChKd			= 0.1*microM;	// arbitrary (value unknown)

	// Ih currents are set based on resting potential and R-in.
	// As noted by Spruston & Johnston, Ih is less in CA3 than CA1.
	// Halliwell and Adams found no effect of muscarine on Ih.
	// Fishan et al. measuremented ACh effects in CA3 but results
	// are somewhat unclear. Parameters here are based on the
	// change in charge transfer resulting from 20 microM muscarine.
	_gIh				= 0.004*mS/cm_2;
	_Ih_AChA			= 0.4f;			// from charge transfer ratio
	_Ih_AChKd			= 1*microM;		// arbitrary (value unknown)

	// Kdr and Ka are set to permit repolarization and burst response.
	// See Klee et al. Table 1 for K-A/Kdr current ratios during
	// development. The ratio here follows the low Ka pattern of P26.
	// Nakajima et al. found an ACh effect in CA1 via changes in
	// voltage sensitivity. This is only roughly approximated here.
	_gKdr				= 3.0*mS/cm_2;
	_gKa				= 1.2*mS/cm_2;
	_Ka_AChA			= -0.3f;
	_Ka_AChKd			= 0.1*microM;

	// Kc is taken as 20% of the peak sustained K+ conductance following
	// Klee et al. Because the density of Kc decreases with distance from
	// the soma, this is a rough approximation.
	_gKc				= 0.8*mS/cm_2;

	// Km is component of the medium AHP. It also affects resting input
	// resistance and resists accommodation. Kahp is set for a slower AHP 
	// current following a burst.Kahp could be modulated directly by ACh
	// or indirectly by reduction of Ca++ currents. Reduction of Ca++
	// is sufficient to inhibit slow AHP when [ACh]>0.
	_gKm				= 0.3*mS/cm_2;
	_Km_AChA			= -1;
	_Km_AChKd			= 64*microM;

	_gKahp				= 0.5*mS/cm_2;
	_Kahp_AChA			= 0;			// effect via Ca++ reduction only
	_Kahp_AChKd			= 3*microM;

	// Ca++ R/N and L channels conductances were estimated by Magee & Johnston
	// for CA1 pyramidal cells using barium salts. To compensate for the 
	// difference between Ba++ and Ca++, conductances are reduced by a 
	// factor of 2.5 to 3 (Johnston & Wu 1995). To compensate for the 
	// difference between room temp and body temp, the conductance is 
	// doubled, (Brown et al.) leaving us roughly back where we started.
	// See Menschik for ACh modulation as well as Fischer & Johnston.	
	_gCaN				= 1.5*mS/cm_2;
	_CaN_AChA			= -0.4f;
	_CaN_AChKd			= 1.6*microM;

	_gCaL				= 1.2*mS/cm_2;
	_CaL_AChA			= -0.7f;
	_CaL_AChKd			= 1.7*microM;

	// Ca-T conductances for CA3 pyramidal cell. For CA1, Magee & Johnston
	// estimated 1.0 mS/cm_2. Fisher found that CA3 patches have 2 to 3 
	// times as many Ca-T channels as CA1 patches and that Ca-T channels
	// are found in the soma as well as dendrites. For these channels, 
	// Ba++ has a relative permissivity of 0.85 (Takahashi, 1991). Body 
	// versus room temp could also affect currents by a factor of two. 
	// ACh up-regulates these currents but also stimulates interneuros 
	// resulting in GABA release which down-regulates the same channels. 
	// Nor-adrenergic moduluation also down-regulates Ca-T channels 
	// (Fischer & Johnston) though this is not modelled here. The
	// conductance and ACh modulation are obviously estimates only.

	_gCaTSoma			= 2.0*mS/cm_2;
	_gCaTDendrite		= 2.0*mS/cm_2;
	_CaT_AChA			= 1.0f;			// see Fisher & Johnston
	_CaT_AChKd			= 1.7*microM;	// estimate only (from Ca-L)

	// Location dependency parameters ---------------------------------

	// Distance dependent ion channel conductance parameters
	_maxDensityDist		= 9999*micron;			// max dist for setting conductances
	_maxBlendDist		= 150*micron;			// max dist for blended channels
	_maxCaDepKDist		= 150*micron;			// max dist for proximal Kc

	_gLeakSlope			= 0.0/(350*micron);		// change in leak density
	_gKdrSlope			= 0.14/(350*micron);	// change in Kdr
	_gKaSlope			= 5.2/(350*micron);		// change in K-A
	_gIhSlope			= 5.8/(350*micron);		// change in Ih

	// Set properties for oblique (terminal) dendrites. 
	// Evidence for reduced excitability in obliques is somewhat 
	// indirect as in Frick et al. The formula for conductance is:
	// g_oblique = g_normal*(radius/obliqueRadius)^power
	_obliqueRadius		=  0.4*micron;
	_obliqueKaPower		= -1.5;

	// Synaptic parameters --------------------------------------------

	// Synapse peak conductances per action potential event. 
	// These are conductances per individual synapse.
	// GABAb is set to approximate total charge transfer ratio
	// with GABAa (ca 1:6, Scanziani 2000) at low frequencies.
	// Given model time constants, this works out to a 1:734 ratio
	// in conductance. GABAb is enhanced at higher frequencies
	// by presynaptic plasticity rules. The value of gGABAb is thus
	// not derived from the conductance of a single ion channel.

	_gAMPA_PP			= 1.0 * nanoS;
	_gAMPA_AC			= 1.0 * nanoS;
	_gAMPA_MF			= 2.5 * nanoS;
	_gNR2A_PP			= 0.2 * nanoS;	// 20% gAMPA_PP
	_gNR2A_AC			= 0.2 * nanoS;	// 20% gAMPA_MF
	_gNR2A_MF			= 0.2 * nanoS;	// 20% gAMPA_MF x 50% Zn++ inact
	_gGABAa				= 0.4 * nanoS;	// value for perisomatic inhibition
	_gGABAas			= 0.4 * nanoS;	// use GABAa value as estimate
	_gGABAb				= 0.5 * picoS;	// see note above
	
	// Initialization parameters --------------------------------------

	// The following parameters control voltage initialization.
	// These values are set to correspond with the equilibrium
	// state reached (e.g. after 3-5 sec settling time) and are 
	// determined empirically after other parameters are set.
	_Vinit				= -69.2*mV;				// Initial value at soma (ACh=0)
	_VinitDendriteSlope	=  0.7*mV/(550*micron);	// Change in Vinit wrt distance
	_VinitAxonSlope		= -4*mV/(315*micron);	// Change wrt distance for axon
	_VinitAChAdjust		=  3.8*mV;				// [ACh]>0 depolarization at soma
}



// -----------------------------------------------
// L56aPyramidalCell class body
// -----------------------------------------------



// Constructors and a destructor
L56aPyramidalCell::L56aPyramidalCell(Model* m, bool doInit) 
: CA3PyramidalCell(m,false) 
{
	if (doInit) 
		createCell();
}

L56aPyramidalCell::~L56aPyramidalCell() {}

// Set the morphology table used in building the cell
void L56aPyramidalCell::setMorphology()
{
	// Pick the morphology to use
	morphology(
		cell_l56a_50_micron()
		// cell_l56a_25_micron()
		// cell_l56a_5_micron()
	);

	// Set an orientation vector
	_orientationX = 0;
	_orientationY = 1;
	_orientationZ = 0;
}

