// Neuron for Rallpack Benchmark (rallpack3)
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: neuron_rallpack.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the declarations defining neuron
// objects for the rallpack1 and rallpack3 benchmarks. Even though
// a neuron is not explicitly used in the benchmark, one is needed 
// to contain the necessary compartments.

#include "neuron_rallpack.h"

using namespace std;
using namespace BNSF;

using namespace RALLPACK;


// -----------------------------------------------
// RPNeuron1 class body
// -----------------------------------------------


// Constructors and a destructor
RPNeuron1::RPNeuron1(int nComp, Number totalLen)
{
	using namespace UOM;

	// Constants
	const Number		Rm = 40*kohm*cm_2;
	const Number		Ri = 100*ohm*cm;
	const Number		Cm = 1*microF/cm_2;
	const Number		Vleak = -65*mV;

	const Number		radius = 0.5*micron;
	const Number		len = totalLen/nComp;

	// Working variables
	Compartment*		pcomp;
	char				nameBuf[32];
	int					k;

	// Allocate compartments
	for (k=0; k<nComp; k++) {

		pcomp = new Compartment(radius,len,Cm,Ri,Rm);

		// Number compartments for debug and reporting.
		// Note that numbering starts with 1 not 0.
		sprintf(nameBuf,"C%04d",k+1);
		pcomp->componentName(nameBuf);

		// Set the leak reversal and initial voltage
		pcomp->Vinit(Vleak);
		pcomp->Vleak(Vleak);
		
		// Add the compartment to the neuron's collection
		add(pcomp);

		// Make electrical connections between compartments
		if (k!=0) {
			new ElectricalCoupling(pcomp,_compartments[k-1]);
		}
	}

	// Set compartment distance from soma where the
	// first compartment is the root node (ie soma).
	setDistFromSoma();
}

RPNeuron1::~RPNeuron1() {}

// Create a default ODE solver for this neuron.
// This can be overridden if a different solver is needed.
ODESolver* RPNeuron1::defaultSolver()
{
	return new NeuronSolver;
}




// -----------------------------------------------
// RPNeuron2 class body
// -----------------------------------------------



// Constructors and a destructor
RPNeuron2::RPNeuron2()
{
	using namespace UOM;

	// Constants
	const int			nLevel = 10;
	const Number		Rm = 40*kohm*cm_2;
	const Number		Ri = 100*ohm*cm;
	const Number		Cm = 1*microF/cm_2;
	const Number		Vleak = -65*mV;

	typedef struct {
		int				level;
		int				size;
		double			lenInMicrons;
		double			diamInMicrons; 
	}	AllocEntryType;

		AllocEntryType		allocTable[nLevel] = {

			0,		1,		32.0,	16.0,
			1,		2,		25.4,	10.08,
			2,		4,		20.16,	6.35,
			3,		8,		16.0,	4.0,
			4,		16,		12.7,	2.52,
			5,		32,		10.08,	1.587,
			6,		64,		8.0,	1.0,
			7,		128,	6.35,	0.63,
			8,		256,	5.04,	0.397,
			9,		512,	4.0,	0.25 };

	// Working variables
	Compartment*		pcomp;
	ElectricalJunction*	pjunct;
	char				nameBuf[32];
	int					n,k;
	int					parentIndex = 0;


	// Allocate each of the desired levels in a
	// branching binary tree of compartments.
	for (n=0; n<nLevel; n++) {

		// Allocate compartments at this level
		for (k=0; k<allocTable[n].size; k++) {

			pcomp = new Compartment(
				Number( allocTable[n].diamInMicrons/2 * micron ),	// radius
				Number( allocTable[n].lenInMicrons * micron ),		// length
				Cm,
				Ri,
				Rm);

			// Number compartments for debug and reporting.
			// For names, numbering starts with 1 not 0.
			sprintf(nameBuf,"C%04d-%04d",n+1,k+1);
			pcomp->componentName(nameBuf);

			// Set the leak reversal (and initial voltage)
			pcomp->Vleak(Vleak);
			pcomp->Vinit(Vleak);

			// Add the compartment to the neuron's collection
			add(pcomp);

			// Add this compartment to the junction in process
			// if this is not the first compartment.
			if (n>0) {

				// Allow two child compartments per branch
				if (k%2==0) {
					pjunct = new ElectricalJunction(_compartments[parentIndex++]);
				}

				// Add this compartment to the junction in process
				pjunct->add(pcomp);
			}
		}
	}

	// Set compartment distance from soma where the
	// first compartment is the root node (ie soma).
	setDistFromSoma();
}

RPNeuron2::~RPNeuron2() {}

// Create a default ODE solver for this neuron.
// This can be overridden if a different solver is needed.
ODESolver* RPNeuron2::defaultSolver()
{
	return new NeuronSolver;
}



// -----------------------------------------------
// RPNeuron3 class body
// -----------------------------------------------



// Constructors and a destructor
RPNeuron3::RPNeuron3(int nComp) : RPNeuron1(nComp)
{
	using namespace UOM;

	// Constants
	const Number		gNa = 120*mS/cm_2;
	const Number		gK = 36*mS/cm_2;

	// Working variables
	Compartment*		pcomp;
	int					k;

	// Populate compartments with ion channels
	for (k=0; k<_compartments.size(); k++) {

		pcomp = _compartments[k];

		// Populate compartment with ion channels
		pcomp->add(new Na_channel(gNa));
		pcomp->add(new K_channel(gK));
	}
}

RPNeuron3::~RPNeuron3() {}
