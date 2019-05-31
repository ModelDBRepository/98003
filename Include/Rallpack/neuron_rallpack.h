// Neuron for Rallpack Benchmark (rallpack3)
//
// Copyright 2006 John L Baker. All rights reserved.
//
// File: neuron_rallpack.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the declarations defining a neuron
// object for the rallpack3 benchmark. Even though a neuron is
// not explicitly used in the benchmark, one is needed to contain
// the necessary compartments.

 

// Only include this header once
#ifndef __NEURON_RALLPACK_H_
#define __NEURON_RALLPACK_H_

#include "bnsf.h"
#include "ionchan_na_rallpack.h"
#include "ionchan_k_rallpack.h"

using namespace std;
using namespace BNSF;


// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace RALLPACK {

	class RPNeuron1 : public Neuron {

	public:

		// Constructors and destructor
		RPNeuron1(int nComp=1000, Number totalLen=1*UOM::mm);
		virtual ~RPNeuron1();
		
		// Accessors
		Compartment*			nearComp()	{ return _compartments.front(); }
		Compartment*			farComp()	{ return _compartments.back(); }

		// Use NeuronSolver for this neuron
		virtual ODESolver*		defaultSolver();
	};

	class RPNeuron2 : public Neuron {

	public:

		// Constructors and destructor
		RPNeuron2();
		virtual ~RPNeuron2();
		
		// Accessors
		Compartment*			nearComp()	{ return _compartments.front(); }
		Compartment*			farComp()	{ return _compartments.back(); }

		// Use NeuronSolver for this neuron
		virtual ODESolver*		defaultSolver();
	};


	class RPNeuron3 : public RPNeuron1 {

	public:

		// Constructors and destructor
		RPNeuron3(int nComp=1000);
		virtual ~RPNeuron3();
	};
};


#endif // #ifndef __NEURON_RALLPACK_H_
