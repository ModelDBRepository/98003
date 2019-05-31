// K-AHP Ion Channel Dynamics
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_ahp_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K-AHP channel definitions. This follows the formalism
// used in Migliore but adjusted for a half activation value
// of 600 nanoM [Ca++]-in (see Hirschberg et al.)
//
// Note that Lancaster & Adams suggest a temperature sensivity
// that would imply a Q10 for the current of around 16, but
// temperature effects on calcium currents may be a factor. 
// A Q10 for this channel alone is not known and thus no
// Q10 is included in this implementation pending further data.
//
// References:
//
// Hirschberg B, Maulie J, Adelman JP, Marrion NV (1999). Gating
// properties of single SK channels in hippocampal CA1 pyramidal
// neurons. Biophysical Journal 77, 1905-1913.
//
// Lancaster B & Adams PR (1986). Calcium-dependent current
// generating the afterhyperpolatization of hippocampal neurons.
// J Neurophysiology 55, 1268-1282.
//
// Migliore M, Cook EP, Jaffe DB, Turner DA, and Johnston D. (1995).
// Computer simulations of morphologically reconstructed CA3
// hippocampal neurons. J. Neurophysiol. 73(3), 1157-1168.

 
// Only include this header once
#ifndef __IONCHAN_K_AHP_BAKER_2003_H_
#define __IONCHAN_K_AHP_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace BAKER_2003 {


	// --------------------------------------------
	// K_AHP Ion Channel classes
	// --------------------------------------------

	class K_AHP_channel : public HHIonChannel {

	public:
		// Constructors and destructor
		K_AHP_channel(Number gSpVal=0, CalciumPool* pool=NULL);
		virtual ~K_AHP_channel();

		// State vector label functions
		virtual const char*	componentName() {return "K_AHP"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "q" }; return sl; }

		// Alpha/beta Hodgkin Huxley functions
		virtual Number	alpha();
		virtual Number	beta();

		// Reversal potential for K
		virtual Number	Vrev() { return _Vrev; }

		// Conductance for the whole channel
		virtual Number	conductance();

	protected:

		static const	Number		_Vrev;

	};
};


#endif // #ifndef __IONCHAN_K_AHP_BAKER_2003_H_
