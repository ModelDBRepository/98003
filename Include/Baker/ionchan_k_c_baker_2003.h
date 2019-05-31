// K-C Ion Channel Dynamics
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_c_baker_2003.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the classes used to implement
// the K-C channel. Parameters are adapted from Moczydlowski 
// and Latorre (1982). A similar model is also distributed 
// as an example with Neuron. See also Migliori et al. 1995.
//
// Starting with values from Moczydlowski and Latorre (1982),
// calcium voltage constants (k1&k2, d1&d2) were rescaled to 
// approximate voltage sensitivity as in Gong et al. that is, 
// vhalf=0mV and  k=17mV for [Ca++]in = 2 microM. The vhalf 
// value is adjusted to the form of Boltzmann equation used 
// here rather than the variant appearing in Gong et al.
//
// Because original measurements were done at 22-24 degC, a Q10
// value typical of voltage gated channels is assumed to allow 
// corrections for temp difference. Lancaster and Adams have a
// much higher Q10 for other Ca-dep K+ currents.
//
// References:
//
// Destexhe A, Contreras D, Sejnowski TJ, Steriade M. (1994).
// A model of spindle rhythmicity in isolated thamalic reticular
// nucleus. J. Neurophysiol. 72(2), 803-181.
//
// Gong L-W, Gao T-M, Huang H, Tong Z (2001). Properties of large
// conductance calcium-activated potassium channels in pyramidal
// neurons from the hippocampal CA1 region of adult rats.
// Japanese Journal of Physiology 51, 725-731.
//
// Lancaster B and Adams PR (1986). Calcium-dependent current
// generating the afterhyperpolarization of hippocampal neurons.
// J. Neurophysiol. 55, 1268-1282. 
//
// Migliore M, Cook EP, Jaffe DB, Turner DA, and Johnston D. (1995).
// Computer simulations of morphologically reconstructed CA3
// hippocampal neurons. J. Neurophysiol. 73(3), 1157-1168.
//
// Moczydlowski E and Latorre R (1983). Gating kinetics of Ca++
// activated K+ channels from rat muscle incorporated into planar
// lipid bilayers. J. Gen. Physiol. 82, 511-542.
//
// Yoshida A, Oda M, and Ikemoto Y (1991). Kinetics of the Ca++
// activatived K+ channel in rat hippocampal neurons. Japanese
// Journal of Physiology 41: 297-315.


// Only include this header once
#ifndef __IONCHAN_K_C_BAKER_2003_H_
#define __IONCHAN_K_C_BAKER_2003_H_

#include "bnsf.h"

using namespace std;
using namespace BNSF;

namespace BAKER_2003 {

	// --------------------------------------------
	// K-C Ion Channel class
	// --------------------------------------------


	class K_C_channel : public HHIonChannel {

	public:

		// Constructors and destructor
		K_C_channel(Number gScaled=0, CalciumPool* pool=NULL);
		virtual ~K_C_channel();

		// Assumed temperature sensitivity
		virtual Number			ratedTempC() { return 22; }
		virtual Number			Q10() { return 4; }

		// Alpha beta computations
		virtual Number			alpha();
		virtual Number			beta();

		// Required functions
		virtual Number			conductance();
		virtual Number			Vrev() { return _Vrev; }

		// State vector label functions
		virtual const char* componentName() {return "K_C"; }
		virtual const char** stateLabels() { 
			static const char* sl[] = { "c" }; return sl; }

	protected:
		static const Number		_Vrev;
	};
};


#endif // #ifndef __IONCHAN_K_C_BAKER_2003_H_
