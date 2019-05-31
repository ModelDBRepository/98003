// K-Dr Ion Channel Dynamics (Baker 2003)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: ionchan_k_dr_baker_2003.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file contains the classes used to implement
// the K-Dr channel based on data from Hoffman et al.
// Only the activation component is modeled here.
//
// The K+ reversal potential is chosen to fit behavior. See Bekkers
// for a measurement of -85mV in neocortical cells. Martina et al.
// claimed a measurement of -96mV for CA1. It seems that Hoffman et al.
// used -80mV though the article is not specific. Klee et al. appear
// to have assumed (as opposed to measured) a value of -78mV.
//
// References:
// 
// Bekkers JM, 2000. Properties of voltage-gated potassium currents in
// nucleated patches from large laryer 5 cortical pyramidal neurons
// of the rat. J. Physiology 525.3, 593-609.
//
// Borg-Graham LJ, 1998. Interpretations of Data and Mechanisms for 
// Hippocampal Cell Models, in Cerebral Cortex vol 13. New York: Plenum Press.
// Also available via Surf-Hippo web site.
//
// Hoffman DA, Magee JC, Colbert CM, Johnston D, 1997.
// K+ channel regulation of signal propagation in dendrites
// of hippocampal cells. Nature 387(6636), 869-875.
//
// Klee R, Ficker E, Heinemann U, 1995. Comparison of voltage-dependent 
// potassium currents in rat pyramidal neurons acutely isolated from 
// hippocampal regions CA1 and CA3.
//
// Martina M, Schultz JH, Ehmke H, Monyer H, Jonas P 1998. Functional
// and molecular differences between voltage-gated K+ channels of
// fast-spiking interneourns and pyramidal neurons of rat hippocampus.


#include "ionchan_k_dr_baker_2003.h"

using namespace std;
using namespace BNSF;
using namespace UOM;
using namespace BAKER_2003;


// ==============================================
// K_DR class bodies
// ==============================================

AlphaBetaEntry* K_DR_channel::_abTable = NULL;
const Number K_DR_channel::_Vrev = -85*mV; 

K_DR_channel::K_DR_channel(Number gSpVal)
{
	gSpecific(gSpVal);
}
