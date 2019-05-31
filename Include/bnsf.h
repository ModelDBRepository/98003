// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf.h
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This header file contains the basic classes used as a
// framework for building biological network simulations.
// As is normal for frameworks, many of these classes are
// abstract and subclassed as needed for the simulation.
//
// The following predefined symbols are tested for
// conditional compilation:
//
// WIN32 - when defined, any Microsoft specific code is used
// BNSF_DOUBLE - when defined forces use of double precision



// Only include this header once per compilation unit
#ifndef __BSNF_H_
#define __BSNF_H_

// ====================================================================
// MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================
#ifdef WIN32

// Tell compiler to read this file only once
#pragma once

// Disable warning C4786: symbol greater than 255 character,
#pragma warning( disable: 4786)

#endif
// ====================================================================
// END OF MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================


// BNSF header files to include automatically
#include "bnsf_base.h"
#include "bnsf_math.h"
#include "bnsf_sim.h"
#include "bnsf_nmod.h"

#endif // #ifndef __BSNF_H_
