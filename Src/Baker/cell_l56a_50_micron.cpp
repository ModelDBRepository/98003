// Mammalian CA3 cell morphology
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: cell_l56a_50_micron.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// Returns CA3 cell morphology in which compartment size is limited 
// to less than 50 microns (approx).
//
// SWC conversion parameters are as follows:
//
//	const float		maxCompLen = 50;	// Max compartment size (a few may be larger)
//	const float		rChgRelTol = 0.25f;	// Fractional change allowed in radius along branch
//	const float		rChgTolLen = 5;		// Minimum branch length when dividing for radius change
//	const float		rMaxForDendrite = 5.0;	// Max dendrite size next to soma (larger are merged into soma)	
//
//	const bool		useYAxisForSomaArea = false;	// Whether to assume soma aligned on Y axis
//	const bool		skipAxon = true;				// Skip writing the axon
//	const bool		debugBranchRadiusChg = false;	// Notify on each change
//
//	const float		maxXJump = 30;		// Maximum jump in x (microns)
//	const float		maxYJump = 30;		// Maximum jump in y (microns)
//	const float		maxZJump = 30;		// Maximum jump in z (microns)
//
// Cell geometry is taken from the Southampton archives.
// This file was created by the utility program swc_reader.
//
// References:
//
// See the Duke/Southamptom Cell Archive at http://neuron.duke.edu/cells/
// for information about SWC files (this was previously located at
// http://www.neuro.soton.ac.uk/cells/). The README file associated
// with the program CVAPP has further information about the SWC format.

#include "bnsf.h"
using namespace BNSF;

// Declare the morphology function as part of the namespace
namespace BAKER_2003 {
	MorphologyEntry* cell_l56a_50_micron();
}

// Provide the body of the function
MorphologyEntry* BAKER_2003::cell_l56a_50_micron() 
{
	
static MorphologyEntry cellMorphology[] = {

// ORIGINAL_SOURCE Neurolucida
// CREATURE rat F344
// REGION Hippocampus
// FIELD/LAYER CA3
// TYPE CA3b Pyramidal Cell in vivo young
// CONTRIBUTOR Buzsaki_G & Turner_DA
// REFERENCE J. Comp. Neurol. 356: 580-594, 1995
// RAW l56a.asc
// EXTRAS Turner_P.CA3
// SOMA_AREA 1.05E3
// SHRINKAGE_CORRECTION 1.33 1.33 2.5
// VERSION_NUMBER 2.0
// VERSION_DATE 1998-03-27
// *********************************************
// SCALE 1.33  1.33  2.5  
// 
{ 1,   0,  -1,   -1,   4.125,  63.473,    0.000,   -0.6,  -1.3,   1.0,   0.0,  0.0,  0.0},
{ 3,   1,   0,   11,   1.441,  29.615,   14.808,   -2.2,  -8.3,  -5.9,  -3.2,-14.1,-13.8},
{ 3,   2,   1,   11,   1.441,  29.615,   44.423,   -5.4, -22.4, -19.6,  -3.2,-14.1,-13.8},
{ 3,   3,   2,   18,   1.009,  38.155,   78.308,  -14.0, -38.1, -34.9, -14.0,-17.2,-16.8},
{ 3,   4,   3,   20,   0.720,  17.561,  106.166,  -21.4, -48.0, -34.9,  -0.9, -2.6, 16.8},
{ 3,   5,   4,   42,   0.400,  39.445,  134.669,  -33.9, -50.5, -27.8, -24.0, -2.4, -2.6},
{ 3,   6,   5,   42,   0.400,  39.445,  174.114,  -57.9, -52.9, -30.4, -24.0, -2.4, -2.6},
{ 3,   7,   6,   57,   0.247,  47.424,  217.549,  -71.3, -65.9, -31.3,  -2.7,-23.7,  0.8},
{ 3,   8,   7,   57,   0.247,  47.424,  264.973,  -73.9, -89.6, -30.5,  -2.7,-23.7,  0.8},
{ 3,   9,   8,   57,   0.247,  47.424,  312.396,  -76.6,-113.3, -29.7,  -2.7,-23.7,  0.8},
{ 3,  10,   4,   59,   0.560,   6.433,  118.163,  -24.6, -50.7, -26.5,  -5.5, -2.8,  0.0},
{ 3,  11,  10,   67,   0.560,  23.955,  133.358,  -28.3, -60.1, -27.8,  -1.9,-16.1, -2.5},
{ 3,  12,  11,   69,   0.524,  16.258,  153.465,  -32.0, -74.8, -32.0,  -5.5,-13.3, -6.0},
{ 3,  13,  12,   72,   0.400,   8.275,  165.731,  -36.2, -85.2, -35.0,  -3.0, -7.4,  0.0},
{ 3,  14,  13,   75,   0.400,  14.488,  177.113,  -40.2, -94.1, -30.8,  -4.9,-10.5,  8.5},
{ 3,  15,  14,  106,   0.240,  41.870,  205.292,  -44.2,-111.3, -27.8,  -3.1,-24.0, -2.5},
{ 3,  16,  15,  106,   0.240,  41.870,  247.162,  -47.3,-135.3, -30.3,  -3.1,-24.0, -2.5},
{ 3,  17,  16,  106,   0.240,  41.870,  289.032,  -50.4,-159.3, -32.8,  -3.1,-24.0, -2.5},
{ 3,  18,  17,  106,   0.240,  41.870,  330.903,  -53.5,-183.3, -35.3,  -3.1,-24.0, -2.5},
{ 3,  19,  14,  141,   0.240,  39.007,  203.860,  -45.4,-108.6, -30.4,  -5.7,-18.7, -7.9},
{ 3,  20,  19,  141,   0.240,  39.007,  242.868,  -51.1,-127.3, -38.3,  -5.7,-18.7, -7.9},
{ 3,  21,  20,  141,   0.240,  39.007,  281.875,  -56.8,-145.9, -46.2,  -5.7,-18.7, -7.9},
{ 3,  22,  21,  141,   0.240,  39.007,  320.882,  -62.4,-164.6, -54.1,  -5.7,-18.7, -7.9},
{ 3,  23,  13,  144,   0.400,   9.883,  174.810,  -36.4, -92.3, -33.8,   2.8, -7.0,  2.5},
{ 3,  24,  23,  170,   0.245,  44.024,  201.764,  -31.4,-109.9, -38.3,   7.2,-28.2,-11.5},
{ 3,  25,  24,  170,   0.245,  44.024,  245.788,  -24.2,-138.1, -49.8,   7.2,-28.2,-11.5},
{ 3,  26,  25,  170,   0.245,  44.024,  289.812,  -17.0,-166.3, -61.3,   7.2,-28.2,-11.5},
{ 3,  27,  11,  183,   0.240,  32.862,  161.767,  -35.7, -77.5, -27.8, -12.9,-18.7,  2.5},
{ 3,  28,  27,  219,   0.398,  45.449,  200.922,  -47.0,-100.0, -29.2,  -9.7,-26.2, -5.3},
{ 3,  29,  28,  219,   0.398,  45.449,  246.371,  -56.8,-126.2, -34.5,  -9.7,-26.2, -5.3},
{ 3,  30,  29,  219,   0.398,  45.449,  291.819,  -66.5,-152.4, -39.8,  -9.7,-26.2, -5.3},
{ 3,  31,  30,  219,   0.398,  45.449,  337.268,  -76.3,-178.7, -45.1,  -9.7,-26.2, -5.3},
{ 3,  32,  10,  222,   0.523,  16.769,  129.765,  -30.8, -55.7, -24.3,  -6.8, -7.2,  4.5},
{ 3,  33,  32,  223,   0.400,   0.488,  138.393,  -34.3, -59.6, -22.0,  -0.2, -0.4,  0.0},
{ 3,  34,  33,  225,   0.400,  15.444,  146.359,  -32.3, -63.4, -15.5,   4.2, -7.2, 13.0},
{ 3,  35,  34,  253,   0.242,  41.098,  174.630,  -33.4, -81.9,  -8.3,  -6.6,-29.9,  1.4},
{ 3,  36,  35,  253,   0.242,  41.098,  215.727,  -40.0,-111.7,  -6.9,  -6.6,-29.9,  1.4},
{ 3,  37,  36,  253,   0.242,  41.098,  256.825,  -46.5,-141.6,  -5.5,  -6.6,-29.9,  1.4},
{ 3,  38,  37,  253,   0.242,  41.098,  297.923,  -53.1,-171.5,  -4.1,  -6.6,-29.9,  1.4},
{ 3,  39,  38,  253,   0.242,  41.098,  339.020,  -59.6,-201.3,  -2.7,  -6.6,-29.9,  1.4},
{ 3,  40,  33,  279,   0.400,  41.090,  159.182,  -41.8, -74.9, -19.1, -14.7,-30.2,  5.9},
{ 3,  41,  40,  279,   0.400,  41.090,  200.272,  -56.5,-105.0, -13.2, -14.7,-30.2,  5.9},
{ 3,  42,  41,  279,   0.400,  41.090,  241.361,  -71.2,-135.2,  -7.3, -14.7,-30.2,  5.9},
{ 3,  43,  42,  279,   0.400,  41.090,  282.451,  -85.9,-165.4,  -1.5, -14.7,-30.2,  5.9},
{ 3,  44,   3,  299,   0.240,  47.859,  121.315,  -22.0, -61.7, -32.1,  -1.9,-30.0, 22.4},
{ 3,  45,  44,  299,   0.240,  47.859,  169.174,  -23.9, -91.7,  -9.7,  -1.9,-30.0, 22.4},
{ 3,  46,  45,  316,   0.240,  49.977,  218.092,  -15.7,-128.1,   0.4,  18.1,-42.8, -2.2},
{ 3,  47,  46,  316,   0.240,  49.977,  268.069,    2.4,-170.9,  -1.8,  18.1,-42.8, -2.2},
{ 3,  48,  47,  316,   0.240,  49.977,  318.046,   20.5,-213.7,  -4.0,  18.1,-42.8, -2.2},
{ 3,  49,  45,  323,   0.400,  28.519,  207.363,  -26.0,-118.1,   4.6,  -2.3,-22.8,  6.3},
{ 3,  50,  49,  342,   0.247,  48.954,  246.099,  -25.6,-148.9,   4.7,   3.1,-38.7, -5.9},
{ 3,  51,  50,  342,   0.247,  48.954,  295.053,  -22.5,-187.6,  -1.2,   3.1,-38.7, -5.9},
{ 3,  52,  51,  342,   0.247,  48.954,  344.006,  -19.4,-226.4,  -7.1,   3.1,-38.7, -5.9},
{ 3,  53,   2,  349,   0.784,  30.145,   74.303,   -2.4, -32.2, -35.1,   9.2, -5.5,-17.1},
{ 3,  54,  53,  349,   0.784,  30.145,  104.448,    6.8, -37.8, -52.2,   9.2, -5.5,-17.1},
{ 3,  55,  54,  355,   0.565,  34.182,  136.612,   13.0, -48.1, -64.1,   3.2,-15.0, -6.8},
{ 3,  56,  55,  361,   0.717,  36.995,  172.200,   15.6, -60.9, -62.0,   2.1,-10.7, 11.0},
{ 3,  57,  56,  363,   0.560,  23.083,  202.239,   13.1, -69.4, -66.8,  -7.1, -6.3,-20.5},
{ 3,  58,  57,  381,   0.560,  49.058,  238.309,    7.1, -85.8, -78.8,  -4.9,-26.4, -3.5},
{ 3,  59,  58,  381,   0.560,  49.058,  287.368,    2.2,-112.2, -82.3,  -4.9,-26.4, -3.5},
{ 3,  60,  57,  410,   0.400,  49.943,  238.752,   20.3, -81.6, -87.5,  21.4,-18.0,-20.9},
{ 3,  61,  60,  410,   0.400,  49.943,  288.694,   41.7, -99.6,-108.4,  21.4,-18.0,-20.9},
{ 3,  62,  61,  410,   0.400,  49.943,  338.637,   63.1,-117.6,-129.3,  21.4,-18.0,-20.9},
{ 3,  63,  56,  418,   0.560,  34.218,  207.806,   21.9, -71.2, -62.0,  10.5, -9.8,-11.0},
{ 3,  64,  63,  450,   0.400,  46.009,  247.920,   38.0, -84.7, -74.3,  21.7,-17.2,-13.7},
{ 3,  65,  64,  450,   0.400,  46.009,  293.930,   59.7,-101.9, -88.0,  21.7,-17.2,-13.7},
{ 3,  66,  65,  450,   0.400,  46.009,  339.939,   81.4,-119.1,-101.7,  21.7,-17.2,-13.7},
{ 3,  67,  66,  450,   0.400,  46.009,  385.948,  103.1,-136.2,-115.4,  21.7,-17.2,-13.7},
{ 3,  68,  67,  450,   0.400,  46.009,  431.958,  124.9,-153.4,-129.1,  21.7,-17.2,-13.7},
{ 3,  69,  63,  482,   0.400,  42.920,  246.375,   33.7, -85.5, -74.5,  13.0,-18.6,-14.1},
{ 3,  70,  69,  482,   0.400,  42.920,  289.296,   46.7,-104.1, -88.6,  13.0,-18.6,-14.1},
{ 3,  71,  70,  482,   0.400,  42.920,  332.216,   59.7,-122.7,-102.6,  13.0,-18.6,-14.1},
{ 3,  72,  71,  482,   0.400,  42.920,  375.136,   72.8,-141.3,-116.7,  13.0,-18.6,-14.1},
{ 3,  73,  72,  482,   0.400,  42.920,  418.056,   85.8,-160.0,-130.7,  13.0,-18.6,-14.1},
{ 4,  74,   0,  506,   2.150,  13.005,    6.503,  -12.2,  45.1,   8.5,  -1.5,  3.3, 12.5},
{ 4,  75,  74,  509,   2.384,  14.513,   20.262,  -15.9,  51.6,  11.9,  -5.9,  9.6, -5.8},
{ 4,  76,  75,  510,   1.510,   0.925,   27.981,  -19.0,  56.8,   9.1,  -0.2,  0.9,  0.3},
{ 4,  77,  76,  522,   0.901,  34.402,   45.645,  -26.9,  67.8,  13.8, -15.8, 21.2,  9.1},
{ 4,  78,  77,  522,   0.901,  34.402,   80.047,  -42.7,  89.0,  22.9, -15.8, 21.2,  9.1},
{ 4,  79,  78,  525,   0.880,  18.775,  106.636,  -54.5, 107.6,  29.4,  -7.8, 15.9,  3.8},
{ 4,  80,  79,  535,   0.400,  29.839,  130.942,  -64.3, 126.4,  32.1, -11.8, 21.8,  1.6},
{ 4,  81,  80,  535,   0.400,  29.839,  160.781,  -76.1, 148.2,  33.7, -11.8, 21.8,  1.6},
{ 4,  82,  79,  538,   0.400,   7.619,  119.832,  -59.9, 115.6,  34.4,  -3.0,  0.2,  6.3},
{ 4,  83,  82,  561,   0.242,  42.019,  144.652,  -68.5, 118.1,  36.1, -14.2,  4.7, -2.8},
{ 4,  84,  83,  561,   0.242,  42.019,  186.671,  -82.7, 122.8,  33.4, -14.2,  4.7, -2.8},
{ 4,  85,  84,  561,   0.242,  42.019,  228.690,  -96.9, 127.5,  30.6, -14.2,  4.7, -2.8},
{ 4,  86,  78,  564,   0.560,  11.152,  102.824,  -48.4, 100.8,  31.6,   4.4,  2.4,  8.3},
{ 4,  87,  86,  566,   0.696,  14.317,  115.559,  -44.6, 104.2,  38.8,   3.2,  4.4,  6.0},
{ 4,  88,  87,  569,   0.720,  15.796,  130.616,  -41.9, 107.2,  40.8,   2.1,  1.6, -2.0},
{ 4,  89,  88,  576,   0.593,  41.068,  159.048,  -40.9, 116.1,  47.1,  -0.0, 16.2, 14.8},
{ 4,  90,  89,  584,   0.560,  49.978,  204.571,  -44.3, 132.5,  62.0,  -6.8, 16.5, 15.0},
{ 4,  91,  90,  603,   0.403,  44.853,  251.986,  -47.3, 146.9,  75.6,   0.7, 12.3, 12.2},
{ 4,  92,  91,  603,   0.403,  44.853,  296.839,  -46.6, 159.1,  87.8,   0.7, 12.3, 12.2},
{ 4,  93,  89,  609,   0.560,  31.010,  195.087,  -42.6, 138.0,  51.3,  -3.4, 27.4, -6.5},
{ 4,  94,  93,  619,   0.400,  37.333,  229.258,  -45.7, 162.7,  56.8,  -2.6, 22.1, 17.7},
{ 4,  95,  94,  642,   0.244,  49.996,  272.922,  -49.7, 193.3,  74.3,  -5.4, 39.2, 17.3},
{ 4,  96,  95,  642,   0.244,  49.996,  322.918,  -55.2, 232.6,  91.6,  -5.4, 39.2, 17.3},
{ 4,  97,  93,  674,   0.400,  44.446,  232.815,  -38.5, 163.2,  55.9,  11.8, 23.0, 15.8},
{ 4,  98,  97,  674,   0.400,  44.446,  277.261,  -26.6, 186.2,  71.6,  11.8, 23.0, 15.8},
{ 4,  99,  98,  674,   0.400,  44.446,  321.707,  -14.8, 209.2,  87.4,  11.8, 23.0, 15.8},
{ 4, 100,  88,  689,   0.240,  45.713,  161.370,  -46.2, 115.2,  47.7, -10.6, 14.3, 15.9},
{ 4, 101, 100,  689,   0.240,  45.713,  207.083,  -56.8, 129.5,  63.6, -10.6, 14.3, 15.9},
{ 4, 102, 101,  721,   0.400,  42.161,  251.020,  -71.2, 142.8,  77.7, -18.2, 12.4, 12.3},
{ 4, 103, 102,  721,   0.400,  42.161,  293.181,  -89.4, 155.2,  90.0, -18.2, 12.4, 12.3},
{ 4, 104, 103,  721,   0.400,  42.161,  335.342, -107.7, 167.6, 102.3, -18.2, 12.4, 12.3},
{ 4, 105, 101,  723,   0.318,   6.714,  233.296,  -61.5, 138.2,  68.8,   1.1,  3.2, -5.5},
{ 4, 106, 105,  759,   0.400,  45.814,  259.560,  -62.5, 150.9,  71.4,  -3.0, 22.2, 10.9},
{ 4, 107, 106,  759,   0.400,  45.814,  305.375,  -65.5, 173.0,  82.3,  -3.0, 22.2, 10.9},
{ 4, 108, 107,  759,   0.400,  45.814,  351.189,  -68.5, 195.2,  93.2,  -3.0, 22.2, 10.9},
{ 4, 109, 108,  759,   0.400,  45.814,  397.003,  -71.5, 217.4, 104.0,  -3.0, 22.2, 10.9},
{ 4, 110,  87,  762,   0.481,  10.520,  127.978,  -44.8, 109.4,  41.4,  -3.6,  5.9, -0.8},
{ 4, 111, 110,  764,   0.560,  12.372,  139.424,  -48.5, 115.2,  42.6,  -3.7,  5.9,  3.3},
{ 4, 112, 111,  770,   0.560,  27.247,  159.233,  -54.7, 129.7,  43.6,  -8.6, 23.0, -1.3},
{ 4, 113, 112,  815,   0.400,  43.491,  194.602,  -67.6, 152.9,  36.4, -17.2, 23.2,-13.1},
{ 4, 114, 113,  815,   0.400,  43.491,  238.093,  -84.8, 176.1,  23.3, -17.2, 23.2,-13.1},
{ 4, 115, 114,  815,   0.400,  43.491,  281.584, -102.0, 199.3,  10.1, -17.2, 23.2,-13.1},
{ 4, 116, 115,  815,   0.400,  43.491,  325.074, -119.2, 222.6,  -3.0, -17.2, 23.2,-13.1},
{ 4, 117, 116,  815,   0.400,  43.491,  368.565, -136.4, 245.8, -16.2, -17.2, 23.2,-13.1},
{ 4, 118, 112,  826,   0.560,  43.812,  194.763,  -66.0, 157.6,  40.3, -14.1, 32.8, -5.5},
{ 4, 119, 118,  826,   0.560,  43.812,  238.575,  -80.1, 190.4,  34.8, -14.1, 32.8, -5.5},
{ 4, 120, 119,  829,   0.560,  14.024,  267.493,  -87.3, 213.8,  32.0,  -0.2, 14.0,  0.0},
{ 4, 121, 120,  857,   0.240,  45.383,  297.197,  -95.4, 231.7,  31.3, -16.0, 21.9, -1.3},
{ 4, 122, 121,  857,   0.240,  45.383,  342.579, -111.4, 253.6,  30.0, -16.0, 21.9, -1.3},
{ 4, 123, 122,  857,   0.240,  45.383,  387.962, -127.4, 275.5,  28.7, -16.0, 21.9, -1.3},
{ 4, 124, 123,  857,   0.240,  45.383,  433.344, -143.4, 297.4,  27.4, -16.0, 21.9, -1.3},
{ 4, 125, 120,  860,   0.560,  18.766,  283.889,  -85.8, 229.8,  32.0,   3.1, 17.9,  0.0},
{ 4, 126, 125,  863,   0.684,  20.315,  303.429,  -83.2, 248.7,  32.0,   2.2, 20.0,  0.0},
{ 4, 127, 126,  908,   0.400,  47.236,  337.204,  -87.5, 273.3,  31.5, -10.9, 29.2, -1.0},
{ 4, 128, 127,  908,   0.400,  47.236,  384.441,  -98.4, 302.5,  30.5, -10.9, 29.2, -1.0},
{ 4, 129, 128,  908,   0.400,  47.236,  431.677, -109.3, 331.7,  29.5, -10.9, 29.2, -1.0},
{ 4, 130, 129,  908,   0.400,  47.236,  478.913, -120.2, 360.9,  28.5, -10.9, 29.2, -1.0},
{ 4, 131, 130,  908,   0.400,  47.236,  526.149, -131.1, 390.1,  27.5, -10.9, 29.2, -1.0},
{ 4, 132, 131,  924,   0.240,  37.701,  568.617, -133.1, 418.2,  25.1,   7.0, 26.9, -3.8},
{ 4, 133, 131,  939,   0.240,  45.875,  572.704, -143.9, 424.4,  28.6, -14.6, 39.4,  3.3},
{ 4, 134, 133,  939,   0.240,  45.875,  618.579, -158.5, 463.8,  31.9, -14.6, 39.4,  3.3},
{ 4, 135, 126,  967,   0.565,  42.189,  334.681,  -77.9, 276.1,  33.2,   8.3, 34.7,  2.4},
{ 4, 136, 135,  967,   0.565,  42.189,  376.870,  -69.6, 310.8,  35.6,   8.3, 34.7,  2.4},
{ 4, 137, 136,  967,   0.565,  42.189,  419.059,  -61.3, 345.5,  38.0,   8.3, 34.7,  2.4},
{ 4, 138, 137, 1005,   0.400,  50.709,  465.508,  -53.0, 379.1,  39.0,   8.4, 32.4, -0.5},
{ 4, 139, 138, 1005,   0.400,  50.709,  516.216,  -44.6, 411.5,  38.5,   8.4, 32.4, -0.5},
{ 4, 140, 139, 1005,   0.400,  50.709,  566.925,  -36.2, 443.9,  38.0,   8.4, 32.4, -0.5},
{ 4, 141, 140, 1005,   0.400,  50.709,  617.634,  -27.8, 476.3,  37.5,   8.4, 32.4, -0.5},
{ 4, 142, 137, 1046,   0.560,  42.364,  461.335,  -57.8, 372.5,  36.8,  -1.3, 19.2, -4.8},
{ 4, 143, 142, 1046,   0.560,  42.364,  503.699,  -59.1, 391.7,  32.0,  -1.3, 19.2, -4.8},
{ 4, 144, 143, 1046,   0.560,  42.364,  546.064,  -60.5, 410.9,  27.2,  -1.3, 19.2, -4.8},
{ 4, 145, 144, 1046,   0.560,  42.364,  588.428,  -61.8, 430.2,  22.4,  -1.3, 19.2, -4.8},
{ 4, 146, 119, 1056,   0.400,  45.324,  283.143,  -95.1, 210.1,  36.4, -15.9,  6.7,  8.8},
{ 4, 147, 146, 1056,   0.400,  45.324,  328.466, -111.0, 216.8,  45.1, -15.9,  6.7,  8.8},
{ 4, 148, 147, 1066,   0.240,  38.289,  370.273, -132.3, 232.0,  45.3, -26.8, 23.6, -8.5},
{ 4, 149, 147, 1068,   0.323,   6.837,  354.547, -122.2, 220.3,  49.5,  -6.7,  0.3,  0.0},
{ 4, 150, 149, 1074,   0.240,  35.973,  375.952, -141.5, 225.2,  46.6, -31.9,  9.4, -5.8},
{ 4, 151, 111, 1081,   0.560,  32.041,  161.631,  -57.4, 122.5,  36.4, -14.2,  8.7,-15.8},
{ 4, 152, 151, 1116,   0.240,  42.029,  198.666,  -75.5, 130.1,  28.0, -22.0,  6.5, -1.1},
{ 4, 153, 152, 1116,   0.240,  42.029,  240.695,  -97.5, 136.6,  26.9, -22.0,  6.5, -1.1},
{ 4, 154, 153, 1116,   0.240,  42.029,  282.724, -119.4, 143.1,  25.8, -22.0,  6.5, -1.1},
{ 4, 155, 154, 1116,   0.240,  42.029,  324.753, -141.4, 149.6,  24.6, -22.0,  6.5, -1.1},
{ 4, 156, 155, 1116,   0.240,  42.029,  366.782, -163.4, 156.1,  23.5, -22.0,  6.5, -1.1},
{ 4, 157, 151, 1118,   0.483,  18.157,  186.730,  -64.2, 128.0,  36.8,   0.6,  2.3, 16.5},
{ 4, 158, 157, 1141,   0.248,  39.766,  215.691,  -71.6, 133.4,  40.2, -15.2,  8.4, -9.6},
{ 4, 159, 158, 1141,   0.248,  39.766,  255.457,  -86.8, 141.8,  30.6, -15.2,  8.4, -9.6},
{ 4, 160,  76, 1158,   0.716,  48.266,   52.576,  -19.4,  75.7,  12.5,  -0.6, 36.9,  6.5},
{ 4, 161, 160, 1158,   0.716,  48.266,  100.842,  -20.0, 112.5,  19.0,  -0.6, 36.9,  6.5},
{ 4, 162, 161, 1174,   0.400,  45.998,  147.973,  -19.1, 143.3,  25.8,   2.4, 24.7,  7.0},
{ 4, 163, 162, 1174,   0.400,  45.998,  193.971,  -16.7, 168.0,  32.8,   2.4, 24.7,  7.0},
{ 4, 164, 163, 1183,   0.248,  31.339,  232.639,   -9.0, 188.4,  31.1,  12.9, 16.1,-10.3},
{ 4, 165, 164, 1198,   0.240,  31.188,  263.903,   -5.1, 201.5,  26.9,  -5.1, 10.1,  1.8},
{ 4, 166, 165, 1198,   0.240,  31.188,  295.091,  -10.2, 211.6,  28.6,  -5.1, 10.1,  1.8},
{ 4, 167, 166, 1201,   0.240,   7.567,  314.468,  -13.7, 216.9,  25.9,  -1.9,  0.4, -7.3},
{ 4, 168, 166, 1235,   0.240,  44.934,  333.152,  -11.4, 232.2,  27.8,   2.5, 30.9, -3.5},
{ 4, 169, 168, 1235,   0.240,  44.934,  378.086,   -8.9, 263.1,  24.3,   2.5, 30.9, -3.5},
{ 4, 170, 169, 1235,   0.240,  44.934,  423.019,   -6.4, 294.0,  20.8,   2.5, 30.9, -3.5},
{ 4, 171, 170, 1235,   0.240,  44.934,  467.953,   -3.8, 325.0,  17.3,   2.5, 30.9, -3.5},
{ 4, 172, 164, 1243,   0.240,  40.310,  268.464,    2.8, 213.6,  28.8,  10.6, 34.4,  5.5},
{ 4, 173, 172, 1262,   0.240,  44.858,  311.049,    8.2, 250.3,  29.2,   0.2, 39.0, -4.6},
{ 4, 174, 173, 1262,   0.240,  44.858,  355.907,    8.4, 289.4,  24.6,   0.2, 39.0, -4.6},
{ 4, 175, 172, 1272,   0.240,  40.064,  308.651,    6.6, 240.2,  32.9,  -2.9, 18.7,  2.8},
{ 4, 176, 175, 1272,   0.240,  40.064,  348.716,    3.7, 258.9,  35.6,  -2.9, 18.7,  2.8},
{ 4, 177, 176, 1272,   0.240,  40.064,  388.780,    0.8, 277.5,  38.4,  -2.9, 18.7,  2.8},
{ 4, 178, 161, 1283,   0.400,  40.174,  145.061,  -17.9, 148.2,  20.0,   4.8, 34.4, -4.5},
{ 4, 179, 178, 1283,   0.400,  40.174,  185.235,  -13.2, 182.5,  15.5,   4.8, 34.4, -4.5},
{ 4, 180, 179, 1294,   0.400,  33.193,  221.918,  -13.5, 208.3,  17.8,  -5.4, 17.1,  9.1},
{ 4, 181, 180, 1294,   0.400,  33.193,  255.111,  -18.9, 225.3,  26.9,  -5.4, 17.1,  9.1},
{ 4, 182, 181, 1297,   0.533,  15.716,  279.565,  -14.0, 235.1,  31.1,  15.3,  2.5, -0.8},
{ 4, 183, 182, 1300,   0.526,  20.024,  297.435,   -5.4, 243.5,  35.0,   1.9, 14.3,  8.5},
{ 4, 184, 183, 1305,   0.560,  40.551,  327.722,   -5.0, 263.3,  33.9,  -1.0, 25.3,-10.8},
{ 4, 185, 184, 1313,   0.407,  30.403,  363.200,   -6.3, 288.3,  31.1,  -1.7, 24.6,  5.3},
{ 4, 186, 185, 1313,   0.407,  30.403,  393.603,   -7.9, 312.8,  36.4,  -1.7, 24.6,  5.3},
{ 4, 187, 186, 1327,   0.400,  46.013,  431.811,   -4.2, 334.1,  43.9,   9.1, 18.1,  9.8},
{ 4, 188, 187, 1379,   0.244,  43.432,  476.534,    2.3, 356.5,  51.9,   3.9, 26.6,  6.3},
{ 4, 189, 188, 1379,   0.244,  43.432,  519.966,    6.2, 383.0,  58.2,   3.9, 26.6,  6.3},
{ 4, 190, 189, 1379,   0.244,  43.432,  563.398,   10.1, 409.6,  64.5,   3.9, 26.6,  6.3},
{ 4, 191, 190, 1379,   0.244,  43.432,  606.830,   14.0, 436.2,  70.8,   3.9, 26.6,  6.3},
{ 4, 192, 191, 1379,   0.244,  43.432,  650.262,   18.0, 462.8,  77.1,   3.9, 26.6,  6.3},
{ 4, 193, 186, 1390,   0.400,  28.954,  423.282,   -9.6, 334.9,  45.2,  -1.7, 19.5, 12.4},
{ 4, 194, 193, 1390,   0.400,  28.954,  452.236,  -11.3, 354.4,  57.6,  -1.7, 19.5, 12.4},
{ 4, 195, 194, 1415,   0.244,  40.664,  487.045,  -17.6, 376.8,  68.9, -11.1, 25.2, 10.3},
{ 4, 196, 195, 1415,   0.244,  40.664,  527.709,  -28.7, 402.0,  79.2, -11.1, 25.2, 10.3},
{ 4, 197, 196, 1415,   0.244,  40.664,  568.373,  -39.8, 427.2,  89.5, -11.1, 25.2, 10.3},
{ 4, 198, 197, 1415,   0.244,  40.664,  609.036,  -50.8, 452.4,  99.8, -11.1, 25.2, 10.3},
{ 4, 199, 179, 1451,   0.240,  48.732,  229.688,   -9.7, 218.0,  14.1,   2.2, 36.6,  1.7},
{ 4, 200, 199, 1451,   0.240,  48.732,  278.420,   -7.5, 254.6,  15.8,   2.2, 36.6,  1.7},
{ 4, 201, 200, 1451,   0.240,  48.732,  327.152,   -5.4, 291.2,  17.5,   2.2, 36.6,  1.7},
{ 4, 202, 201, 1451,   0.240,  48.732,  375.884,   -3.2, 327.8,  19.2,   2.2, 36.6,  1.7},
{ 4, 203,  74, 1455,   1.621,  21.144,   23.577,  -11.6,  55.4,  16.8,   2.5, 17.2,  4.0},
{ 4, 204, 203, 1467,   1.195,  34.298,   51.298,   -8.8,  78.5,  14.3,   3.1, 29.2, -9.0},
{ 4, 205, 204, 1467,   1.195,  34.298,   85.596,   -5.8, 107.7,   5.3,   3.1, 29.2, -9.0},
{ 4, 206, 205, 1470,   0.400,  11.419,  108.455,   -0.5, 124.1,   4.1,   7.4,  3.7,  6.7},
{ 4, 207, 206, 1487,   0.400,  31.792,  130.061,    9.7, 127.8,   5.9,  13.0,  3.6, -3.0},
{ 4, 208, 207, 1487,   0.400,  31.792,  161.853,   22.7, 131.4,   2.9,  13.0,  3.6, -3.0},
{ 4, 209, 208, 1502,   0.400,  32.268,  193.883,   38.3, 137.6,  -4.3,  18.1,  8.8,-11.5},
{ 4, 210, 209, 1502,   0.400,  32.268,  226.151,   56.4, 146.4, -15.8,  18.1,  8.8,-11.5},
{ 4, 211, 210, 1515,   0.240,  32.484,  258.527,   79.4, 149.5, -22.5,  28.0, -2.6, -2.0},
{ 4, 212, 211, 1515,   0.240,  32.484,  291.011,  107.4, 146.8, -24.5,  28.0, -2.6, -2.0},
{ 4, 213, 210, 1519,   0.400,  27.042,  255.806,   68.6, 155.4, -12.2,   6.3,  9.4, 18.8},
{ 4, 214, 213, 1552,   0.250,  38.480,  288.567,   80.2, 169.2,  -8.5,  16.9, 18.2,-11.4},
{ 4, 215, 214, 1552,   0.250,  38.480,  327.047,   97.2, 187.4, -19.9,  16.9, 18.2,-11.4},
{ 4, 216, 215, 1552,   0.250,  38.480,  365.527,  114.1, 205.6, -31.3,  16.9, 18.2,-11.4},
{ 4, 217, 208, 1558,   0.400,  21.480,  188.490,   35.7, 139.5,   0.3,  12.9, 12.6, -2.3},
{ 4, 218, 217, 1600,   0.243,  42.002,  220.231,   55.8, 149.1,  -2.4,  27.4,  6.7, -3.2},
{ 4, 219, 218, 1600,   0.243,  42.002,  262.233,   83.2, 155.9,  -5.6,  27.4,  6.7, -3.2},
{ 4, 220, 219, 1600,   0.243,  42.002,  304.235,  110.6, 162.6,  -8.8,  27.4,  6.7, -3.2},
{ 4, 221, 220, 1600,   0.243,  42.002,  346.237,  137.9, 169.4, -12.0,  27.4,  6.7, -3.2},
{ 4, 222, 206, 1602,   0.400,   3.152,  115.740,    2.3, 127.0,   6.8,  -1.7,  2.2, -1.3},
{ 4, 223, 222, 1620,   0.240,  48.448,  141.540,   15.3, 125.6,   0.2,  27.6, -5.1,-12.0},
{ 4, 224, 223, 1620,   0.240,  48.448,  189.988,   42.9, 120.5, -11.8,  27.6, -5.1,-12.0},
{ 4, 225, 222, 1622,   0.400,   2.965,  118.799,    0.5, 127.9,   6.9,  -1.9, -0.4,  1.5},
{ 4, 226, 225, 1631,   0.240,  28.411,  134.487,    5.4, 133.4,   0.7,  11.6, 11.4,-14.0},
{ 4, 227, 226, 1631,   0.240,  28.411,  162.897,   17.0, 144.8, -13.3,  11.6, 11.4,-14.0},
{ 4, 228, 227, 1649,   0.392,  40.161,  197.183,   32.7, 158.1, -19.5,  19.7, 15.0,  1.6},
{ 4, 229, 228, 1649,   0.392,  40.161,  237.344,   52.4, 173.1, -17.9,  19.7, 15.0,  1.6},
{ 4, 230, 229, 1649,   0.392,  40.161,  277.505,   72.1, 188.2, -16.3,  19.7, 15.0,  1.6},
{ 4, 231, 225, 1652,   0.368,  12.813,  126.687,    1.6, 131.9,  10.7,   4.0,  8.5,  6.0},
{ 4, 232, 231, 1677,   0.240,  41.068,  153.628,    8.1, 146.6,   5.5,   9.0, 20.9,-16.3},
{ 4, 233, 232, 1677,   0.240,  41.068,  194.696,   17.2, 167.5, -10.8,   9.0, 20.9,-16.3},
{ 4, 234, 233, 1677,   0.240,  41.068,  235.764,   26.2, 188.4, -27.1,   9.0, 20.9,-16.3},
{ 4, 235, 234, 1680,   0.321,  17.372,  264.984,   26.8, 202.7, -31.0,  -7.8,  7.7,  8.5},
{ 4, 236, 235, 1686,   0.240,  44.827,  296.084,   22.3, 221.4, -38.0,  -1.2, 29.7,-22.5},
{ 4, 237, 205, 1693,   0.880,  29.333,  117.412,   -2.4, 136.0,   0.7,   3.6, 27.4, -0.1},
{ 4, 238, 237, 1703,   0.400,  30.023,  147.090,    7.8, 157.3,   2.7,  16.9, 15.2,  4.0},
{ 4, 239, 238, 1725,   0.240,  36.246,  180.225,   26.0, 167.3,   1.8,  19.3,  4.9, -5.8},
{ 4, 240, 239, 1725,   0.240,  36.246,  216.470,   45.3, 172.2,  -4.1,  19.3,  4.9, -5.8},
{ 4, 241, 240, 1725,   0.240,  36.246,  252.716,   64.7, 177.1,  -9.9,  19.3,  4.9, -5.8},
{ 4, 242, 238, 1748,   0.400,  36.371,  180.287,   22.3, 174.8,   4.8,  12.0, 19.8,  0.1},
{ 4, 243, 242, 1748,   0.400,  36.371,  216.658,   34.3, 194.6,   4.9,  12.0, 19.8,  0.1},
{ 4, 244, 243, 1781,   0.240,  38.708,  254.198,   46.6, 216.0,  -0.7,  12.6, 23.0,-11.3},
{ 4, 245, 244, 1781,   0.240,  38.708,  292.906,   59.2, 239.0, -11.9,  12.6, 23.0,-11.3},
{ 4, 246, 245, 1781,   0.240,  38.708,  331.615,   71.8, 262.0, -23.2,  12.6, 23.0,-11.3},
{ 4, 247, 243, 1802,   0.400,  51.887,  260.787,   43.9, 222.9,   9.6,   7.3, 36.8,  9.4},
{ 4, 248, 247, 1802,   0.400,  51.887,  312.674,   51.2, 259.6,  19.0,   7.3, 36.8,  9.4},
{ 4, 249, 237, 1808,   0.720,  29.616,  146.887,    0.3, 163.7,   0.7,   1.9, 28.1,  0.0},
{ 4, 250, 249, 1810,   0.720,   0.890,  162.140,    1.3, 178.2,   0.7,   0.0,  0.9,  0.0},
{ 4, 251, 250, 1839,   0.240,  50.759,  187.965,   -1.4, 192.7,   2.3,  -5.5, 28.1,  3.1},
{ 4, 252, 251, 1839,   0.240,  50.759,  238.724,   -6.9, 220.8,   5.4,  -5.5, 28.1,  3.1},
{ 4, 253, 252, 1839,   0.240,  50.759,  289.483,  -12.4, 248.9,   8.5,  -5.5, 28.1,  3.1},
{ 4, 254, 253, 1839,   0.240,  50.759,  340.242,  -17.8, 277.0,  11.6,  -5.5, 28.1,  3.1},
{ 4, 255, 250, 1842,   0.720,   8.163,  166.667,    2.4, 182.4,   1.1,   2.3,  7.6,  0.8},
{ 4, 256, 255, 1849,   0.240,  21.590,  181.543,    0.6, 191.1,   7.4,  -5.9,  9.8, 12.0},
{ 4, 257, 255, 1854,   0.720,  22.892,  182.195,    5.9, 195.9,  -0.8,   4.7, 19.4, -4.5},
{ 4, 258, 257, 1895,   0.240,  49.252,  218.267,    9.9, 221.1,   1.1,   3.2, 31.0,  8.3},
{ 4, 259, 258, 1895,   0.240,  49.252,  267.518,   13.1, 252.1,   9.3,   3.2, 31.0,  8.3},
{ 4, 260, 259, 1895,   0.240,  49.252,  316.770,   16.3, 283.1,  17.6,   3.2, 31.0,  8.3},
{ 4, 261, 260, 1895,   0.240,  49.252,  366.022,   19.6, 314.1,  25.8,   3.2, 31.0,  8.3},
{ 4, 262, 257, 1897,   0.720,   4.760,  196.021,    9.0, 206.8,  -4.7,   1.5,  2.4, -3.3},
{ 4, 263, 262, 1925,   0.240,  40.063,  218.432,    2.4, 220.4,  -2.5, -14.7, 24.8,  7.6},
{ 4, 264, 263, 1925,   0.240,  40.063,  258.496,  -12.3, 245.2,   5.1, -14.7, 24.8,  7.6},
{ 4, 265, 264, 1925,   0.240,  40.063,  298.559,  -26.9, 270.0,  12.8, -14.7, 24.8,  7.6},
{ 4, 266, 265, 1925,   0.240,  40.063,  338.623,  -41.6, 294.8,  20.4, -14.7, 24.8,  7.6},
{ 4, 267, 262, 1936,   0.560,  31.350,  214.076,   12.6, 219.9,  -5.6,   5.7, 23.9,  1.5},
{ 4, 268, 267, 1956,   0.405,  42.041,  250.771,   20.4, 248.9,  -7.9,  10.0, 33.9, -6.2},
{ 4, 269, 268, 1956,   0.405,  42.041,  292.812,   30.4, 282.8, -14.2,  10.0, 33.9, -6.2},
{ 4, 270, 269, 1964,   0.400,  25.792,  326.729,   42.4, 305.2, -22.3,  14.0, 10.9,-10.0},
{ 4, 271, 270, 2002,   0.242,  49.121,  364.185,   57.0, 324.5, -25.1,  15.4, 27.6,  4.4},
{ 4, 272, 271, 2002,   0.242,  49.121,  413.306,   72.4, 352.1, -20.7,  15.4, 27.6,  4.4},
{ 4, 273, 272, 2002,   0.242,  49.121,  462.428,   87.8, 379.7, -16.3,  15.4, 27.6,  4.4},
{ 4, 274, 273, 2019,   0.240,  35.062,  504.519,   98.9, 401.3, -18.5,   6.9, 15.6, -9.0},
{ 4, 275, 274, 2019,   0.240,  35.062,  539.581,  105.8, 416.9, -27.5,   6.9, 15.6, -9.0},
{ 4, 276, 273, 2048,   0.240,  38.068,  506.022,   96.7, 407.8, -15.1,   2.3, 28.6, -2.0},
{ 4, 277, 276, 2048,   0.240,  38.068,  544.090,   99.0, 436.4, -17.0,   2.3, 28.6, -2.0},
{ 4, 278, 277, 2048,   0.240,  38.068,  582.157,  101.3, 465.0, -19.0,   2.3, 28.6, -2.0},
{ 4, 279, 269, 2062,   0.400,  40.375,  334.020,   33.0, 316.3, -22.4,  -4.8, 33.0,-10.3},
{ 4, 280, 279, 2062,   0.400,  40.375,  374.395,   28.2, 349.2, -32.7,  -4.8, 33.0,-10.3},
{ 4, 281, 280, 2102,   0.240,  44.424,  416.795,   23.8, 382.3, -39.2,  -4.0, 33.1, -2.8},
{ 4, 282, 281, 2102,   0.240,  44.424,  461.218,   19.9, 415.4, -42.0,  -4.0, 33.1, -2.8},
{ 4, 283, 282, 2102,   0.240,  44.424,  505.642,   15.9, 448.5, -44.8,  -4.0, 33.1, -2.8},
{ 4, 284, 283, 2102,   0.240,  44.424,  550.065,   11.9, 481.6, -47.6,  -4.0, 33.1, -2.8},
{ 4, 285, 280, 2146,   0.240,  44.954,  417.060,   18.7, 380.9, -30.9, -14.3, 30.5, 13.9},
{ 4, 286, 285, 2146,   0.240,  44.954,  462.014,    4.3, 411.4, -17.0, -14.3, 30.5, 13.9},
{ 4, 287, 286, 2146,   0.240,  44.954,  506.968,  -10.0, 441.8,  -3.1, -14.3, 30.5, 13.9},
{ 4, 288, 287, 2146,   0.240,  44.954,  551.922,  -24.4, 472.3,  10.8, -14.3, 30.5, 13.9},
{ 4, 289, 249, 2149,   0.345,   6.761,  165.075,   -0.7, 179.2,  17.4,  -3.9,  2.9, 33.3},
{ 4, 290, 289, 2157,   0.400,  30.023,  183.467,   -0.4, 190.4,  36.8,   4.4, 19.6,  5.5},
{ 4, 291, 290, 2177,   0.249,  44.766,  220.862,    6.0, 214.0,  45.6,   8.6, 27.4, 12.2},
{ 4, 292, 291, 2177,   0.249,  44.766,  265.628,   14.6, 241.4,  57.8,   8.6, 27.4, 12.2},
{ 4, 293, 292, 2177,   0.249,  44.766,  310.394,   23.1, 268.8,  69.9,   8.6, 27.4, 12.2},
{ 3, 294,   0, 2183,   1.040,  33.088,   16.544,   -5.3,  -8.3, -14.0,  -8.5, -7.5,-30.0},
{ 3, 295, 294, 2200,   0.599,  47.187,   56.681,  -17.7, -22.8, -21.8, -16.3,-21.5, 14.4},
{ 3, 296, 295, 2200,   0.599,  47.187,  103.868,  -34.0, -44.4,  -7.4, -16.3,-21.5, 14.4},
{ 3, 297, 296, 2204,   0.432,  40.482,  147.702,  -47.5, -59.7,   1.8, -10.6, -9.1,  4.0},
{ 3, 298, 297, 2213,   0.240,  33.120,  184.504,  -55.8, -72.0,   1.4,  -6.0,-15.4, -4.8},
{ 3, 299, 298, 2213,   0.240,  33.120,  217.624,  -61.8, -87.5,  -3.4,  -6.0,-15.4, -4.8},
{ 3, 300, 299, 2245,   0.399,  47.369,  257.869,  -65.2,-116.0,  -3.0,  -0.8,-41.7,  5.5},
{ 3, 301, 300, 2245,   0.399,  47.369,  305.237,  -66.0,-157.8,   2.5,  -0.8,-41.7,  5.5},
{ 3, 302, 301, 2245,   0.399,  47.369,  352.606,  -66.8,-199.5,   8.0,  -0.8,-41.7,  5.5},
{ 3, 303, 297, 2256,   0.240,  51.237,  193.562,  -69.6, -77.2,   0.8, -33.7,-25.9, -6.0},
{ 3, 304, 303, 2276,   0.399,  42.797,  240.579,  -97.7,-103.0,  -1.8, -22.6,-25.7,  0.8},
{ 3, 305, 304, 2276,   0.399,  42.797,  283.376, -120.3,-128.8,  -1.0, -22.6,-25.7,  0.8},
{ 3, 306, 305, 2276,   0.399,  42.797,  326.172, -142.9,-154.5,  -0.2, -22.6,-25.7,  0.8},
{ 3, 307,   0, 2280,   0.720,  13.930,    6.965,   -0.5, -11.0,  -3.5,   2.0, -6.8, -9.0},
{ 3, 308, 307, 2386,   0.720,  41.082,   34.471,    9.2, -26.7, -10.8,  17.4,-24.6, -5.5},
{ 3, 309, 308, 2390,   0.720,  16.394,   63.209,   21.0, -45.3, -11.5,   6.3,-12.6,  4.0},
{ 3, 310, 309, 2393,   0.400,  14.555,   78.684,   28.2, -51.3, -14.8,   8.0,  0.7,-10.5},
{ 3, 311, 310, 2407,   0.240,  41.110,  106.517,   46.2, -53.0, -24.9,  27.9, -4.1, -9.9},
{ 3, 312, 311, 2407,   0.240,  41.110,  147.627,   74.2, -57.2, -34.8,  27.9, -4.1, -9.9},
{ 3, 313, 310, 2416,   0.560,  51.211,  111.567,   40.9, -65.2, -18.3,  17.4,-28.5,  3.5},
{ 3, 314, 313, 2423,   0.560,  46.863,  160.604,   61.8, -81.9, -19.6,  24.4, -4.8, -6.3},
{ 3, 315, 314, 2427,   0.410,  22.050,  195.061,   83.1, -88.3, -22.8,  18.2, -8.0,  0.0},
{ 3, 316, 315, 2447,   0.400,  37.703,  224.937,  102.3, -95.7, -26.2,  20.2, -6.7, -6.9},
{ 3, 317, 316, 2447,   0.400,  37.703,  262.640,  122.5,-102.4, -33.1,  20.2, -6.7, -6.9},
{ 3, 318, 315, 2449,   0.400,   9.217,  210.694,   95.8, -94.8, -24.1,   7.2, -5.0, -2.8},
{ 3, 319, 318, 2467,   0.247,  43.392,  236.999,  108.4,-101.5, -24.9,  18.0, -8.3,  1.2},
{ 3, 320, 319, 2467,   0.247,  43.392,  280.391,  126.4,-109.8, -23.8,  18.0, -8.3,  1.2},
{ 3, 321, 320, 2467,   0.247,  43.392,  323.783,  144.4,-118.2, -22.6,  18.0, -8.3,  1.2},
{ 3, 322, 313, 2472,   0.560,  32.693,  153.519,   53.9, -90.5, -11.1,   8.5,-22.0, 10.8},
{ 3, 323, 322, 2494,   0.405,  51.580,  195.656,   70.8,-120.7,  -9.0,  25.5,-38.5, -6.5},
{ 3, 324, 323, 2494,   0.405,  51.580,  247.236,   96.3,-159.2, -15.5,  25.5,-38.5, -6.5},
{ 3, 325, 324, 2494,   0.405,  51.580,  298.816,  121.8,-197.7, -22.0,  25.5,-38.5, -6.5},
{ 3, 326, 309, 2505,   0.569,  33.875,   88.344,   27.0, -64.0,  -7.8,   5.5,-24.8,  3.5},
{ 3, 327, 326, 2510,   0.400,  22.399,  116.481,   36.7, -84.0,  -7.6,  14.0,-15.2, -3.3},
{ 3, 328, 327, 2537,   0.244,  50.171,  152.766,   50.6,-106.3, -12.0,  13.8,-29.4, -5.5},
{ 3, 329, 328, 2537,   0.244,  50.171,  202.937,   64.4,-135.7, -17.5,  13.8,-29.4, -5.5},
{ 3, 330, 329, 2537,   0.244,  50.171,  253.109,   78.2,-165.1, -23.0,  13.8,-29.4, -5.5},
{ 3, 331, 330, 2537,   0.244,  50.171,  303.280,   92.0,-194.5, -28.5,  13.8,-29.4, -5.5},
{ 3, 332, 326, 2541,   0.443,   7.787,  109.175,   30.9, -79.3,  -4.8,   2.3, -5.7,  2.5},
{ 3, 333, 332, 2575,   0.400,  44.036,  135.086,   32.0,-100.2,  -1.4,  -0.1,-36.3,  4.1},
{ 3, 334, 333, 2575,   0.400,  44.036,  179.122,   32.0,-136.5,   2.7,  -0.1,-36.3,  4.1},
{ 3, 335, 334, 2575,   0.400,  44.036,  223.158,   31.9,-172.8,   6.8,  -0.1,-36.3,  4.1},
{ 3, 336, 335, 2575,   0.400,  44.036,  267.194,   31.9,-209.1,  10.9,  -0.1,-36.3,  4.1},
{ 3, 337, 308, 2578,   0.560,  12.547,   61.286,   22.8, -41.8, -11.9,   9.9, -5.7,  3.3},
{ 3, 338, 337, 2586,   0.426,  47.120,   91.119,   41.4, -51.7, -19.5,  27.1,-14.1,-18.5},
{ 3, 339, 338, 2590,   0.552,  28.107,  128.733,   57.6, -59.6, -26.4,   5.3, -1.6,  4.6},
{ 3, 340, 339, 2590,   0.552,  28.107,  156.840,   62.8, -61.2, -21.8,   5.3, -1.6,  4.6},
{ 3, 341, 340, 2593,   0.560,  22.805,  182.296,   71.9, -66.1, -24.4,  12.9, -8.0, -9.8},
{ 3, 342, 341, 2619,   0.400,  49.662,  218.530,   90.0, -77.5, -36.8,  23.4,-14.9,-15.2},
{ 3, 343, 342, 2619,   0.400,  49.662,  268.191,  113.4, -92.5, -52.0,  23.4,-14.9,-15.2},
{ 3, 344, 343, 2619,   0.400,  49.662,  317.853,  136.7,-107.4, -67.3,  23.4,-14.9,-15.2},
{ 3, 345, 344, 2619,   0.400,  49.662,  367.515,  160.1,-122.4, -82.4,  23.4,-14.9,-15.2},
{ 3, 346, 345, 2619,   0.400,  49.662,  417.177,  183.4,-137.3, -97.6,  23.4,-14.9,-15.2},
{ 3, 347, 341, 2621,   0.483,   5.876,  196.637,   77.9, -71.6, -31.6,  -0.8, -3.1, -4.8},
{ 3, 348, 347, 2636,   0.400,  39.721,  219.435,   88.1, -78.6, -44.6,  21.2,-10.9,-21.1},
{ 3, 349, 348, 2636,   0.400,  39.721,  259.156,  109.3, -89.5, -65.7,  21.2,-10.9,-21.1},
{ 3, 350, 349, 2650,   0.250,  45.737,  301.885,  135.1,-103.5, -84.6,  30.3,-17.3,-16.8},
{ 3, 351, 350, 2650,   0.250,  45.737,  347.621,  165.4,-120.8,-101.4,  30.3,-17.3,-16.8},
{ 3, 352, 340, 2663,   0.400,  44.336,  193.062,   77.8, -56.8, -37.8,  24.7, 10.4,-36.5},
{ 3, 353, 352, 2677,   0.400,  47.641,  239.050,  102.4, -47.5, -53.0,  24.5,  8.4,  6.0},
{ 3, 354, 353, 2677,   0.400,  47.641,  286.690,  126.9, -39.1, -47.0,  24.5,  8.4,  6.0},
{ 3, 355, 352, 2685,   0.400,  42.805,  236.632,   95.0, -48.6, -57.7,   9.6,  6.1, -3.4},
{ 3, 356, 355, 2685,   0.400,  42.805,  279.437,  104.6, -42.6, -61.1,   9.6,  6.1, -3.4},

// end marker
{-1,  -1,  -1,   -1,   0.000,   0.000,    0.000,    0.0,   0.0,   0.0,   0.0,  0.0,  0.0}

};	// end of table

return cellMorphology;

}	// end of function
