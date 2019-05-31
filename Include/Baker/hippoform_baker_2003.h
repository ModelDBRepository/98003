// Provide classes for simulating cells of the hippocampal formation
//
// Copyright 2007 John L Baker. All rights reserved.
// This software is provided under the terms of the Open Source MIT License.
// See http://www.opensource.org/licenses/mit-license.php.
//
// File: hippoform_baker_2003.h
//
// Release:		1.0.1
// Author:		John Baker
// Updated:		6 March 2007
//
// Description:
//
// This header file declares classes used in the phenomenological 
// simulation of afferent place cells. Only EC, DG, and CA3 are 
// included here. Principal cells for EC, DG, and CA3 are included 
// while interneurons are included for CA3 only.
//
// Interneuron firing rates and phase relationships come largely from
// Klausberger et al. Overall phase is based on the scheme in Fox et al.
// in which peak CA3 pyramidal activity occurs at 19-deg wrt theta.
// Klausberger puts peak CA1 pyramidal activity at 20-deg which is
// pretty much the same thing.
//
// References:
//
// Amaral DG, Witter MP (1989) The three-dimensional organization 
// of the hippocampal formation: a review of anatomical data. 
// Neuroscience 31: 571-591.
//
// Brazhnik ES, Muller RU, Fox SE (2003) Muscarinic blockade slows 
// and degrades the location-specific firing of hippocampal pyramidal
// cells. Journal of Neuroscience 23: 611-621.
//
// Fox SE, Wolfson S, Ranck JB Jr. (1986) Hippocampal theta rhythm
// and firing of neurons in walking and urethane anesthetized rats.
// Experimental Brain Research 62: 495-508.
//
// Jung MW , McNaughton BL (1993) Spatial selectivity of unit 
// activity in the hippocampal granular layer. 
// Hippocampus 3: 165-182.
//
// Klausberger T, Magill PJ, Márton LF, Roberts JDB, Cobden PM, 
// Buzsaki G, Smogoyi P (2003) Brain-state- and cell-type-specific
// firing of hippocampal interneurons in vivo. Nature 421: 844-848.
//
// Klausberger T, Márton LF, Baude A, Roberts JDB, Magill PJ, 
// Smogyi P (2004) Spike timing of dendrite-targeting bistratified
// cells during hippocampal network oscillations in vivo. 
// Nature Neuroscience 7: 41-47.
//
// Miles R, Tóth K, Gulyás AI, Hájos N, Freund TF (1996) 
// Differences between somatic and dendritic inhibition in the
// hippocampus. Neuron 16: 815-823.
//
// Muller RU, Kubie JL, Ranck Jr. JB (1987) Spatial firing patterns of 
// hippocampal complex-spike cells in a fixed environment. 
// Journal of Neuroscience 7: 1935-1950.
//
// Nakazawa K, Quick MC, Chitwood RA, Watanabe M, Yeckel MF, Sun LD, 
// Kata A, Carr CA, Johnston D, Wilson MA, Tonegawa S (2002) Requirement
// for hippocampal CA3 NMDA receptors in associative memory recall.
// Science 297: 211-218.
//
// O'Keefe J, Recce ML (1993) Phase relationship between hippocampal
// place units and the EEG theta rhythm. Hippocampus 3(3): 317-330.
//
// Quirk GJ, Muller RU, Kubie JL, Ranck Jr. JB (1992) The positional
// firing properties of medial entorhinal neurons: description and 
// comparison with hippocampal place cells. 
// Journal of Neuroscience 12: 1945-1963.
//
// Skaggs WE, McNaughton BL, Wilson MA, Barnes C (1996) Theta phase
// precession in hippocampal populations and the compression of 
// temporal sequences. Hippocampus 6:149-172.
//
// Stewart M, Quirk GJ, Barry M, Fox SE (1992) Firing relations
// of medial entorhinal neurons to the hippocampal theta rhythm
// in urethane anesthetized and walking rats. 
// Experimental Brain Research 90: 21-28.
//
// Wilson MA, McNaughton BL (1993) Dynamics of the hippocampal 
// ensemble code for space. Science 261(5124): 1055-1058.



// --------------------------------------------------------------------
// Only include the definitions in this header once
// --------------------------------------------------------------------

#ifndef __HIPPOFORM_BAKER_2003_H_
#define __HIPPOFORM_BAKER_2003_H_


// --------------------------------------------------------------------
// MICROSOFT SPECIFIC DECLARATIONS
// --------------------------------------------------------------------
#ifdef WIN32

// Disable warning C4786: symbol greater than 255 character,
#pragma warning( disable: 4786)

#endif
// --------------------------------------------------------------------
// END OF MICROSOFT SPECIFIC DECLARATIONS
// --------------------------------------------------------------------


#include "bnsf.h"
#include "bnsf_liaf.h"
#include "subject_baker_2003.h"
#include "placecell_baker_2003.h"

using namespace std;
using namespace BNSF;

// Declare a namespace so that different models
// can be intermixed in the same simulation

namespace BAKER_2003 {

	// ----------------------------------------------------------------
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ----------------------------------------------------------------

	class Mouse;							// defined elsewhere

	class CellLayer;						// defined elsewhere
		class PlaceCellLayer;				// defined elsewhere
			class ECLayer;
			class DGLayer;
			class CA3Layer;

		class InterneuronLayer;				// defined elsewhere
			class CA3InterneuronLayer;
				class CA3AxoAxonicLayer;
				class CA3BasketLayer;
				class CA3BistratifiedLayer;
				class CA3OLMLayer;

	class ConnectionPolicy;					// defined elsewhere
		class ECPerforantPath;
		class DGMossyFibers;
		class CA3AssocCollaterals;

		class CA3Inhibition;
			class CA3AxonicInhibition;
			class CA3SomaticInhibition;
			class CA3ProximalInhibition;
			class CA3DistalInhibition;


	// ----------------------------------------------------------------
	// CLASS:	CA3InterneuronLayer
	// EXTENDS:	InterneuronLayer
	// DESC:	Abstract class for defining a group of interneurons.
	// ----------------------------------------------------------------

	class CA3InterneuronLayer : public InterneuronLayer {

	public:
		// Constructors and destructor
		CA3InterneuronLayer(
			Mouse*				sub, 
			int					numCells, 
			int					numericId=0,
			UniformRandom*		spkunif=NULL);

		virtual ~CA3InterneuronLayer();

		// Accessor for the maze subject as a mouse object
		inline Mouse*			mouse() { return _mouse; }

		// Accessor for the baseline firing rates
		inline  Number			baselineRate() { return _baselineRate; }
		virtual void			baselineRate(Number r);

		// Accessor for the aggregate theta phase for this layer.
		// This can be overridden on a cell-by-cell basis in subclasses.
		// Changes to this value are not automatically propagated
		// after allocateCells() has been invoked.
		inline  Number			thetaPhase() { return _thetaPhase; }
		virtual void			thetaPhase(Number tph) { _thetaPhase = tph; }

		// Accessors for input layers and rate multipliers. See below
		// for a description of how the feedforward rates are applied.

		inline  PlaceCellLayer*	ECPlaceCells() { return _ECPlaceCells; }
		inline  PlaceCellLayer*	DGPlaceCells() { return _DGPlaceCells; }
		inline  PlaceCellLayer*	CA3PlaceCells() { return _CA3PlaceCells; }

		virtual void			ECPlaceCells(PlaceCellLayer* layer);
		virtual void			DGPlaceCells(PlaceCellLayer* layer);
		virtual void			CA3PlaceCells(PlaceCellLayer* layer);

		inline  Number			ECFeedforward() { return _ECFeedforward; }
		inline  Number			DGFeedforward() { return _DGFeedforward; }
		inline  Number			CA3Feedforward() { return _CA3Feedforward; }

		virtual void			ECFeedforward(Number r) { _ECFeedforward = r; }
		virtual void			DGFeedforward(Number r) { _DGFeedforward = r; }
		virtual void			CA3Feedforward(Number r) { _CA3Feedforward = r; }

		// Provide accessors to control coupling of target cell firing rate
		// with interneuron firing rate.
		inline  Number			CA3FeedbackMax() { return _CA3FeedbackMax; }
		inline  Number			CA3FeedbackFRh() { return _CA3FeedbackFRh; }
		inline  Number			CA3FeedbackFRk() { return _CA3FeedbackFRk; }
		virtual void			CA3Feedback(Number frm, Number frh, Number frk)				
								{_CA3FeedbackMax=frm; 
								_CA3FeedbackFRh=frh; 
								_CA3FeedbackFRk=frk; }

		// Accessors for values associated with CA3FeedbackRate.
		// By default CA3FeedbackOngoingRate is initialized to zero, 
		// but subclasses (or testcases) may choose other values to 
		// avoid start-up transients in rate estimation.
		inline  Number			CA3FeedbackOngoingRate() 
								{ return _CA3FeedbackOngoingRate; }	
		virtual void			CA3FeedbackOngoingRate(Number fr) 
								{_CA3FeedbackOngoingRate = fr; }

		// An exponential smoother is applied to the target cell
		// firing rate. By default the time constant is 10 sec.
		inline  SimTime			CA3FeedbackTau() { return _CA3FeedbackTau; }
		virtual void			CA3FeedbackTau(SimTime tau) {_CA3FeedbackTau=tau; }

		// Act on an update message from the (CA3 only) input layer
		// to set new overall interneuron layer firing rates. All inputs
		// are considered in setting this rate but CA3 is assumed to
		// be the last layer with new rates and thus the signal to
		// set new rates overall.
		virtual void			updateFrom(ModelComponent* mc, int reason);

		// Act on end of the time step to update feedback firing rate.
		virtual void			timeStepEnded();

	protected:
		Mouse*					_mouse;		// maze subject

		// Provide a baseline firing rate to which other
		// contributions (feedforward, feedback) can be added.
		Number					_baselineRate;
		
		// Aggregate theta phase for cells in this layer (in radians)
		Number					_thetaPhase;

		// For each layer which stimulates interneurons,
		// provide access to the layer and a feedforward rate.
		// This rate is the multiplier of the average per-cell
		// firing rate for the afferent layer which contributes
		// towards interneuron firing. The firing rate of the 
		// interneuron layer is the sum of these contributions 
		// plus the baseline rate. Gamma and theta modulation 
		// are then applied to the aggregate firing rate.

		PlaceCellLayer*			_ECPlaceCells;
		Number					_ECFeedforward;

		PlaceCellLayer*			_DGPlaceCells;
		Number					_DGFeedforward;

		PlaceCellLayer*			_CA3PlaceCells;
		Number					_CA3Feedforward;

		// CA3 feedback allows the firing rate of the target cell to
		// affect interneuron firing rates. Ongoing firing rate of the
		// target cell is used to get the contribution.
		// The increment to interneuron firing rate is the sigmoidal:
		//
		// INrate += max/(1+exp(-(fr-frh)/frk)) 
		//
		// where fr is the target cell ongoing firing rate. Ongoing rate
		// is estimated with an exponential decay kernel applied to the 
		// cell's current firing rate using decay time constant tau.

		SimTime					_CA3FeedbackTau;	// Estimator time const.
		Number					_CA3FeedbackMax;	// Max rate increment
		Number					_CA3FeedbackFRh;	// Half activation rate
		Number					_CA3FeedbackFRk;

		// Estimated target cell firing rate and the time at which
		// the rate was determined (so that we can avoid multiple
		// updates for a given time).
		Number					_CA3FeedbackOngoingRate;
		SimTime					_CA3FeedbackAsOfTime;

		// Update the CA3 feedback rate by smoothing current firing rate.
		virtual void			updateCA3FeedbackRate();

		// Allocate a new interneuron customized for this layer.
		virtual Interneuron*	newInterneuron(int k);
	};	

	// ----------------------------------------------------------------
	// CLASS:	CA3Inhibition
	// EXTENDS:	ConnectionPolicy
	// DESC:	Abstract class for an inhibitory connection policy.
	// ----------------------------------------------------------------

	class CA3Inhibition : public ConnectionPolicy {

	public:
		CA3Inhibition(CA3InterneuronLayer* source, UniformRandom* randomizer);
		virtual ~CA3Inhibition();

		// Access hardcoded numbers of synapses to connect
		inline  int				numSomaSynapses() { return _numSomaSynapses; }
		inline  int				numISSynapses() { return _numISSynapses; }

		virtual void			numSomaSynapses(int n) { _numSomaSynapses=n; }
		virtual void			numISSynapses(int n) { _numISSynapses=n; }

	protected:

		// Numbers of synapses for soma and initial segment
		int						_numSomaSynapses;
		int						_numISSynapses;

		// Connect synapses of the number specified specified
		virtual void			connectWithSoma();
		virtual void			connectWithIS();


	};

	// ----------------------------------------------------------------
	// CLASS:	ECLayer
	// EXTENDS:	PlaceCellLayer
	// DESC:	Defines a layer of place cells corresponding
	//			with the output of entorhinal cortex layer II and III.
	// RESP:
	//		1.	Initialize place cells with EC characteristics
	//
	// NOTE:	See newPlaceCell for option of coarse or fine spatial
	//			tuning. Obviously this is only a rough approximation
	//			of either Quirk et al. coarse tuning or Hafting et al.
	//			CA-like spatial tuning for grid cells. More explicit
	//			grid cell support would be useful once we know what
	//			to assume about grid cells.
	// ----------------------------------------------------------------

	class ECLayer : public PlaceCellLayer {

	public:

		// Constructors and destructor
		ECLayer(MazeSubject* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~ECLayer();

	protected:
		virtual PlaceCell*		newPlaceCell(int k);
	};

	// ----------------------------------------------------------------
	// CLASS:	ECPerforantPath
	// EXTENDS:	ConnectionPolicy
	// DESC:	Provides parameters controlling connections. Other
	//			function are inherited from the superclass.
	// ----------------------------------------------------------------

	class ECPerforantPath : public ConnectionPolicy {

	public:
		ECPerforantPath(PlaceCellLayer* source, UniformRandom* randomizer);
		virtual ~ECPerforantPath();
	};

	// ----------------------------------------------------------------
	// CLASS:	DGLayer
	// EXTENDS:	PlaceCellLayer
	// DESC:	Defines a layer of place cells corresponding
	//			with the output of the dentate gyrus.
	// RESP:
	//		1.	Initialize place cells with DG characteristics
	// ----------------------------------------------------------------

	class DGLayer : public PlaceCellLayer {

	public:

		// Constructors and destructor
		DGLayer(MazeSubject* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~DGLayer();

	protected:
		virtual PlaceCell*		newPlaceCell(int k);
	};
 
	// ----------------------------------------------------------------
	// CLASS:	DGMossyFiber
	// EXTENDS:	ConnectionPolicy
	// DESC:	Provides parameters controlling connections. Other
	//			function are inherited from the superclass.
	// ----------------------------------------------------------------

	class DGMossyFibers : public ConnectionPolicy {

	public:
		DGMossyFibers(PlaceCellLayer* source, UniformRandom* randomizer);
		virtual ~DGMossyFibers();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3Layer
	// EXTENDS:	PlaceCellLayer
	// DESC:	Defines a layer of place cells corresponding
	//			with the output of the CA3 portion of the hippocampus.
	// RESP:
	//		1.	Initialize place cells with CA3 characteristics
	// ----------------------------------------------------------------

	class CA3Layer : public PlaceCellLayer {

	public:

		// Constructors and destructor
		CA3Layer(MazeSubject* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~CA3Layer();

	protected:
		virtual PlaceCell*		newPlaceCell(int k);
	};	
 
	// ----------------------------------------------------------------
	// CLASS:	CA3AssocCollateral
	// EXTENDS:	ConnectionPolicy
	// DESC:	Provides parameters controlling connections. Other
	//			function are inherited from the superclass.
	// ----------------------------------------------------------------

	class CA3AssocCollaterals : public ConnectionPolicy {

	public:
		CA3AssocCollaterals(PlaceCellLayer* source, UniformRandom* randomizer);
		virtual ~CA3AssocCollaterals();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3AxoAxonicLayer
	// EXTENDS:	CA3InterneuronLayer
	// DESC:	Defines a group of interneurons targeting the
	//			initial segment.
	// ----------------------------------------------------------------

	class CA3AxoAxonicLayer : public CA3InterneuronLayer {

	public:
		// Constructors and destructor
		CA3AxoAxonicLayer(Mouse* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~CA3AxoAxonicLayer();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3AxonicInhibition
	// EXTENDS:	CA3Inhibition
	// DESC:	Provides parameters for inhibitory connections
	//			that target the initial segment.
	// ----------------------------------------------------------------

	class CA3AxonicInhibition : public CA3Inhibition {

	public:
		CA3AxonicInhibition(CA3InterneuronLayer* source, UniformRandom* randomizer);
		virtual ~CA3AxonicInhibition();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3BasketLayer
	// EXTENDS:	CA3InterneuronLayer
	// DESC:	Defines a group of interneurons targeting the
	//			soma and nearby proximal dendrites.
	// ----------------------------------------------------------------

	class CA3BasketLayer : public CA3InterneuronLayer {

	public:
		// Constructors and destructor
		CA3BasketLayer(Mouse* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~CA3BasketLayer();

	protected:

		// Allocate a saw-tooth interneuron
		virtual Interneuron*		newInterneuron(int	k);

	};

	// ----------------------------------------------------------------
	// CLASS:	CA3SomaticInhibition
	// EXTENDS:	CA3Inhibition
	// DESC:	Provides parameters for inhibitory connections
	//			that target the soma.
	// ----------------------------------------------------------------

	class CA3SomaticInhibition : public CA3Inhibition {

	public:
		CA3SomaticInhibition(CA3InterneuronLayer* source, UniformRandom* randomizer);
		virtual ~CA3SomaticInhibition();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3BistratifiedLayer
	// EXTENDS:	CA3InterneuronLayer
	// DESC:	Defines a group of interneurons targeting the
	//			proximal and medium dendrites.
	// ----------------------------------------------------------------

	class CA3BistratifiedLayer : public CA3InterneuronLayer {

	public:
		// Constructors and destructor
		CA3BistratifiedLayer(Mouse* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~CA3BistratifiedLayer();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3ProximalInhibition
	// EXTENDS:	CA3Inhibition
	// DESC:	Provides parameters for inhibitory connections
	//			that target proximal apical and basal dendrites.
	// ----------------------------------------------------------------

	class CA3ProximalInhibition : public CA3Inhibition {

	public:
		CA3ProximalInhibition(CA3InterneuronLayer* source, UniformRandom* randomizer);
		virtual ~CA3ProximalInhibition();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3OLMLayer
	// EXTENDS:	CA3InterneuronLayer
	// DESC:	Defines a group of interneurons targeting the
	//			distal apical dendrites.
	// ----------------------------------------------------------------

	class CA3OLMLayer : public CA3InterneuronLayer {

	public:
		// Constructors and destructor
		CA3OLMLayer(Mouse* sub, int numCells, 
				int numericId=0, UniformRandom* spkunif=NULL);
		virtual ~CA3OLMLayer();
	};

	// ----------------------------------------------------------------
	// CLASS:	CA3DistalInhibition
	// EXTENDS:	CA3Inhibition
	// DESC:	Provides parameters for inhibitory connections
	//			that target distal apical dendrites.
	// ----------------------------------------------------------------

	class CA3DistalInhibition : public CA3Inhibition {

	public:
		CA3DistalInhibition(CA3InterneuronLayer* source, UniformRandom* randomizer);
		virtual ~CA3DistalInhibition();
	};

};

#endif // #ifndef
