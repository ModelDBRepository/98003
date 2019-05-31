// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_math.h
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
// This header file is oriented towards mathematical functions.

// Only include this header once per compilation unit
#ifndef __BSNF_MATH_H_
#define __BSNF_MATH_H_

// ====================================================================
// MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================
#ifdef WIN32

// Disable warning C4786: symbol greater than 255 character,
#pragma warning( disable: 4786)

#endif
// ====================================================================
// END OF MICROSOFT SPECIFIC DECLARATIONS
// ====================================================================


// ====================================================================
// Header files included here by default
// ====================================================================

// Standard C++ Template Library headers
#include <vector>
#include <cmath>

// Incorporate all the names from the std library by reference
using namespace std;

// Required BNSF headers
#include "bnsf_base.h"


// ====================================================================
// Primary namespace for the framework
// ====================================================================


namespace BNSF {

	// ====================================================================
	// Prototype declarations to allow forward references.
	// See below for descriptions of the individual classes.
	// ====================================================================

	// Math utility classes
	class SparseMatrix;
	class SparseRowRef;
	class SparseEntry;

	class Interpolator;
		class LinearInterpolator;
		class CubicInterpolator;

	class RandomAlgorithm;
		class InverseCDF;
		class RatioOfUniforms;

	class RandomNumberGenerator;
		class UniformRandom;
			class LF_UniformRandom;
			class WH_UniformRandom;
			class MT19937_UniformRandom;
		class NonuniformRandom;
			class NormalRandom;
			class PoissonRandom;
			class BinomialRandom;

	// ====================================================================
	// Mathematical constants
	// ====================================================================

	static const double Pi =	3.14159265358979323846264338327950288;
	static const double E =		2.71828182845904523536028747135266249;

	// ================================================================
	// Global function declarations
	// ================================================================

	// Quick and dirty exponential calculation -- fast but not very 
	// accurate. Approximate the exp function with a Chebyshev polynomial 
	// for small x forcing QDExp(0.0) == 1.0. For x<0 maximum absolute 
	// error is less than 2.1e-5 and relative error is less than 7.4e-4. 
	// If qdexp is used to evaluate the sigmoidal function 1/(1+exp(x)), 
	// the absolute error is less than 1.3e-5. Note that many modern 
	// processors have the exp function implemented via firmware and that 
	// this approximation may or may not save processing time.

	inline double qdexp(double x) {

		// For large x revert to the standard exp routine.
		// The limit is chosen to minimize discontinuity.
		if (fabs(x)>15.21690426072246)
			return exp(x);

		int n = 0;
		double ex;

		// Scale x to a range of relative accuracy. To minimize,
		// discontinuity, limits are chosen as points where 
		// qdexp(x/2)^2=qdexp(x).
		if (x<0)	while (x<-0.49667923004623) {x /= 2; n++; }
		else		while (x> 0.50361989372701) {x /= 2; n++; }

		// Evaluate a Chebyshev polynomial approximation for exp(x)
		ex = (((
			  0.04210263712312  * x
			+ 0.16928638488963) * x
			+ 0.49997272146577) * x
			+ 0.99983602435444) * x
			+ 1.0;

		// Square to compensate for scaling of x value.
		while (n>0) {ex *= ex; n--; }
		return ex;
	}

	// General note: the functions min, max, fmin & fmax
	// are frequently defined elsewhere (i.e. via math.h
	// cmath, algorithm, or windef.h) such that there is no
	// easy way to use those names here without explicitly
	// qualifying by namespace and, in some cases, not even
	// then. There are some subtleties surrounding handing 
	// of NaN values that are not critical here. Hence the
	// more efficient inline versions are defined with 
	// hopefully unique function names. Such is C++.

	// Utility formula for minimum of two real numbers
	inline double minval(double x, double y) { return x<y ? x:y; }
	inline float  minval(float  x, float  y) { return x<y ? x:y; }

	// Utility formula for maximum of two real numbers
	inline double maxval(double x, double y) { return x>y ? x:y; }
	inline float  maxval(float  x, float  y) { return x>y ? x:y; }

	// Utility function for returning x when x is positive
	// and 0 otherwise, in other words, the rectified value.
	inline double plusval(double x) { return x<0.0 ? 0.0 : x; }
	inline float  plusval(float x)  { return x<0.0f ? 0.0f : x; }

	// Utility formula for Euclidean norm in 2D. Note that this is
	// a fast method but not the most accurate possible.
	inline double norm(double x, double y) {return sqrt(x*x+y*y); }

	// Utility formula for vector Euclidean norm (L-2 norm)
	double norm(int n, double* X);
	double norm(valarray<double> X);
	float  norm(valarray<float> X);

	// Utility formula for vector Euclidean distance (L-2 norm)
	double dist(int n, double* X, double* Y);
	double dist(valarray<double> X, valarray<double> Y);
	float  dist(valarray<float> X, valarray<float> Y);
	
	// Utility formula for vector dot product
	double dot(int n, double* X, double* Y);
	double dot(valarray<double> X, valarray<double>Y);
	float  dot(valarray<float> X, valarray<float>Y);

	// Utility formulas for the gamma function and log thereof
	double gamma(double x);
	double logGamma(double x);

	// Utility formula for the beta function
	double beta(double z, double w);

	// Utility formulas for factorial and natural log thereof
	// Doubles are returned to allow larger numbers.
	double factorial(int n);		// usual factorial
	double logFactorial(int n);		// log of factorial

	// Utility formulas for combinations and natural log thereof
	// Doubles are returned to allow larger numbers.
	double comb(int n, int m);		// n things taken m at a time
	double logComb(int n, int m);	// log of comb

	

	// ====================================================================
	// Math utility classes
	// ====================================================================



	// ----------------------------------------------------------------
	// CLASS:	SparseEntry
	// EXTENDS:	none
	// DESC:	Represent a single value in a Sparse matrix.
	//
	// RESP:
	//		1.	Store value of the entry
	//		2.	Store row and col of the entry
	//		3.	Store link to next entry by row and col.
	//
	// NOTE:	This class is just a structure. A stdlib::vector
	//			can be used to store them, though this is not
	//			done for performance reasons.
	// ----------------------------------------------------------------

	class SparseEntry {
	public:
		double			value;
		int				row;
		int				col;
		int				nextInRow;
		int				nextInCol;
	};

	// ----------------------------------------------------------------
	// CLASS:	SparseMatrix
	// EXTENDS:	none
	// DESC:	Represent the sparse matrix used in solving
	//			ODEs. This is a specialized (and also limited)
	//			implementation of sparse matrixes emphasizing
	//			performance over storage optimization.
	//
	// RESP:
	//		1.	Store matrix data.
	//		2.	Multiply matrix times vector.
	//		3.	Perform LU decomp (without pivoting).
	//		4.	LU back substitution solving A*x=b
	//		5.	Iterate over non-zero elements.
	//
	// NOTE:	This class should be considered provisional pending
	//			availability of a clean, easy to use, object-oriented 
	//			implementation of matrix operations in C++ as, for
	//			example, promised in TNT from NIST (still incomplete).
	//
	//			Entries other than diagonals are stored individually
	//			in SparseEntry objects that are linked together
	//			by row and by column in ascending order of row and
	//			column respectively.
	//
	//			Like all C++ arrays & matrixes, indexes start at 0.
	// ----------------------------------------------------------------

	class SparseMatrix {

	public:

		// Standard constructor.
		SparseMatrix(
			int dim1=0,				// number of rows (can be 0)
			int dim2=0,				// number of columns (can be 0)
			double diagValue=0,		// initial value of diagonal terms
			int capacity=0);		// initial allocation for off diag entries.

		// Copy constructor (performs deep copy)
		SparseMatrix(const SparseMatrix& A);

		// Destructor
		virtual ~SparseMatrix();

		// Overload the assignment operator to do a deep copy
		virtual SparseMatrix&	operator=(SparseMatrix& A);

		// Accessors for dimensions and diagonal size
		inline  int				dim1() { return _dim1; }	// rows
		inline  int				dim2() { return _dim2; }	// columns
		inline  int				diagSize() { return _dim1<_dim2 ? _dim1 : _dim2; }

		// Return the number of off-diagonal data elements currently stored.
		// Elements that have been removed can still be counted here.
		inline  int				dataSize() { return _dataSize; }

		// Accessor for the padding factor. This is a factor determining
		// extra storage allocated when additional entries are stored.
		inline  double			paddingFactor() { return _paddingFactor; }
		virtual void			paddingFactor(double padfact) { _paddingFactor = padfact; } 

		// Return or set the capacity in terms of allocated elements
		// Diagonal elements are not included in this count.
		inline  int				capacity() { return _dataCapacity; }
		virtual void			capacity(int n);

		// Change the size of the matrix. diagValue is inserted in any
		// new diagonal elements created.
		virtual void			resize(int dim1, int dim2, double diagValue=0);

		// Clear the matrix preserving storage but discarding off-diagonal entries
		virtual void			clear(double diagValue=0);

		// Utility functions for resetting the matrix to special values
		inline  void			setToIdentity() { clear(1); }
		inline  void			setToZero() { clear(0); }

		// Swap data with the matrix provided avoiding reallocating storage
		virtual void			swapWith(SparseMatrix& B);

		// Get a reference to a value in a row via the [] operator.
		// For example, A[1][2] returns a reference to row 1, col 2.
		SparseRowRef			operator[](int row);

		// Get a value or reference to a specific entry in the matrix
		// Retrieving a value for an entry that is not already present
		// returns zero, while retrieving a reference to that entry results 
		// in the creation of a new entry with the indicated row and col.
		virtual double			at(int row, int col);
		virtual double&			refAt(int row, int col);

		// Set an entry to zero. If a data entry was allocated at
		// this location it is freed.
		virtual void			zeroAt(int row, int col);

		// Provide direct access to entry data for customized operations.
		// n should intially be set to 0. After each call it is adjusted to
		// reference the next item. Items are returned without ordering. 
		// Return true if an item was returned and false otherwise.
		virtual bool			nextEntry(int &n, int& row, int& col, double& val); 

		// Provide access to entries by row. n should initially be set
		// to 0. After each call, the next row entry is returned in
		// order of increasing row. Return true if an entry was returned
		// and false if not.
		virtual bool			nextInRow(int& n, int row, int& col, double& val);

		// Provide access to entries by column. n should initially be set
		// to 0. After each call, the next column entry is returned in
		// order of increasing col. Return true if an entry was returned
		// and false if not.
		virtual bool			nextInCol(int& n, int& row, int col, double& val);

		// Return the 1-norm of the matrix, that is, maximum column
		// sum of entry absolute values.
		virtual double			norm1();

		// Return the infinity-norm of the matrix, that is, maximum
		// row sum of entry absolute values.
		virtual double			normInf();

		// Scale the current matrix in place.
		virtual void			scaleBy(double c);

		// Add the current matrix to the one provided.
		// The operation is B += A where A is the current
		// matrix and B is the result.
		virtual void			addTo(SparseMatrix& B);

		// Scale the current matrix and add it to another one.
		// The operation is B += c*A where A is the current
		// matrix, c is a scalar, and B is the output.
		virtual void			addScaledByTo(double c, SparseMatrix& B);

		// Multiply a matrix times a column vector on the right.
		// The vector is represented as an array of doubles.
		// The equation is A*x=b where b is the output.
		virtual void			rightMult(double* x, double* b);

		// Multiply a matrix time a row vector on the left.
		// The vector is represented as an array of doubles.
		// The equation is x*A=b when b is the output.
		virtual void			leftMult(double*x, double* b);

		// Perform LU decomposition placing the combined results in 
		// the current matrix. If the operation suceeds return true 
		// and if not return false. tol specifies an absolute value 
		// that can be used to trim the resulting LU decomposition
		// resulting in an incomplete form. Pivoting is not done.
		virtual bool			luDecompWithoutPivot(double abstol=0);

		// Perform back substitution for an L/U combination matrix.
		// The equation is A*x=b where x is the output and b the input.
		virtual void			luBackSubRight(double* x, double* b);

		// Perform back substitution for an L/U combination matrix.
		// The equation is x*A=b where x is the output and b the input.
		virtual void			luBackSubLeft(double* x, double* b);

		// Solve the system A*x=b using the current matrix as A
		// and b as input placing the result in x. A is not changed.
		// Return true if the operation succeeds and false if not.
		// LU decomposition is used to solve the system.
		virtual bool			luSolve(double* x, double* b);

		// Solve the system A*x=b using the current matrix as A
		// and b as input placing the result in x. A is not changed.
		// Return true if the operation succeeds and false if not.
		// Preconditioned bi-conjugate gradient method is used to 
		// solve the system. If usePinv is false, then P should be
		// the LU decomposition of a matrix approximating A. If
		// usePinv is true, then P approximates A^-1 instead
		// and is not in LU decomposition format. If P is an empty
		// matrix (dim1 = dim2 = 0) or not supplied then it is not used.
		virtual bool			pbcgSolve(
			double* x,					// Solution output <out>
			double* b,					// Desired value input <in or inout>
			double* pErr=NULL,			// Error (norm-inf of residual) <out-optional>
			int*	pIterUsed=NULL,		// Iterations used <out-optional>
			double	tol=1e-15,			// Absolute tolerance for convergence <in>
			int		maxIter=0,			// Maximum iterations allowed (0 means unlimited) <in>
			double* x0=NULL,			// Initial guess for solution value <in-optional>
			SparseMatrix& P=_dummyMatrix,	// Precondition matrix(s) <in>
			bool	usePinv=false);		// P is approx A^-1 instead of P approx A <in>
		
	protected:
		int					_dim1;			// row dimension
		int					_dim2;			// column dimension
		int					_dataSize;		// entries in _data array
		int					_dataCapacity;	// allocated size of _data
		double				_paddingFactor;	// factor of extra storage allocated
		
		double*				_diag;			// diagonal elements data
		int*				_rowHead;		// head ptr for row list
		int*				_colHead;		// head ptr for col list
		int					_freeHead;		// free entry list header

		SparseEntry*		_data;			// off-diagonal data

		int					_lastRef;		// Last lookup cached index

		static SparseMatrix _dummyMatrix;	// Dummy empty matrix

		// Add a new entry to the data array expanding it if necessary.
		// Return the index of the entry added.
		virtual int			addDataEntry();

		// Copy another matrix (B) into the current one
		virtual void		copyFrom(const SparseMatrix& B);
	};

	// ----------------------------------------------------------------
	// CLASS:	SparseRowReference
	// EXTENDS:	none
	// DESC:	Allow use of A[i][j] notation is accessing
	//			references to entries in matrix A. Because these
	//			are references, an entry is always found
	//			for any valid i and j (or one is created).
	//
	// RESP:
	//		1.	Store matrix and row on creation.
	//		2.	Access value by column via operator overload.
	//
	// NOTE:	These references are considered lightweight
	//			and should not be allowed to outlast the
	//			data they are referencing. 
	//
	//			Blind use of this facility can populate the matrix
	//			with 0 values. "if (A[i][j]==0)" should be avoided.
	//			C++ matrix notation is basically unfixable,
	//			but this is a bit of syntactic sugar coating.
	// ----------------------------------------------------------------

	class SparseRowRef {
	public:

		// Constructor and destructor
		SparseRowRef(SparseMatrix* matrix, int row) : _matrix(matrix),_row(row) {}
		~SparseRowRef() {}

		// Overload [] to access reference to a value
		virtual double& operator[](int col) { return _matrix->refAt(_row,col); }

	protected:
		SparseMatrix*			_matrix;
		int						_row;
	};
	
	// ----------------------------------------------------------------
	// CLASS:	Interpolator
	// EXTENDS:	none
	// DESC:	Abstract class for defining various forms
	//			of interpolation. The common protocol is
	//			oriented towards one dimensional interpolation.
	//
	//			As a matter of notation, we interpolate the
	//			relation y=f(x). Given x we attempt to estimate
	//			a value for y given x using the data provided.
	//
	// RESP:
	//		1.	Declare common protocols
	//		2.	Provide accessors for 1-D interpolation
	//		3.	Provide utility lookup functions for subclasses
	//
	// NOTE:	double is used to provide maximum precision
	//			in the interpolation. Simple arrays of
	//			doubles are used for values rather than something
	//			fancier but less standard and harder to initialize.
	// ----------------------------------------------------------------

	class Interpolator {

	public:

		// Constructors and destructor.
		Interpolator();
		virtual ~Interpolator();

		// Get/set min and max interpolated values returned.
		// Values outside this range are quietly clipped to conform.
		inline  double		yMin() { return _yMin; }
		virtual void		yMin(double y) { _yMin = y; }
		inline  double		yMax() { return _yMax; }
		virtual void		yMax(double y) { _yMax = y; }

		// Get/set min and max for range of x values allowed.
		// Values outside this range cause an error.
		inline  double		xMin() { return _xMin; }
		virtual void		xMin(double x) { _xMin = x; }
		inline  double		xMax() { return _xMax; }
		virtual void		xMax(double x) { _xMax = x; }

		// Set data values for this interpolator
		virtual void		data(int n, double* x, double* y);	// from two arrays
		virtual void		data(int n, double* xy[2]);			// from n x 2 matrix

		// Interpolate a value for y given x
		virtual double		yAtX(double x);		

	protected:
		int					_nData;		// number of data values
		double*				_xData;		// array of x values
		double*				_yData;		// array of y values
		double				_yMin;		// min value returned
		double				_yMax;		// max value returned
		double				_xMin;		// min range value allowed
		double				_xMax;		// max range value allowed
		int					_lastFound;	// index into table for last findX

		// Prepare to load data into arrays
		virtual	void		preDataLoad();

		// Sort data into ascending sequence
		virtual void		sortData();

		// Look up the value x in the data table and return an associated
		// index of the largest value smaller than x. Return -1 if no
		// such value is found. Values of xData are assumed to be sorted
		// into ascending order.
		virtual int			findX(double x);

		// Subclass responsibilities --------------------------

		// Do whatever is necessary, such as load tables, following
		// loading of data arrays.
		virtual void		postDataLoad() {}

		// Get interpolated value of y given x
		virtual double		rawYAtX(double x) = 0;
	};

	// ----------------------------------------------------------------
	// CLASS:	LinearInterpolator
	// EXTENDS:	none
	// DESC:	Interpolate a value using a linear fit between
	//			adjacent values of x.
	//
	//			While statistically dubious, extrapolation is done.
	//			For x<xdata[0] or x>xdata[n], the fit is based
	//			on x[0] and x[1] or else x[n-1] and x[n], where
	//			n is the maximum index in the data array.
	//
	// RESP:
	//		1.	Compute interpolated values
	//
	// ----------------------------------------------------------------

	class LinearInterpolator : public Interpolator {

	public:

		// Constructor and destructor
		LinearInterpolator();

		LinearInterpolator(
			int			n,		// number of values
			double*		x,		// x data array
			double*		y);		// y data array

		LinearInterpolator(
			int			n,		// number of values
			double*		xy[2]);	// data matrix

		virtual ~LinearInterpolator();

	protected:

		// Get the interpolated value
		virtual double	rawYAtX(double x);
	};

	// ----------------------------------------------------------------
	// CLASS:	CubicInterpolator
	// EXTENDS:	Interpolator
	// DESC:	Interpolate a value fitting a piecewise cubic
	//			equation to the data points and an estimate of
	//			the first derivative. The fit has a continuous
	//			first derivative.
	//
	// RESP:
	//		1.	Compute interpolated values
	//
	// NOTE:	This is similar to the familiar spline interpolation,
	//			but is not continuous in the second derivative.
	//			However, this method is less likely to go far
	//			wrong if the data is from a noisy source.
	// ----------------------------------------------------------------

	class CubicInterpolator : public Interpolator {

	public:

		// Constructor and destructor
		CubicInterpolator();

		CubicInterpolator(
			int			n,		// number of values
			double*		x,		// x data array
			double*		y);		// y data array

		CubicInterpolator(
			int			n,		// number of values
			double*		xy[2]);	// data matrix

		virtual ~CubicInterpolator();

	protected:
		double*			_yPrime;	// y prime value for each data point
		
		// Get the interpolated value
		virtual double	rawYAtX(double x);

		// Compute the yPrime values
		virtual void	postDataLoad();

		// Estimate a derivative based on values at three points
		// Value returned is quadratic estimate of dy/dx at x1.
		// The values of x do not necessarily need to be ascending,
		// but accuracy should be better when x1 is the intermediate point.
		virtual double	yDotFit(
			double x0, double y0,
			double x1, double y1,
			double x2, double y2);
	};

	// ----------------------------------------------------------------
	// CLASS:	RandomAlgorithm
	// EXTENDS:	none
	// DESC:	Abstract class for a random number generation
	//			algorithm using a base random number generator as 
	//			the source of a density and other values.
	// RESP:
	//		1.	Save base random object.
	//		2.	Pass pass along next and inext calls (via subclass)
	// ----------------------------------------------------------------

	class RandomAlgorithm {
	
	public:
		// Constructor and destructor
		RandomAlgorithm(NonuniformRandom* baseRandom);
		virtual ~RandomAlgorithm();

		// Accessors
		inline  NonuniformRandom*	base() { return _base; }

		// Return a random number generated with the algorithm
		virtual double			next();		// random value as a double
		virtual int				inext();	// random value as an integer

	protected:
		NonuniformRandom*		_base;		// Base random object
	};

	// ----------------------------------------------------------------
	// CLASS:	InverseCDF
	// EXTENDS:	RandomAlgorithm
	// DESC:	Generate random number for a discrete distribution.
	// RESP:
	//		1.	Build CDF table by summing densities as needed.
	//		2.	Look up value in table.
	// ----------------------------------------------------------------

	class InverseCDF : public RandomAlgorithm {
	
	public:
		// Constructor and destructor
		InverseCDF(NonuniformRandom* base, int maxSize, int minValue=0);
		virtual ~InverseCDF();

		virtual int				inext();	// Return a random value
		virtual int				cdfSize();	// Size of the CDF table
		virtual double			cdfMax();	// Maximum value in CDF table

	protected:
		vector<double>			_cdf;		// partial CDF table
		int						_maxSize;	// maximum table size
		int						_minValue;	// minimum value returned
	};

	// ----------------------------------------------------------------
	// CLASS:	RatioOfUniforms
	// EXTENDS:	RandomAlgorithm
	// DESC:	Generate random number for a T-concave distribution
	//			using a generalized ratio of uniforms method.
	// RESP:
	//		1.	Generate a random value from a specified density
	// NOTE:
	//			For information about ratio-of-uniforms
	//			methods of generating arbitrary distributions
	//			see Leydold J (2001) ACM TOMS 27(1), 66-82.
	//			Some improvements in the published algorithm are
	//			included using corrections from Leydold.
	// ----------------------------------------------------------------

	class RatioOfUniforms : public RandomAlgorithm {
	
	public:
		//Constructors and destructor
		RatioOfUniforms(NonuniformRandom* base);
		virtual ~RatioOfUniforms();

		// Return a random value
		virtual double			next();		// real number value
		virtual int				inext();	// integer number value

	protected:
		// Basic parameter values
		double					_mu;		// mode value
		double					_um;		// max u value
		double					_vmin;		// min v value
		double					_vmax;		// max v value

		// Squeeze points
		double					_squpl;		// squeeze value of u for v>=0;
		double					_sqvpl;		// squeeze value of v for v>=0;
		double					_squmn;		// squeeze value of u for v<0
		double					_sqvmn;		// squeeze value of v for v<0;

		// Ratios used for squeeze
		double					_sqr1pl;	// ratio from (0,0), v>=0
		double					_sqr2pl;	// ratio from (ur,0), v>=0
		double					_sqr1mn;	// ratio from (0,0), v<0
		double					_sqr2mn;	// ratio from (ul,0), v<0
	};

	// ----------------------------------------------------------------
	// CLASS:	RandomNumberGenerator
	// EXTENDS:	none
	// DESC:	Abstract class for random number generation.
	// RESP:
	//		1.	Provide next random number (via subclass)
	//		2.	Provide density information for generic methods (ROU)
	//		3.	Get and set internal state
	//
	// NOTE:	Interfaces are provided for getting and setting
	//			internal generator states so that a previous state
	//			can be returned to. Internal states are, of course,
	//			the responsibility of the subclass as appropriate.
	//
	//			Some components of the common protocol are
	//			here to enable use of ROU methods. See
	//			the RatioOfUniforms class for background.
	// ----------------------------------------------------------------

	class RandomNumberGenerator {

	public:
		// Typedefs for this class
		enum MethodType {defaultMethod,inverseCDF,normalApprox,ROUMethod};

		// Constructors and destructor
		RandomNumberGenerator();
		virtual ~RandomNumberGenerator();
		
		// Subclass responsibilities ----------------------------------

		// Get next random value of the appropriate type
		virtual int				inext();			// get the next value as an int
		virtual	double			next();				// get the next value as a float

		// Accessor for getting and setting internal generator state.
		// Invoker is responsible for supplying an array for storing state values.
		// If this object does not have internal state, a state size of 0
		// is returned and the get and set functions are no-ops.
		// SO FAR THIS IS NOT GENERALLY IMPLEMENTED.
		virtual int				stateSize() { return 0; }
		virtual void			getState(double* stateValues) {}	// stateValues are output
		virtual void			setState(double* stateValues) {}	// stateValues are input

		// Functions used by generic methods
		virtual double			density(double x);	// Probability density at x
		virtual double			density(int k);		// Probability density at k
		virtual double			modeDensity();		// PD at mode value
		virtual double			mode();				// mode value for cont dist
		virtual int				imode();			// mode value for discrete dist

		// Return the density at k+dir where dir = +1 or -1 and the
		// density at k is already known and is provided in p.
		virtual double			adjDensity(int k, double p, int dir);
	};

	// ----------------------------------------------------------------
	// CLASS:	UniformRandom
	// EXTENDS:	RandomNumberGenerator
	// DESC:	Abstract uniform random number generator.
	// RESP:
	//		1.	Hold a default generator
	//		2.	Access random value from default generator
	//		3.	Seed default generator
	//		4.	Generate multiple streams via stride and offset	
	//
	// NOTES:	To get a random uniform using the default
	//			generator use UniformRandom::value().
	//			To use an alternative generator instance,
	//			a non-abstract class must be used.
	//			The default generator is MT19937.
	//
	//			Stride and offset are intended for use in
	//			a parallel implementation so that independent
	//			streams can be independently generated.
	//
	//			There is no good entropy generator here.
	//			Date and time are used for uniqueness,
	//			but in most cases it would be wise to seed
	//			the generator with a fixed value so that
	//			results can be replicated from one run to
	//			another.
	// ----------------------------------------------------------------

	class UniformRandom : public RandomNumberGenerator {

	public:
		UniformRandom();
		virtual ~UniformRandom();

		// Accessors
		inline  bool			isSeeded() { return _isSeeded; }
		inline  int				stride() { return _stride; }
		inline  int				offset() { return _offset; }
		virtual void			setStride(int stride, int offset);

		// Set seed from various sources
		virtual void			setSeed(int seedCount, unsigned int* seedValues);
		virtual void			setSeed(unsigned int seedValue);
		virtual void			setSeed(double seedValue); // seedValue in (0,1)
		virtual void			setSeed(); // set from date, time, and seedCount

		// Get a uniformly distributed random value
		virtual double			next();							// in (0,1)
		virtual double			next(double min, double max);	// in (min,max)

		// Cause the random number generator to advance by n
		// numbers (ignoring stride specifications).
		virtual void			advance(int n);

		// Static functions ---------------------------------

		// Get values using the static default generator
		static double			value();
		static double			value(double min, double max);

		// Choose an integer using the default uniform generator
		static  int				choice(int n);				// choose in [0,n-1]
		static  int				choice(int min, int max);	// choose in [min,max-1]

		// Set seed for static default generator
		static void				setDefaultSeed(int seedCount, unsigned int* seedValues);
		static void				setDefaultSeed(unsigned int seedValue);
		static void				setDefaultSeed();	// set from date and time

		// Get and set default generator.
		// Setting a new generator deletes the old one
		static inline UniformRandom*	defaultGen() { return _defaultGen; }
		static void						defaultGen(UniformRandom* def);


	protected:
		bool					_isSeeded;		// has a seed been provided.
		int						_stride;		// difference between indep. streams
		int						_offset;		// offset of this instance

		static unsigned int		_seedCount;		// number of time seeding was done
		static UniformRandom*	_defaultGen;	// default instance

		// Subclass responsibilities -----------------------------

		virtual double			basicValue()=0;		// return random value

		// Functions associated with setting seeds
		virtual void			resetSeed();				// clear the current seed
		virtual void			addSeed(unsigned int s);	// add next seed value
	};

	// ----------------------------------------------------------------
	// CLASS:	WH_UniformRandom
	// EXTENDS:	UniformRandom
	// DESC:	Uniform random number generator based on
	//			a Wichmann & Hill uniform generator.
	// RESP:
	//		1.	Generate uniform random numbers.
	//
	// NOTE:	Algorithm taken from Gentle, CSI 801 notes.
	//			See Wichmann BA & Hill ID (1982) 
	//			Applied Statistics 31, 188-190.
	//			corrections 1984 ibid 33, 123.
	//
	//			This is a reasonably quick generator with
	//			good uniformity. It is formed by combining
	//			three linear congruence generators.
	// ----------------------------------------------------------------

	class WH_UniformRandom : public UniformRandom {

	public:
		WH_UniformRandom();
		virtual ~WH_UniformRandom();

	protected:
		unsigned int		_P[3];		// Prime modulus values
		unsigned int		_A[3];		// Multipliers
		unsigned int		_X[3];		// State values

		virtual double		basicValue();
		virtual void		resetSeed();
		virtual void		addSeed(unsigned int);

		virtual int			stateSize() { return 3; }
		virtual void		getState(double* stateValues);	// stateValues are output
		virtual void		setState(double* stateValues);	// stateValues are input
	};

	// ----------------------------------------------------------------
	// CLASS:	LF_UniformRandom
	// EXTENDS:	UniformRandom
	// DESC:	Uniform random number generator based on
	//			an additive lagged-Fibonacci generator.
	// RESP:
	//		1.	Generate uniform random numbers.
	//		2.	Load starting values from an internal generator.
	//
	// NOTE:	Generated numbers use the (521,32) irreducible
	//			polynomial. This is intended to be a fast generator
	//			when a large number of values are to be obtained.
	//
	//			The default generator for filling state values is provided
	//			via the constructor. A new initialization generator can 
	//			be supplied at any time but if possible should be done before 
	//			the  first random number is generated. Any new generator is 
	//			owned by this object and is deleted during this object's 
	//			destruction or whenever a new generator is supplied.
	//
	//			To allow setting of parameters for the initialization
	//			generator, the state vector is only filled when the
	//			first value is accessed, either initially or after a
	//			change of initialization generators.
	// ----------------------------------------------------------------

	class LF_UniformRandom : public UniformRandom {

	public:

		// Constructor and destructor
		LF_UniformRandom(UniformRandom* initGen=new WH_UniformRandom);
		virtual ~LF_UniformRandom();

		// Get/set generator used to fill initial state
		inline  UniformRandom*	initGen() { return _initGen; }
		virtual void			initGen(UniformRandom* initializer);

		// Seed the initialization generator
		virtual void			setSeed(int seedCount, unsigned int* seedValues);
		virtual void			setSeed(unsigned int seedValue);
		virtual void			setSeed();

		// Force filling of the state vector with new random values
		virtual void			fillX();

	protected:
		int						_offset;		// Offset into X for current value
		double					_X[521];		// Random state values
		UniformRandom*			_initGen;		// Generator to fill X with values

		virtual double			basicValue();	// Return a random value

		// Static constants for this generator
		static const int		_L;				// i.e. 521
		static const int		_K;				// i.e. 32
	};

	// ----------------------------------------------------------------
	// CLASS:	MT19937_UniformRandom
	// EXTENDS:	UniformRandom
	// DESC:	Generates random numbers using the Mersenne Twister 
	//			algorithm by Takuji Nishimura and Makoto Matsumoto.
	// RESP:
	//		1.	Generate uniform random numbers.
	//		2.	Initialize from a seed value or array of values.
	//
	// NOTE:	This class is based on the C program by Nishimura and 
	//			Matsumoto. Primary changes here are relocation of 
	//			state variables mt and mti to be instance data in 
	//			objects of this class as opposed to being globals.
	// ----------------------------------------------------------------

	class MT19937_UniformRandom : public UniformRandom {

	public:

		// Constructor and destructor
		MT19937_UniformRandom();
		virtual ~MT19937_UniformRandom();

		// Interface functions ----------------------------------------

		// Set seed different in ways
		virtual void			setSeed() { UniformRandom::setSeed(); } // seed from time
		virtual void			setSeed(int seedCount, unsigned int* seedValues);
		virtual void			setSeed(unsigned int seedValue);

		// Functions provided in the original package -----------------
		void					init_genrand(unsigned long s);
		void					init_by_array(unsigned long init_key[], int key_length);
		unsigned long			genrand_int32(void);
		long					genrand_int31(void);
		double					genrand_real1(void);
		double					genrand_real2(void);
		double					genrand_real3(void);
		double					genrand_res53(void); 

	protected:

		// State variables. The size of the mt array must agree with 
		// the value defined for N in the constants below.
		unsigned long			mt[624];		// state vector
		int						mti;			// state index

		// Get a single random value per framework interfaces
		virtual double			basicValue() { return genrand_res53(); }

		// Static constants -------------------------------------------
		static const int			N;			// 624
		static const int			M;			// 397
		static const unsigned long	MATRIX_A;	// 0x9908b0dfUL /* constant vector a */
		static const unsigned long	UPPER_MASK;	// 0x80000000UL /* most significant w-r bits */
		static const unsigned long	LOWER_MASK;	// 0x7fffffffUL /* least significant r bits */
	};

	// ----------------------------------------------------------------
	// CLASS:	NonuniformRandom
	// EXTENDS:	RandomNumberGenerator
	// DESC:	Abstract class for random number generators based on
	//			on a uniform random number generator (most non-uniform
	//			number generators).
	// RESP:
	//		1.	Access either the default or a saved uniform random.
	//		2.	Generate a uniform random number (as a utility)
	//
	// NOTE:	The uniform generator is assumed to be shared across
	//			multiple random number generators to allow a common
	//			stream of random numbers. Hence the uniform generator
	//			is not owned by instances of this class.
	// ----------------------------------------------------------------

	class NonuniformRandom : public RandomNumberGenerator {
	
	public:

		// Constructor and destructor
		NonuniformRandom(UniformRandom* unif=NULL);
		virtual ~NonuniformRandom();

		// Accessors for uniform random number generator
		// To avoid interactions among streams of uniform random numbers,
		// the source generator cannot be changed once it is set.
		inline  UniformRandom*	uniformRandom() { return _uniformRandom; }
		virtual void			uniformRandom(UniformRandom* unifRand);

		// Utility access to uniform random numbers -------------------

		// Get a uniform random number in the range (0,1) using
		// the provided uniform generator or, if none, the default generator.
		static inline  double	runif(UniformRandom* unif)
		{ return unif==NULL ? UniformRandom::value() : unif->next(); }

		// Get a uniform random number in the range (0,1) using
		// the saved uniform generator or, if none, the default generator.
		inline	double			runif() { return runif(_uniformRandom); }

	protected:
		UniformRandom*			_uniformRandom;		// Base uniform random or NULL
	};

	// ----------------------------------------------------------------
	// CLASS:	NormalRandom
	// EXTENDS:	NonuniformRandom
	// DESC:	Generator for normal density
	// RESP:
	//		1.	Generate a random value with the required density.
	//
	// NOTE:	This uses the Box-Muller method for generating values
	//			in N(0,1) from pairs of uniform random values.
	// ----------------------------------------------------------------

	class NormalRandom : public NonuniformRandom {

	public:
		// Constructors and destructor
		NormalRandom(double mean=0,double sdev=1, UniformRandom* unif=NULL);
		virtual ~NormalRandom();

		// Accessors
		inline  double			mean() { return _mean; }
		inline  double			sdev() { return _sdev; }

		// Return a normally distributed random value
		// using the mean and stdev previously supplied
		virtual double			next();

		// Return a random value using parameters and the uniform
		// random stream provided or, if none, the default stream.
		static  double			value(double mean=0, double sdev=1, 
				UniformRandom*	unif=NULL);

	protected:
		double					_mean;	// Mean of values
		double					_sdev;	// Standard deviation of values

		int						_nX;	// Number of cached values
		double					_X[2];	// Cached pairs of values
	};

	// ----------------------------------------------------------------
	// CLASS:	ExponentialRandom
	// EXTENDS:	NonuniformRandom
	// DESC:	Generator for exponential density
	// RESP:
	//		1.	Generate a random value with the required density.
	//
	// NOTE:	This distribution is sometimes stated in terms of
	//			a parameter lambda = rate = 1/mean.
	// ----------------------------------------------------------------

	class ExponentialRandom : public NonuniformRandom {

	public:

		// Constructors and destructor
		ExponentialRandom(double mean=1, UniformRandom* unif=NULL);
		virtual ~ExponentialRandom();

		// Accessors
		inline	double				mean() { return _mean; }
		inline  double				rate() { return 1/_mean; }

		// Return an exponentially distributed random value
		virtual double				next();

		// Static public functions ------------------------
		static double				value(double mean=1, UniformRandom* unif=NULL);

	protected:
		double						_mean;		// Mean of values generated
	};

	// ----------------------------------------------------------------
	// CLASS:	PoissonRandom
	// EXTENDS:	NonuniformRandom
	// DESC:	Generator for Poisson density random numbers
	// RESP:
	//		1.	Generate a random value with the required density.
	//
	// NOTE:	Since this distribution does not simply scale,
	//			an instance is used to generate the random
	//			value. See functions next and inext.
	//
	//			This distribution is frequently stated in terms
	//			of a parameter lambda (= mean).
	//
	//			InverseCDF method uses a table of size proportional
	//			to the mean. Ratio-of-Uniforms (ROU) does not use a
	//			table but takes longer for each number generated.
	//			The default method selection is to use inverse CDF for
	//			means up to 50 and ROU over that.
	// ----------------------------------------------------------------

	class PoissonRandom : public NonuniformRandom {

	public:

		// Constructors and destructor
		PoissonRandom(double mean,
			UniformRandom*		unif=NULL,					// Source of uniform rand
			MethodType			method = defaultMethod);	// Generation method
		virtual ~PoissonRandom();

		// Accessors
		inline  MethodType		method() { return _method; } // method in use
		inline	double			mean() { return _mean; } // mean of this dist
		virtual double			density(int k);			// probability density at k

		// Return a Poisson distributed random value
		virtual int				inext();		// returned as integer

		// Return a Poisson distributed random value given the mean
		// using the uniform random stream given or, if none, the default.
		static	int				ivalue(double mean, UniformRandom* unif=NULL);

		// Interface for Ratio-of-Uniforms and CDF lookup support
		virtual double			modeDensity();	// PD at mode value
		virtual int				imode();		// mode value

	protected:	
		MethodType				_method;		// which method to use
		double					_mean;			// mean value
		double					_logMean;		// log of mean
		double					_modeDensity;	// density at mode
		InverseCDF*				_cdfLookup;		// Lookup algorithm
		RandomAlgorithm*		_generator;		// Generation algorithm
	};

	// ----------------------------------------------------------------
	// CLASS:	BinomialRandom
	// EXTENDS:	NonuniformRandom
	// DESC:	Generator for binomial density random numbers.
	// RESP:
	//		1.	Generate a random value with the required density.
	//
	// NOTE:	Since this distribution does not simply scale,
	//			an instance is used to generate the random
	//			value. See functions next and inext.
	//
	//			The default method uses a CDF method when the sample
	//			size is up to 50. Otherwise a ratio-of-uniforms (ROU)
	//			method is used. If accuracy is not critical, a fast 
	//			normal approximation can be used.
	// ----------------------------------------------------------------

	class BinomialRandom : public NonuniformRandom {

	public:

		// Constructors and destructor
		BinomialRandom(int size, double probability,
			UniformRandom*		unif=NULL,					// Source of uniform rand
			MethodType			method = defaultMethod);	// Generation method

		virtual ~BinomialRandom();

		// Accessors
		inline	double			prob() { return _prob; }	// Probability
		inline  double			size() { return _N; }		// Sample size
		inline  MethodType		method() { _method; }		// Method used

		// Return a binomial distributed random value
		virtual int				inext();		// returned as integer

		// Get a binomial value using the uniform random stream provided
		// or if none, the default uniform stream.
		static  int				ivalue(int size, double prob, UniformRandom* unif=NULL);

		// Return the probability density at a point
		virtual double			density(int k);

		// Return the density at k+dir where dir = +1 or -1 and the
		// density at k is already known and is provided in p.
		virtual double			adjDensity(int k, double p, int dir);

		// Ratio-of-Uniforms and CDF lookup support
		virtual double			modeDensity();	// PD at mode value
		virtual int				imode();		// mode value

	protected:
		MethodType				_method;		// which method to use
		double					_prob;			// probability on each trial
		double					_logP;			// log of probability
		double					_logQ;			// log of 1-probability
		int						_mode;			// mode cached
		double					_modeDensity;	// mode density cached
		int						_N;				// number of trials
		RandomAlgorithm*		_generator;		// Generation algorithm to use
	};

}; // end of namespace

#endif // #ifndef __BSNF_MATH_H_
