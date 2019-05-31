// Basic Neural Simulation Framework (BSNF)
//
// Copyright 2007 John L Baker. All rights reserved.
//
// This software is provided AS IS under the terms of the Open Source
// MIT License. See http://www.opensource.org/licenses/mit-license.php.
//
// File: bsnf_math.cpp
//
// Release:		1.0.0
// Author:		John Baker
// Updated:		14 July 2006
//
// Description:
//
// This file provides the body of BNSF classes.
// The sequence is generally the same as the definitions
// in the associated header file, but a good PDE is
// almost a necessity for any reasonable development effort.
//
// Class bodies included here are for mathematical functions.


#include "bnsf_math.h"
#include <ctime>

using namespace std;
using namespace BNSF;


// ====================================================================
// Global functions
// ====================================================================


// Compute the Euclidean norm of a vector in several forms
double BNSF::norm(valarray<double> X)
{
	const int	n = X.size();
	int			k;
	double		normsq=0;

	for (k=0;k<n;k++) {
		normsq += X[k]*X[k];
	}
	return sqrt(normsq);
}

float BNSF::norm(valarray<float> X)
{
	const int	n = X.size();
	int			k;
	double		normsq=0;

	for (k=0;k<n;k++) {
		normsq += X[k]*X[k];
	}
	return sqrt(normsq);
}

double BNSF::norm(int n, double* X)
{
	int			k;
	double		normsq=0;

	for (k=0;k<n;k++) {
		normsq += X[k]*X[k];
	}
	return sqrt(normsq);
}

// Compute the Euclidean vector distance
double BNSF::dist(valarray<double> X, valarray<double> Y)
{
	if (X.size()!=Y.size()) {
		FatalError("(dist) Input vectors are of different sizes.");
		return 0;
	}

	const int	n = X.size();
	int			k;
	double		dsq=0;

	for (k=0;k<n;k++) {
		dsq += (X[k]-Y[k])*(X[k]-Y[k]);
	}
	return sqrt(dsq);
}

float BNSF::dist(valarray<float> X, valarray<float> Y)
{
	if (X.size()!=Y.size()) {
		FatalError("(dist) Input vectors are of different sizes.");
		return 0;
	}

	const int	n = X.size();
	int			k;
	double		dsq=0;

	for (k=0;k<n;k++) {
		dsq += (X[k]-Y[k])*(X[k]-Y[k]);
	}
	return sqrt(dsq);
}

double BNSF::dist(int n, double* X, double* Y)
{
	int			k;
	double		dsq=0;

	for (k=0;k<n;k++) {
		dsq += (X[k]-Y[k])*(X[k]-Y[k]);
	}
	return sqrt(dsq);
}

// Compute the dot product of two vectors (in several forms)
double BNSF::dot(valarray<double> X, valarray<double> Y)
{
	if (X.size()!=Y.size()) {
		FatalError("(dot) Input vectors are of different sizes.");
		return 0;
	}

	const int	n = X.size();
	int			k;
	double		dp = 0;

	for (k=0;k<n;k++) {
		dp += X[k]*Y[k];
	}
	return dp;
}

float BNSF::dot(valarray<float> X, valarray<float> Y)
{
	if (X.size()!=Y.size()) {
		FatalError("(dot) Input vectors are of different sizes.");
		return 0;
	}

	const int	n = X.size();
	int			k;
	double		dp = 0;

	for (k=0;k<n;k++) {
		dp += X[k]*Y[k];
	}
	return float(dp);
}

double BNSF::dot(int n, double* X, double* Y)
{
	int			k;
	double		dp = 0;

	for (k=0;k<n;k++) {
		dp += X[k]*Y[k];
	}
	return dp;
}

// Formula for gamma function adapted from
// Press st al. Numerical Recipies in C 2nd ed.
double BNSF::gamma(double x)
{
	return exp(logGamma(x));
}

// Formula for log of gamma function adapted from
// Press et al. Numerical Recipies in C 2nd ed.
// Computing the log avoids overflow in many cases.
double BNSF::logGamma(double x)
{
	static const double cof[7] = {
		1.000000000190015,
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5 };

	const double logSqrtTwoPi = 0.9189385471184401400;

	double temp,series;

	temp=x+5.5;
	temp -= (x+0.5)*log(temp);
	series = cof[0] 
		+ cof[1]/(x+1)
		+ cof[2]/(x+2)
		+ cof[3]/(x+3)
		+ cof[4]/(x+4)
		+ cof[5]/(x+5)
		+ cof[6]/(x+6);

	return -temp+logSqrtTwoPi+log(series/x);
}

// Formula for the beta function from 
// Press et al. Numerical Recipies in C 2nd ed.
double BNSF::beta(double z, double w)
{
	return exp(logGamma(z)+logGamma(w)-logGamma(z+w));
}

// Formula for factorial computed the obvious way
// for small values and using an approximation otherwise.
double BNSF::factorial(int n) 
{
	// table of factorial values of small numbers
	static double	ftbl[31] = {1,1};
	static int		nftbl = 2;

	// Treat factorial of negative n as 0
	if (n<0) {
		return 0;
	}

	// For small values compute exactly
	if (n<=30) {
		// Explicitly multiply 2,3,...,n
		for (;nftbl<=n;nftbl++) {
			ftbl[nftbl] = nftbl * ftbl[nftbl-1];
		}
		return ftbl[n];
	}
	
	// Otherwise, use approximation formula
	return exp(logFactorial(n));
}


// Compute log of a factorial using either an explicit
// computation of factorial or Sterling's approximation.
double BNSF::logFactorial(int n)
{

	// Table of log factorial values of small numbers
	static double	lftbl[31] = {0,0};
	static int		nlftbl = 2;

	// Use the exact function for small n and build a table
	// of values to avoid doing this twice.
	if (n<=30) {
		for (;nlftbl<=n;nlftbl++) {
			lftbl[nlftbl] = log(factorial(nlftbl));
		}
		return lftbl[n];
	}

	double				x = n;
	const double		logSqrtTwoPi = 0.9189385471184401400;

	// Stirling's formula. This improvement over the standard Stirling
	// formula can be found in Feller, An Intoduction to Probability
	// vol I, 3rd edition, 1957, p. 54 footnote 17. More complex
	// and precise formulas can be found in the standard literature.
	// The accuracy here is about the same as for the gamma function.

	return logSqrtTwoPi+(x+0.5)*log(x)-x+(1-1/(30*x*x))/(12*x);
}

// Combinations computed more or less the obvious way
// for small n and using logarithms for larger n.
// A double is returned to reduce overflows.
double BNSF::comb(int n, int m)
{
	double	c;
	int		k;

	// See if we want to use logarithms or not
	if (n>30) {
		return exp(logComb(n,m));
	}

	// Use the smallest possible m
	if (m>n/2) {
		m = n-m;
	}

	// Multiply one term at a time to avoid overflow
	c = 1;
	for (k=1;k<=m;k++) {
		c *= double(n-k+1)/k;
	}

	return c;
}

// Log of combinations. Use the exact formula for
// small n and the logFactorials expansion for larger n.
double BNSF::logComb(int n, int m)
{
	if (n<=30) {
		return log(comb(n,m));
	}

	return logFactorial(n)-logFactorial(m)-logFactorial(n-m);
}



// ====================================================================
// SparseMatrix class body
// ====================================================================



// Static data values
SparseMatrix SparseMatrix::_dummyMatrix;	// allows default empty matrix in params

// Standard Constructor
SparseMatrix::SparseMatrix(int dim1, int dim2, double initValue, int initDataCap)
{
	int		diagSize = dim1<dim2 ? dim1 : dim2;

	// Set initial defaults and values
	_paddingFactor = 0.5;
	_dim1 = 0;
	_dim2 = 0;
	_freeHead = -1;
	_dataSize = 0;
	_dataCapacity = 0;
	_diag = NULL;
	_rowHead = NULL;
	_colHead = NULL;
	_data = NULL;
	_lastRef = -1;

	// Allocate storage for the data arrays
	capacity(initDataCap);

	// Let resize build the rest
	resize(dim1,dim2,initValue);
}

// Copy constructor
SparseMatrix::SparseMatrix(const SparseMatrix& A)
{
	copyFrom(A);
}

// Destructor
SparseMatrix::~SparseMatrix()
{
	// Free all explicitly allocated arrays
	delete[] _rowHead;
	delete[] _colHead;
	delete[] _diag;
	delete[] _data;
}

// Overloaded assignment operator
SparseMatrix& SparseMatrix::operator =(SparseMatrix& A)
{
	// Let resize dispose of old data
	resize(0,0);

	// Copy in the new data
	copyFrom(A);

	// Mandatory return (of current object)
	return *this;
}

// Change the size of the matrix. diagValue is inserted in any
// new diagonal elements created. Old data is preserved if possible.
void SparseMatrix::resize(int dim1, int dim2, double diagValue)
{
	int			i,j,k;
	int			newDiagSize = dim1<dim2 ? dim1 : dim2;
	int			oldDiagSize = _dim1<_dim2 ? _dim1 : _dim2;

	double*		newDiag = NULL;
	int*		newRowHead = NULL;
	int*		newColHead = NULL;

	// Reallocate and copy the diagnonal array (if changed)
	if (newDiagSize != oldDiagSize) {
		if (newDiagSize>0) {
			newDiag = new double[newDiagSize];
		}
		for (i=0;i<newDiagSize && i<oldDiagSize; i++) {
			newDiag[i]=_diag[i];
		}
		for (;i<newDiagSize;i++) {
			newDiag[i] = diagValue;
		}
		delete[] _diag;
		_diag = newDiag;
	}

	// Reallocate and copy row list header and tail pointers
	if (dim1!=_dim1) {
		if (dim1>0) {
			newRowHead = new int[dim1];
		}
		for (i=0;i<dim1 && i<_dim1; i++) {
			newRowHead[i] = _rowHead[i];
		}
		for (;i<dim1;i++) {
			newRowHead[i] = -1;
		}
		delete[] _rowHead;
		_rowHead = newRowHead;
	}

	// Reallocate and copy col list header and tail pointers
	if (dim2!=_dim2) {
		if (dim2>0) {
			newColHead = new int[dim2];
		}
		for (j=0;j<dim2 && j<_dim2; j++) {
			newColHead[j] = _colHead[j];
		}
		for (;j<dim2;j++) {
			newColHead[j] = -1;
		}
		delete[] _colHead;
		_colHead = newColHead;
	}

	// If the size has been reduced, scan for obsolete entries
	// and flag them as free. While we are here, reorder the
	// free list in entry order. If size is reduced to zero,
	// free the data storage, otherwise reuse what is already owned.
	_freeHead = -1;
	if (dim1<_dim1 || dim2<_dim2) {

		// Do quick delete if nothing is saved
		if (dim1==0 || dim2==0) {
			delete[] _data;
			_data = NULL;
			_dataSize = 0;
			_dataCapacity = 0;
		}
		else {
			// Scan old data array and release obsolete entries
			for (k=_dataSize-1;k>=0;k--) {
				if (_data[k].row>dim1 || 
					_data[k].col>dim2 ||
					_data[k].row == -1) {
					// Free the entry
					_data[k].row = -1;
					_data[k].col = -1;
					_data[k].nextInRow = _freeHead;
					_freeHead = k;
				}
			}
		}
	}

	// Set the new dimensions
	_dim1 = dim1;
	_dim2 = dim2;

	// Clear any cached lookups
	_lastRef = -1;
}

// Clear the current matrix of off-diagonal data without
// releasing storage.
void SparseMatrix::clear(double diagValue)
{
	int		i,j;
	int		ds=diagSize();

	// Copy in diagonal values
	for (i=0;i<ds;i++) {
		_diag[i] = diagValue;
	}

	// Clear row and col heads
	for (i=0;i<_dim1;i++) {
		_rowHead[i]=-1;
	}
	for (j=0;j<_dim2;j++) {
		_colHead[j]=-1;
	}

	// Discard data entries and clear any cached data
	_dataSize = 0;
	_freeHead = -1;
	_lastRef = -1;
}

// Swap data between the current matrix and the one provided.
// A tedious process but avoids reallocating and copying data arrays.
void SparseMatrix::swapWith(SparseMatrix& B)
{
	SparseMatrix&	A = *this;	// The current matrix 
	SparseMatrix	C;			// Temporary work area

	// Copy data and array pointers from B to C
	C._dim1 = B._dim1;
	C._dim2 = B._dim2;
	C._dataSize = B._dataSize;
	C._dataCapacity = B._dataCapacity;
	C._paddingFactor =B._paddingFactor;
	C._freeHead = B._freeHead;
	C._lastRef = B._lastRef;
	C._rowHead = B._rowHead;
	C._colHead = B._colHead;
	C._diag = B._diag;
	C._data = B._data;

	// Copy data and array pointers from A to B
	B._dim1 = A._dim1;
	B._dim2 = A._dim2;
	B._dataSize = A._dataSize;
	B._dataCapacity = A._dataCapacity;
	B._paddingFactor =A._paddingFactor;
	B._freeHead = A._freeHead;
	B._lastRef = A._lastRef;
	B._rowHead = A._rowHead;
	B._colHead = A._colHead;
	B._diag = A._diag;
	B._data = A._data;

	// Copy data and array pointers from C to A
	A._dim1 = C._dim1;
	A._dim2 = C._dim2;
	A._dataSize = C._dataSize;
	A._dataCapacity = C._dataCapacity;
	A._paddingFactor =C._paddingFactor;
	A._freeHead = C._freeHead;
	A._lastRef = C._lastRef;
	A._rowHead = C._rowHead;
	A._colHead = C._colHead;
	A._diag = C._diag;
	A._data = C._data;
}

// Perform a deep copy of the matrix A into the current one.
// All data in the current matrix must be released before
// the copy is done.
void SparseMatrix::copyFrom(const SparseMatrix& B)
{
	int		i,j,k;

	// Copy basic data
	_dim1 = B._dim1;
	_dim2 = B._dim2;
	_dataSize = B._dataSize;
	_dataCapacity = B._dataCapacity;
	_paddingFactor =B._paddingFactor;
	_freeHead = B._freeHead;

	_lastRef = B._lastRef;

	// Initialize pointers to NULL
	_diag = NULL;
	_rowHead = NULL;
	_colHead = NULL;
	_data = NULL;

	// Allocate and copy necessary arrays (where size>0)
	if (diagSize()>0) {
		_diag = new double[diagSize()];
		for (k=0;k<diagSize();k++) {
			_diag[k]=B._diag[k];
		}
	}
	if (_dim1>0) {
		_rowHead = new int[_dim1];
		for (i=0;i<_dim1;i++) {
			_rowHead[i] = B._rowHead[i];
		}
	}
	if (_dim2>0) {
		_colHead = new int[_dim2];
		for (j=0;j<_dim2;j++) {
			_colHead[j] = B._colHead[j];
		}
	}
	if (_dataCapacity>0) {
		_data = new SparseEntry[_dataCapacity];
		// Only copy entries with data (even if freed)
		for (k=0;k<_dataSize;k++) {
			_data[k]=B._data[k];
		}
	}
}

// Access an entry from the matrix referenced via an iterator.
// Iteratively incrementing n gives access to the entire matrix.
// Return true if a good entry was found.
bool SparseMatrix::nextEntry(int& n, int& row,int& col, double& value)
{
	int		k;
	int		ds = diagSize();

	// First do the diagonals
	if (n<ds) {
		row = n;
		col = n;
		value = _diag[n];
		n++;
		return true;
	}

	// Otherwise iterate through data elements in no special order.

	// Find the next data element that is not deleted
	for (k=n-ds; k<_dataSize && _data[k].row==-1; k++) {}

	// If a good entry was found, return it.
	// Otherwise, that is all there is.
	if (k<_dataSize) {
		row = _data[k].row;
		col = _data[k].col;
		value = _data[k].value;
		n=k+1+ds;
		return true;
	}
	else {
		return false;
	}
}

// Provide access to entries by row. n should initially be set
// to 0. After each call, the next row entry is returned in
// order of increasing row. Return true if an entry was returned
// and false if not.
bool SparseMatrix::nextInRow(int& n, int row, int& col, double& val)
{
	int		cur,next;


	// Bounds check the row value
	if (row>=_dim1) {
		FatalError("(SparseMatrix::nextInRow) row is out of bounds");
		return false;
	}

	// n would be the index of the next entry returned, 
	// except that _data[0] is valid and n==0 is reserved for 
	// the start-up case. We use n=k+2, where k is the index 
	// to be returned. If a diagonal is to be returned next,
	// the sign of n is reversed, and the index of the next 
	// entry to return is indicated via -n=next+2.
	// When the list is exhausted, n is set to 1.

	// Handle special cases
	if (n<2) {
		// Should we stop now
		if (n==1) {
			return false;
		}

		// See if this is the first call for the row
		// and if so, set next accordingly
		if (n==0) {
			cur = _rowHead[row];
			if (cur == -1 || _data[cur].col>row) {
				// Set up to return diagonal now
				n = -(cur+2);
			}
		}
		else {
			// Otherwise decode n the usual way
			cur = -n-2;
		}
		// See if a diagonal being returned now
		// and if so return it.
		if (n<0) {
			col = row;
			val = _diag[row];
			n = cur+2;
			return true;
		}
	}
	else {
		cur = n-2;
	}

	// Return the entry indicated
	col = _data[cur].col;
	val = _data[cur].value;

	// See if the next entry crosses the diagonal
	next = _data[cur].nextInRow;
	if (col<row && (next==-1 || _data[next].col>row)) {
		n = -(next+2);
	}
	else {
		n = next+2;
	}
	return true;	
}


// Provide access to entries by column. n should initially be set
// to 0. After each call, the next column entry is returned in
// order of increasing column. Return true if an entry was returned
// and false if not.
bool SparseMatrix::nextInCol(int& n, int& row, int col, double& val)
{
	int	cur,next;

	// Bounds check the col value
	if (col>=_dim2) {
		FatalError("(SparseMatrix::nextInCol) col is out of bounds");
		return false;
	}

	// n would be the index of the next entry returned, 
	// except that _data[0] is valid and n==0 is reserved for 
	// the start-up case. We use n=k+2, where k is the index 
	// to be returned. If a diagonal is to be returned next,
	// the sign of n is reversed, and the index of the next 
	// entry to return is indicated via -n=next+2.
	// When the list is exhausted, n is set to 1.

	// Handle special cases
	if (n<2) {
		// Should we stop now
		if (n==1) {
			return false;
		}

		// See if this is the first call for the row
		// and if so, set next accordingly
		if (n==0) {
			cur = _colHead[col];
			if (cur == -1 || _data[cur].row>col) {
				// Set up to return diagonal now
				n = -(cur+2);
			}
		}
		else {
			// Otherwise decode n the usual way
			cur = -n-2;
		}
		// See if a diagonal being returned now
		// and if so return it.
		if (n<0) {
			row = col;
			val = _diag[col];
			n = cur+2;
			return true;
		}
	}
	else {
		cur = n-2;
	}

	// Return the entry indicated
	row = _data[cur].row;
	val = _data[cur].value;

	// See if the next entry crosses the diagonal
	next = _data[cur].nextInCol;
	if (row<col && (next==-1 || _data[next].row>col)) {
		n = -(next+2);
	}
	else {
		n = next+2;
	}
	return true;	
}

// Access a value at an entry specified by row and column.
// If there is no such entry, return 0.
double SparseMatrix::at(int row, int col)
{
	static char	boundsMsg[] = "(SparseMatrix::at) Index is out of range for matrix size";

	int		k;

	// Check row and col ranges
	if (row<0 || row>=_dim1 || col<0 || col>=_dim2) {
		FatalError(boundsMsg);
		return 0;
	}

	// Check for a diagonal entry
	if (row==col) {
		return _diag[row];
	}

	// Pick smaller of row or column to choose which list to search
	if (row<=col) {
		// Search the row list for a matching entry
		k = _rowHead[row];
		while (k!=-1 && _data[k].col<col) {
			k = _data[k].nextInRow;
		}
		if (k!=-1 && _data[k].col == col) {
			return _data[k].value;
		}
	}
	else {
		// Search the col list for a matching entry
		k = _colHead[col];
		while (k!=-1 && _data[k].row<row) {
			k = _data[k].nextInCol;
		}
		if (k!=-1 && _data[k].row == row) {
			return _data[k].value;
		}
	}

	// If nothing found, return 0
	return 0;
}

// Get a reference to the value of an entry specified by 
// row and column. If there is no such entry, create one.
// Optimize adding entries by increasing col within row.
double& SparseMatrix::refAt(int row, int col)
{
	static char	boundsMsg[] = "(SparseMatrix::refAt) Index is out of range for matrix size";
	static double unusedDouble = 0;


	// Make sure row and col are in a valid range
	if (row<0 || row>=_dim1 || col<0 || col>=_dim2) {
		FatalError(boundsMsg);
		return unusedDouble; // Never gets here, but keep the compiler happy anyway
	}

	// Check for diagonal entries and quickly handle them
	if (row==col) {
		return _diag[row];
	}

	int		curForRow = _rowHead[row];
	int		curForCol = _colHead[col];
	int		prevForRow = -1;
	int		prevForCol = -1;

	// By default, pick whichever search direction probably is shorter
	bool searchByRow = (col<row);

	// See if the previous lookup can serve as a starting point
	if (_lastRef != -1) {

		// Search by either row or col if either is the same as the
		// previously accessed reference and the forward list is not
		// already passed the desired index.
		if (_data[_lastRef].row==row && _data[_lastRef].col < col) {
			prevForRow = _lastRef;
			curForRow = _data[_lastRef].nextInRow;
			searchByRow = true;
		}
		else if (_data[_lastRef].col==col && _data[_lastRef].row < row) {
			prevForCol = _lastRef;
			curForCol = _data[_lastRef].nextInCol;
			searchByRow = false;
		}
	}

	// Do the search by either row or col as indicated
	if (searchByRow) {

		// Search over the row list and stop if match found
		while (curForRow != -1 && _data[curForRow].col<=col) {

			// Check for match on column and return if matched
			if (_data[curForRow].col==col) {
				_lastRef = curForRow;
				return _data[curForRow].value;
			}

			// Otherwise go on to the next in the row list
			prevForRow=curForRow;
			curForRow = _data[curForRow].nextInRow;
		}
	}
	else {
		// Search over the row list and stop if match found
		while (curForCol != -1 && _data[curForCol].row<=row) {

			// Check for match on row and return if matched
			if (_data[curForCol].row==row) {
				_lastRef = curForCol;
				return _data[curForCol].value;
			}

			// Otherwise go on to the next in the row list
			prevForCol=curForCol;
			curForCol = _data[curForCol].nextInCol;
		}
	}

	// If we get here, then there was no match and a new entry is needed.
	int newEnt = addDataEntry();

	// Splice the new entry into the row list
	while (curForRow != -1 && _data[curForRow].col<col) {
		prevForRow=curForRow;
		curForRow = _data[curForRow].nextInRow;
	}
	if (prevForRow == -1) {
		_data[newEnt].nextInRow = _rowHead[row];
		_rowHead[row]=newEnt;
	}
	else {
		_data[newEnt].nextInRow = _data[prevForRow].nextInRow;
		_data[prevForRow].nextInRow = newEnt;
	}

	// Splice the new entry into the column list
	while (curForCol != -1 && _data[curForCol].row<row) {
		prevForCol=curForCol;
		curForCol = _data[curForCol].nextInCol;
	}
	if (prevForCol == -1) {
		_data[newEnt].nextInCol = _colHead[col];
		_colHead[col]=newEnt;
	}
	else {
		_data[newEnt].nextInCol = _data[prevForCol].nextInCol;
		_data[prevForCol].nextInCol = newEnt;
	}

	// Fill in the rest of the new entry data
	_data[newEnt].row = row;
	_data[newEnt].col = col;
	_data[newEnt].value = 0;

	// Cache the entry just added for lookup
	_lastRef = newEnt;

	// Return a reference to the value
	return _data[newEnt].value;
}

// Set an entry to zero. If a data entry was allocated at
// this location it is freed.
void SparseMatrix::zeroAt(int row, int col)
{
	static char	boundsMsg[] = "(SparseMatrix::zeroAt) Index is out of range for matrix size";

	int				cur,prev;

	// Make sure row and col are in a valid range
	if (row<0 || row>=_dim1 || col<0 || col>=_dim2) {
		FatalError(boundsMsg);
		return;
	}

	// Test for a diagonal entry, which is never removed
	if (row==col) {
		_diag[row]=0;
		return;
	}

	// Remove from row list
	prev = -1;
	cur = _rowHead[row];
	while (cur!=-1 && _data[cur].col<col) {
		prev = cur;
		cur = _data[cur].nextInRow;
	}
	if (cur!=-1 && _data[cur].col == col) {
		// Found a match
		if (prev == -1) {
			_rowHead[row] = _data[cur].nextInRow;
		}
		else {
			_data[prev].nextInRow = _data[cur].nextInRow;
		}
	}
	else {
		// Did not find a match
		return;
	}

	// Remove from col list
	prev = -1;
	cur = _colHead[col];
	while (cur!=-1 && _data[cur].row<row) {
		prev = cur;
		cur = _data[cur].nextInCol;
	}
	if (cur!=-1 && _data[cur].row == row) {
		// Found a match
		if (prev == -1) {
			_colHead[col] = _data[cur].nextInCol;
		}
		else {
			_data[prev].nextInCol = _data[cur].nextInCol;
		}
	}
	else {
		// Did not find a match -- lists are broken
		FatalError("(SparseMatrix::zeroAt) Column list inconsistent with row list");
		return;
	}

	// Clear the current entry and put it in the free list
	_data[cur].value = 0;
	_data[cur].row = -1;
	_data[cur].col = -1;
	_data[cur].nextInCol = -1;
	_data[cur].nextInRow = _freeHead;
	_freeHead = cur;

	// If necessary, clear the lookup cache
	if (cur == _lastRef) {
		_lastRef = -1;
	}
}

// Set the data capacity to a new and possibly higher value
void SparseMatrix::capacity(int newCap) {

	int				k;
	SparseEntry*	newData;

	// See if additional capacity is needed.
	if (_dataCapacity<newCap) {

		// Allocate the new data array
		newData = new SparseEntry[newCap];

		// Copy data from old array to new
		for (k=0; k<_dataSize; k++) {
			newData[k]=_data[k];
		}

		// Update the data pointer and dispose of the old array
		delete[] _data;
		_data = newData;
		_dataCapacity = newCap;
	}
}

// Expand data to hold a new entry and return the index.
int SparseMatrix::addDataEntry()
{
	int				k, pad, newCap, maxCap;

	// See if there is an element in the free list
	// If there is, take it off the list and allocate it.
	if (_freeHead != -1) {
		k = _freeHead;
		_freeHead = _data[_freeHead].nextInRow;
		return k;
	}

	// If necessary, make room for the new entry
	if (_dataSize+1>=_dataCapacity) {

		// Compute a new capacity counting the pad factor
		// and rounding up to a multiple of 64 entries.
		// However, there is a maximum capacity needed for a
		// full matrix and there is no point going beyond this.
		pad = int( _dataSize*paddingFactor() );
		newCap = 64*((_dataSize+pad)/64+1);
		maxCap = dim1()*dim2()-diagSize();
		capacity(newCap<maxCap ? newCap : maxCap);
	}

	// Return the location of the new entry allocated
	return _dataSize++;
}


// Return the 1-norm of the matrix, that is, maximum column
// sum of entry absolute values.
double SparseMatrix::norm1()
{
	int		ds=diagSize();
	int		j,colIdx;
	double	maxSum = 0;
	double	colSum;

	for (j=0;j<_dim2;j++) {

		// Start with diagonal if it exists in this column
		if (j<ds) {
			colSum = fabs(_diag[j]);
		}
		else {
			colSum = 0;
		}

		// Add up absolute values in column
		colIdx = _colHead[j];
		while (colIdx!= -1) {
			colSum += fabs(_data[colIdx].value);
			colIdx = _data[colIdx].nextInCol;
		}

		// Set the maximum so far
		if (colSum>maxSum) {
			maxSum=colSum;
		}
	}

	return maxSum;
}

// Return the infinity-norm of the matrix, that is, maximum
// row sum of entry absolute values.
double SparseMatrix::normInf()
{
	int		ds=diagSize();
	int		i,rowIdx;
	double	maxSum = 0;
	double	rowSum;

	for (i=0;i<_dim1;i++) {

		// Start with diagonal if it exists in this row
		if (i<ds) {
			rowSum = fabs(_diag[i]);
		}
		else {
			rowSum = 0;
		}

		// Add up absolute values in row
		rowIdx = _rowHead[i];
		while (rowIdx!= -1) {
			rowSum += fabs(_data[rowIdx].value);
			rowIdx = _data[rowIdx].nextInRow;
		}

		// Set the maximum so far
		if (rowSum>maxSum) {
			maxSum=rowSum;
		}
	}

	return maxSum;
}

// Scale the current matrix in place
void SparseMatrix::scaleBy(double c) {

	int		n,k;
	int		ds = diagSize();

	// First scale the diagonal entries
	for (n=0;n<ds;n++) {
		_diag[n] *= c;
	}

	// Now scale the data entries
	for (k=0;k<_dataSize;k++) {
		_data[k].value *= c;
	}
}

// Add the current matrix to the one provided
// This computes B += A where A is the current matrix.
void SparseMatrix::addTo(SparseMatrix& B)
{
	int			n,k,row,col;
	int			ds = diagSize();

	// Make sure dimensions are correct
	if (B.dim1() != _dim1 || B.dim2() != _dim2) {
		FatalError("(SparseMatrix::addTo) Matrix dimensions do not match");
		return;
	}

	// Always add the matrix with the smallest number of
	// off diagonal terms to the one with more. If B
	// has few off diagonal terms, switch roles.
	if (_dataSize>B.dataSize() ) {

		SparseMatrix	Bcopy;
		
		B.swapWith(Bcopy);	// Quickly save the contents of B
		B= *this;			// Buffer the current matrix in B
		Bcopy.addTo(B);		// Do the add the other way around.
		return;
	}

	// Try to get enough storage allocated in B in case
	// it is starting from scratch. This is obviously a
	// minimum since entries may not correspond exactly.
	B.capacity(dataSize());

	// Accumulate in B starting with diagonals
	for (n=0;n<ds;n++) {
		B.refAt(n,n) += _diag[n];
	}

	// Now add off diagonal terms
	for (k=0;k<_dataSize;k++) {
		// Skip any delete entries
		if ((row=_data[k].row) != -1) {
			col = _data[k].col;
			B.refAt(row,col) += _data[k].value;
		}
	}
}

// Scale the current matrix and add it to the one provided
// This computes B += c*A where A is the current matrix.
void SparseMatrix::addScaledByTo(double c, SparseMatrix& B)
{
	int			n,k,row,col;
	int			ds = diagSize();

	// Make sure dimensions are correct
	if (B.dim1() != _dim1 || B.dim2() != _dim2) {
		FatalError("(SparseMatrix::addScaledByTo) Matrix dimensions do not match");
		return;
	}

	// Always add the matrix with the smallest number of
	// off diagonal terms to the one with more. If B
	// has few off diagonal terms, switch roles.
	if (_dataSize>B.dataSize() ) {

		SparseMatrix	Bcopy;
		
		B.swapWith(Bcopy);	// Quickly save the contents of B
		B= *this;			// Buffer the current matrix in B
		B.scaleBy(c);		// Scale it
		Bcopy.addTo(B);		// Do the add the other way around
		return;
	}

	// Try to get enough storage allocated in B in case
	// it is starting from scratch. This is obviously a
	// minimum since entries may not correspond exactly.
	B.capacity(dataSize());

	// Accumulate in B starting with diagonals
	for (n=0;n<ds;n++) {
		B.refAt(n,n) += c*_diag[n];
	}

	// Now add off diagonal terms
	for (k=0;k<_dataSize;k++) {
		// Skip any deleted entries
		if ((row=_data[k].row) != -1) {
			col = _data[k].col;
			B.refAt(row,col) += c*_data[k].value;
		}
	}
}

// Multiply a matrix times a column vector on the right.
// The vector is represented as an array of doubles.
// The equation is A*x=b where b is the output.
void SparseMatrix::rightMult(double* x, double* b)
{
	int		i,k;
	int		row,col,ds;

	// b will be used as an accumulator.
	// Start with the diagonal terms
	ds = diagSize();
	for (i=0;i<ds;i++) {
		b[i]=x[i]*_diag[i];
	}

	// Clear the rest of b (if any)
	for (;i<_dim1;i++) {
		b[i]=0;
	}

	// Now multiply by the off diagonal terms
	for (k=0;k<_dataSize;k++) {
		row = _data[k].row;
		col = _data[k].col;
		b[row] += _data[k].value * x[col];
	}
}

// Multiply a matrix time a row vector on the left.
// The vector is represented as an array of doubles.
// The equation is x*A=b when b is the output.
void SparseMatrix::leftMult(double* x, double* b)
{
	int		i,k;
	int		row,col,ds;

	// b will be used as an accumulator.
	// Start with the diagonal terms
	ds = diagSize();
	for (i=0;i<ds;i++) {
		b[i]=x[i]*_diag[i];
	}

	// Clear the rest of b (if any)
	for (;i<_dim2;i++) {
		b[i]=0;
	}

	// Now multiply by the off diagonal terms
	for (k=0;k<_dataSize;k++) {
		row = _data[k].row;
		col = _data[k].col;
		b[col] += _data[k].value * x[row];
	}
}

// Perform LU decomposition placing the combined results in 
// the current matrix. If the operation suceeds return true 
// and if not return false. tol specifies an absolute value 
// that can be used to trim the resulting LU decomposition
// resulting in an incomplete form. Pivoting is not done and 
// typically is not needed in these matrixes.
bool SparseMatrix::luDecompWithoutPivot(double abstol)
{
	const double		epsilon = 1e-20;
	
	bool				isSingular = false;
	int					i,j;
	int					colIdx,nextColIdx;
	int					rowIdx,rowIdxUpper;
	double				d,ratio;

	// Reduce column by column.
	// Note that the final column has no reduction to do.
	for (j=0;j<_dim2-1;j++) {

		// Test for zero on the diagonal. This can happen
		// since no pivoting is done. If a zero is found,
		// substitute a small value and keep going.
		if (_diag[j]!=0) {
			d = _diag[j];
		}
		else {
			isSingular = true;
			_diag[j] = epsilon;
			d = epsilon;
		}

		// Locate the starting index of that portion of
		// the row to the right (upper triangle) of the diagonal.
		rowIdxUpper = _rowHead[j];
		while (rowIdxUpper != -1 && _data[rowIdxUpper].col<j) {
			rowIdxUpper = _data[rowIdxUpper].nextInRow;
		}

		// Loop through rows with an entry in column
		// but skip entries above the diagonal
		colIdx=_colHead[j]; 
		while (colIdx!= -1) {

			// Locate entry in column to use in next pass
			nextColIdx = _data[colIdx].nextInCol;

			// Test if this entry is in the lower triangular part
			i = _data[colIdx].row;
			if (i>j) {

				// Get ratio A(i,j)/A(j,j) and test against tolerance.
				// If abstol>0, we are doing an incomplete LU decomp.
				ratio = _data[colIdx].value/d;

				// See if this row should be reduced or
				// should the value in A(i,j) just be set to 0.
				// This is only one way to do such a test and a
				// fairly naive way at that (might be improved later).
				if (fabs(ratio)>=abstol) {

					// Do a row reduction: row(i) -= r*row(j)
					// but skipping any embedded L entries.
					rowIdx = rowIdxUpper;

					while (rowIdx!=-1) {
						int k = _data[rowIdx].col;
						double val = _data[rowIdx].value;
						if (k>j && val!=0) {

							// Handle the diagonal case directly, but
							// otherwise let refAt do the dirty work.
							if (i==k) {
								// Diagonal case
								_diag[i] -= ratio*val;
							}
							else {
								// Off-diagonal cases
								refAt(i,k) -= ratio*val;
							}
						}
						rowIdx = _data[rowIdx].nextInRow;
					}

					// Store embedded L value in the entry just vacated
					_data[colIdx].value = ratio;
				}
				else {
					// Zero out the intersection value and continue.
					// Since the next entry is already found,
					// this can be done safely even though the column
					// list will be changed.
					zeroAt(i,j);
				}
			}

			// Move on to the next entry in the column
			colIdx = nextColIdx;
		}
	}

	// Return an indication if a zero diagonal was encountered or not.
	return !isSingular;
}

// Perform back substitution for an L/U combination matrix.
// The equation is A*x=b where x is the output and b the input.
// Note that permutation from pivoting is not supported here.
void SparseMatrix::luBackSubRight(double* x, double* b)
{
	int				i,j,rowIdx;
	double			sum;


	// Only square matrixes are supported for this
	if (_dim1 != _dim2) {
		FatalError("(SparseMatrix::luBackSubRight) Matrix must be square");
		return;
	}

	// Just in case this is a 0x0 matrix return now
	if (_dim1==0) {
		return;
	}

	// Copy data from b to x where the work is done
	// In this case, b and x can be the same vector.
	if (x!=b) {
		for (i=0;i<_dim1;i++) {
			x[i]=b[i];
		}
	}

	// Do backsubstitution using the L component
	for (i=1;i<_dim1;i++) {
		sum = x[i];
		rowIdx = _rowHead[i];
		while (rowIdx != -1 && (j=_data[rowIdx].col)<i) {
			sum -= _data[rowIdx].value*x[j];
			rowIdx = _data[rowIdx].nextInRow;
		}
		x[i]=sum;
	}

	// Do backsubstitution using the U component
	x[_dim1-1] /= _diag[_dim1-1];
	for (i=_dim1-2;i>=0;i--) {
		sum = x[i];
		rowIdx = _rowHead[i];
		while (rowIdx != -1) {
			// Take only the upper part of the matrix
			j = _data[rowIdx].col;
			if (j>i) {
				sum -= x[j]*_data[rowIdx].value;
			}
			rowIdx = _data[rowIdx].nextInRow;
		}
		x[i] = sum/_diag[i];
	}
}

// Perform back substitution for an L/U combination matrix.
// The equation is x*A=b where x is the output and b the input.
// Note that permutation from pivoting is not supported here.
void SparseMatrix::luBackSubLeft(double* x, double* b)
{
	int				i,j,colIdx;
	double			sum;


	// Only square matrixes are supported for this
	if (_dim1 != _dim2) {
		FatalError("(SparseMatrix::luBackSubLeft) Matrix must be square");
		return;
	}

	// Just in case this is a 0x0 matrix return now
	if (_dim1==0) {
		return;
	}

	// Copy data from b to x where the work is done
	// In this case, b and x can be the same vector.
	if (x!=b) {
		for (i=0;i<_dim1;i++) {
			x[i]=b[i];
		}
	}

	// Do backsubstitution using the U component
	x[0] /= _diag[0];
	for (j=1;j<_dim1;j++) {
		sum = x[j];
		colIdx = _colHead[j];
		while (colIdx != -1 && (i=_data[colIdx].row)<j) {
			sum -= x[i]*_data[colIdx].value;
			colIdx = _data[colIdx].nextInCol;
		}
		x[j] = sum/_diag[j];
	}

	// Do backsubstitution using the L component
	for (j=_dim1-2;j>=0;j--) {
		sum = x[j];
		colIdx = _colHead[j];
		while (colIdx != -1) {
			i=_data[colIdx].row;
			if (i>j) {
				sum -= _data[colIdx].value*x[i];
			}
			colIdx = _data[colIdx].nextInCol;
		}
		x[j]=sum;
	}
}

// Solve the system A*x=b using the current matrix as A,
// b as the input, and placing the result in x. A is not changed.
// Return true if the operation succeeds and false if not.
bool SparseMatrix::luSolve(double* x, double* b)
{
	bool			isSingular;

	// Copy the current matrix to a work area
	SparseMatrix	LU_A(*this);

	// Do LU decomp and backsub
	isSingular = LU_A.luDecompWithoutPivot();
	LU_A.luBackSubRight(x,b);

	// Return the singularity status.
	// Note that values in x may still be useful
	// even if the current matrix is singular.
	return isSingular;
}

// Solve the system A*x=b using the current matrix as A
// and b as input placing the result in x. A is not changed.
// Return true if the operation succeeds and false if not.
// Preconditioned bi-conjugate gradient method is used to 
// solve the system. If usePinv is false, then P should be
// the LU decomposition of a matrix approximating A. If
// usePinv is true, then P approximates A^-1 instead
// and is not in LU decomposition format.
// Formulas are taken directly from Press et al.

bool SparseMatrix::pbcgSolve(
	double*	x,				// Solution output 
	double*	b,				// Desired value input
	double*	pErr,			// error (norm-inf of residual)
	int*	pIterUsed,		// iterations used
	double	tol,			// Absolute tolerance for convergence
	int		maxIter,		// Maximum iterations allowed (0 means _dim1)
	double*	x0,				// Initial guess for solution value
	SparseMatrix& P,		// Precondition matrix(s) or empty matrix
	bool	usePinv)		// P is approx A^-1 rather than A	
{
	int				N = _dim1;
	int				i;
	bool			useP = true;

	double*			workArea;
	double*			r;
	double*			rbar;
	double*			z;
	double*			zbar;
	double*			p;
	double*			pbar;
	double*			Ap;			// A*p
	double*			pbarA;		// pbar * A
	
	double			alpha,beta;
	double			rbarDotZAlpha,rbarDotZBeta,pbarDotAp;

	double			err = 0;
	int				iterUsed = 0;

	// Only square matrixes are supported for this
	if (_dim1 != _dim2) {
		FatalError("(SparseMatrix::bicgSolveAxb) Matrix must be square");
		return false;
	}

	// If the preconditioning matrix is empty, it is not used
	if (P.dim1()==0 && P.dim2()==0) {
		useP = false;
	}
	else {
		// Otherwise, preconditioning matrix must be the same size
		if (_dim1!= P.dim1() || _dim2 != P.dim2()) {
			FatalError("(SparseMatrix::bicgSolveAxb) Matrix sizes do not match");
			return false;
		}
	}

	// If x0 is not NULL, copy x0 to x to start
	if (x0!=NULL) {
		for (i=0;i<N;i++) {
			x[i]=x0[i];
		}
	}
	else {
		// Otherwise, start with x=b (A is close to unity)
		for (i=0;i<N;i++) {
			x[i]=b[i];
		}
	}

	// Get a work area for all vectors
	workArea = new double[8*N];

	// Subdivide work area into arrays.
	// Might be able to improve cache hits for large
	// N if r,z,p, etc all allocated as part of a
	// common structure such that for each i, 
	// r[i] and z[i] etc. are located close together.
	r			= workArea;
	z			= r+N;
	p			= z+N;
	Ap			= p+N;
	rbar		= Ap+N;
	zbar		= rbar+N;
	pbar		= zbar+N;
	pbarA		= pbar+N;

	// Get the residual, r=b-A*x0
	rightMult(x,r);
	err = 0;
	for (i=0;i<N;i++) {
		rbar[i] = r[i] = b[i]-r[i];
		if (fabs(r[i])>err) {
			err = fabs(r[i]);
		}
	}

	// Was the first x a lucky guess -- if so stop now
	if (err<=tol) {
		delete[] workArea;
		if (pErr!=NULL) *pErr = err;
		if (pIterUsed!=NULL) *pIterUsed = iterUsed;
		return true;
	}

	if (useP) {
		// Get initial z values via the preconditioning matrix
		if (usePinv) {
			// In this case, we can get z and zbar directly
			P.rightMult(r,z);
			P.leftMult(rbar,zbar);
		}
		else {
			// Solve P*z=r and zbar*P=rbar
			P.luBackSubRight(z,r);
			P.luBackSubLeft(zbar,rbar);
		}
	}
	else {
		// Use diagonal values for preconditioning
		for (i=0;i<N;i++) {
			if (_diag[i]!=0) {
				z[i] = r[i]/_diag[i];
				zbar[i] = rbar[i]/_diag[i];
			}
			else {
				z[i] = r[i];
				zbar[i] = rbar[i];
			}
		}
	}

	// And then initialize the arrays
	for (i=0;i<N;i++) {
		rbar[i] = r[i];
		p[i]=z[i];
		pbar[i]=zbar[i];
	}

	// Do iterations until the error tolerance is met
	// or the maximum number of iterations is exceeded.
	iterUsed = 1;
	while (maxIter<=0 || iterUsed<=maxIter) {

		// Compute alpha
		rbarDotZAlpha=0;
		for (i=0;i<N;i++) {
			rbarDotZAlpha += rbar[i]*z[i];
		}
		rightMult(p,Ap);
		pbarDotAp = 0;
		for (i=0;i<N;i++) {
			pbarDotAp += pbar[i]*Ap[i];
		}
		alpha = rbarDotZAlpha/pbarDotAp;

		// Accumulate the result in x
		for (i=0;i<N;i++) {
			x[i] += alpha*p[i];
		}

		// Get next r and check for convergence
		err = 0;
		for (i=0;i<N;i++) {
			r[i] -= alpha*Ap[i];
			if (fabs(r[i])>err) {
				err = fabs(r[i]);
			}
		}
		if (err<=tol) {
			// Stop the loop now - we are done
			break;
		}

		// Get the next value rbar
		leftMult(pbar,pbarA);
		for (i=0;i<N;i++) {
			rbar[i] -= alpha*pbarA[i];
		}

		if (useP) {
			// Get next z values via the preconditioning matrix
			if (usePinv) {
				// In this case, we can get z and zbar directly
				P.rightMult(r,z);
				P.leftMult(rbar,zbar);
			}
			else {
				// Solve P*z=r and zbar*P=rbar
				P.luBackSubRight(z,r);
				P.luBackSubLeft(zbar,rbar);
			}
		}
		else {
			// Use diagonal values for preconditioning
			for (i=0;i<N;i++) {
				if (_diag[i]!=0) {
					z[i] = r[i]/_diag[i];
					zbar[i] = rbar[i]/_diag[i];
				}
				else {
					z[i] = r[i];
					zbar[i] = rbar[i];
				}
			}
		}

		// Compute beta = (rbar(next)*z(next))/(rbar*z)
		rbarDotZBeta = 0;
		for (i=0;i<N;i++) {
			rbarDotZBeta += rbar[i]*z[i];
		}
		beta = rbarDotZBeta/rbarDotZAlpha;

		// Get new values of p and pbar
		for (i=0;i<N;i++) {
			p[i]=z[i]+beta*p[i];
			pbar[i]=zbar[i]+beta*pbar[i];
		}

		// Count this iteration
		iterUsed++;
	}

	// Release allocated storage
	delete[] workArea;

	// Return results (other than x)
	if (pErr!=NULL) *pErr = err;
	if (pIterUsed!=NULL) *pIterUsed = iterUsed;
	return (err<=tol);
}


// ====================================================================
// SparseRowReference class body
// ====================================================================



// Return a reference to a row.
// This cannot be done inline due to compiler restrictions.
SparseRowRef SparseMatrix::operator [] (int row)
{
	return SparseRowRef(this,row);
}



// ====================================================================
// Interpolator class body
// ====================================================================



// Constructor
Interpolator::Interpolator()
{
	// Set known starting values
	_nData = 0;
	_xData = NULL;
	_yData = NULL;
	_yMin = -numeric_limits<double>::infinity();
	_yMax = numeric_limits<double>::infinity();
	_xMin = -numeric_limits<double>::infinity();
	_xMax = numeric_limits<double>::infinity();
}

// Destructor
Interpolator::~Interpolator() 
{
	// Delete allocated data arrays
	delete[] _xData;
	delete[] _yData;
}

// Prepare to load data into arrays
void Interpolator::preDataLoad()
{
	int n = _nData;

	// Check that there is enough data to do interpolation
	if (n<2) {
		FatalError("(Interpolator::preDataLoad) Insufficient data provided");
	}

	// Delete any old values
	delete[] _xData;
	delete[] _yData;

	// Allocate arrays
	_xData = new double[n];
	_yData = new double[n];

	// Clear results of any prior searches
	_lastFound = 0;
}

// Copy data from two arrays of values
void Interpolator::data(int n, double* x, double* y)
{
	int			k;

	// Allocate arrays
	_nData = n;
	preDataLoad();

	// Copy values
	for (k=0;k<n;k++) {
		_xData[k] = x[k];
		_yData[k] = y[k];
	}

	// Put things in order
	sortData();

	// Let subclass finish up with load
	postDataLoad();
}

// Copy data from n x 2 matrix of values
void Interpolator::data(int n, double* xy[2])
{
	int			k;

	// Allocate arrays
	_nData = n;
	preDataLoad();

	// Copy values
	for (k=0;k<n;k++) {
		_xData[k] = xy[k][0];
		_yData[k] = xy[k][1];
	}

	// Put things in order
	sortData();

	// Let subclass finish up with load
	postDataLoad();
}

// Sort data into ascending sequence by x value
void Interpolator::sortData()
{
	int			i,j;
	double		x,y;

	// Use a garden variety insertion sort
	for (i=1;i<_nData;i++) {

		// Save data to be inserted
		x=_xData[i];
		y=_yData[i];

		// Move entries up to make room for insert
		for (j=i;j>0 && x<_xData[j-1];j--) {
			_xData[j]=_xData[j-1];
			_yData[j]=_yData[j-1];
		}

		// Insert data at the right point
		_xData[j]=x;
		_yData[j]=y;
	}
}

// Do the interpolation and enforce range limits
double Interpolator::yAtX(double x)
{
	// Make sure data was loaded first
	if (_xData == NULL) {
		FatalError("(Interpolator::yAtX) no data loaded");
	}

	// Check for valid x value ranges
	if (x<xMin() || x>xMax() ) {
		FatalError("(Interpolator::yAtX) x value is out of range");
	}

	// Let subclass handle the real work of extrapolation
	double y = rawYAtX(x);

	// Enforce limits on y
	if (y<yMin() ) {
		y = yMin();
	}
	else if (y>yMax() ) {
		y = yMax();
	}

	return y;
}

// Look up the value x in the data table and return an associated
// index of the largest value smaller than x. Return -1 if no
// such value is found. 
int Interpolator::findX(double x)
{
	int			il,ih,im;

	// Check to see if lastFound still works.
	// This speeds up sequential access to the table.
	if (_lastFound == -1) {
		if (x<_xData[0]) {
			return -1;
		}
	}
	else if (_xData[_lastFound]<= x &&
		( _lastFound+1 == _nData || x < _xData[_lastFound+1]) ) {
		return _lastFound;
	}

	// Otherwise, do a binary search
	il = 0;
	ih = _nData;
	while (ih-il>1) {
		im = (il+ih)/2;
		if (x < _xData[im]) {
			ih=im;
		}
		else {
			il= im;
		}
	}

	// Check for case x is below range of data
	if (il==0 && x < _xData[0]) {
		return _lastFound = -1;
	}
	
	// Otherwise, il points to entry <= to x
	return _lastFound = il;
}



// ====================================================================
// Linear interpolator class body
// ====================================================================



// Constructors and destructor
LinearInterpolator::LinearInterpolator() {}

LinearInterpolator::LinearInterpolator(
	int			n,		// number of values
	double*		x,		// x data array
	double*		y)		// y data array
{
	data(n,x,y);
}

LinearInterpolator::LinearInterpolator(
	int			n,		// number of values
	double*		xy[2])	// data matrix
{
	data(n,xy);
}

LinearInterpolator::~LinearInterpolator() {}

// Interpolate a value
double LinearInterpolator::rawYAtX(double x)
{
	// Look up x in the data table
	int i0 = findX(x);
	int	i1 = i0+1;

	// Adjust for extrapolation cases
	if (i0 < 0) {
		i0=0;
		i1=1;
	}
	else if (i1>=_nData) {
		i0--;
		i1--;
	}

	// Simplify notation below
	double y0 = _yData[i0];
	double y1 = _yData[i1];
	double x0 = _xData[i0];
	double x1 = _xData[i1];

	// Return the interpolated value as linear fit
	return y0+(y1-y0)/(x1-x0)*(x-x0);
}



// ====================================================================
// Cubic interpolator class body
// ====================================================================



// Constructors and destructor
CubicInterpolator::CubicInterpolator() {}

CubicInterpolator::CubicInterpolator(
	int			n,		// number of values
	double*		x,		// x data array
	double*		y)		// y data array
{
	_yPrime = NULL;
	data(n,x,y);
}

CubicInterpolator::CubicInterpolator(
	int			n,		// number of values
	double*		xy[2])	// data matrix
{
	_yPrime = NULL;
	data(n,xy);
}

CubicInterpolator::~CubicInterpolator() 
{
	delete[] _yPrime;
}

// Compute yPrime values at each data point
void CubicInterpolator::postDataLoad()
{
	int			i;
	int			n = _nData-1;

	// Dispose of any prior values
	delete[] _yPrime;

	// Will default to linear interpolation with only two points
	if (_nData<3)
		return;

	// Allocate the array of values
	_yPrime = new double[_nData];

	// Handle first and last cases specially
	_yPrime[0] = yDotFit(
		_xData[1],		_yData[1],
		_xData[0],		_yData[0],
		_xData[2],		_yData[2] );
	_yPrime[n] = yDotFit(
		_xData[n-1],	_yData[n-1],
		_xData[n],		_yData[n],
		_xData[n-2],	_yData[n-2] );

	// Load the rest of the table
	for (i=1;i<n;i++) {
		_yPrime[i] = yDotFit(
			_xData[i-1],_yData[i-1],
			_xData[i],	_yData[i],
			_xData[i+1],_yData[i+1] );
	}
}

// Interpolate a value by fitting a piecewise cubic polynomial
double CubicInterpolator::rawYAtX(double x)
{
	double		x0,x1;			// x values at end point of interval used
	double		y0,y1;			// y values corr to x
	double		yp0,yp1;		// y prime (dy/dx) values
	double		a,b,c;			// cubic polynomial coefficients
	double		z,r1,r2;		// work variables
	int			i,n;			// work variables

	// Handle degenerate case of only two points in data table so
	// that we don't have to test for it later.
	if (_nData==2) {
		return _yData[0]+(x-_xData[0])*
			(_yData[1]-_yData[0])/(_xData[1]-_xData[0]);
	}

	// Look up x in the data table
	n = _nData-1;
	i = findX(x);

	// Locate the data points containing x or if none, those closest
	if (i<1) {
		i = 0;
	}
	else if (i==n) {
		i = n-1;
	}
	x0 = _xData[i];
	y0 = _yData[i];
	yp0 = _yPrime[i];
	x1 = _xData[i+1];
	y1 = _yData[i+1];
	yp1 = _yPrime[i+1];

	// Fit a cubic equation to both y and dy values at the ends of
	// the interval containing x. To simplify the fit, a change of
	// variables from x to z (see below) is used.

	z = (x-x0)/(x1-x0);			// z in [0,1]
	c = yp0*(x1-x0);			// c = dy/dz at x=x0
	r1 = y1-(c+y0);				// r1 = a+b
	r2 = (yp1-yp0)*(x1-x0);		// r2 = 3a+2b
	b = 3*r1-r2;
	a = r1-b;

	// Evaluate and return the polynomial using the coefficients found
	return a*z*z*z + b*z*z + c*z + y0;
} 


// Estimate the derivative of a curve y=f(x) given three pairs of x,y
// by fitting a quadratic equation through the points given. 
// Return an estimate of f' at x=x1.
double CubicInterpolator::yDotFit(
		double x0, double y0,
		double x1, double y1,
		double x2, double y2)
{
	double		dy;

	// Expand derivative of polynomial interpolation formula
	dy  = y0*(x1-x2)/((x0-x1)*(x0-x2));
	dy += y1*(2*x1-x0-x2)/((x1-x0)*(x1-x2));
	dy += y2*(x1-x0)/((x2-x1)*(x2-x0));

	return dy;
}



// ====================================================================
// RandomNumberGenerator class body
// ====================================================================



// Constructors and destructor
RandomNumberGenerator::RandomNumberGenerator() {}
RandomNumberGenerator::~RandomNumberGenerator() {}

// Get a random value as an integer (via subclass)
int RandomNumberGenerator::inext()
{
	FatalError("(RandomNumberGenerator::inext) function not provided by subclass");
	return 0;
}

// Get a random value as a double
double RandomNumberGenerator::next()
{
	// Convert integer value to double
	return inext();
}


// Return a probability density at x.
// This must be implemented by a subclass if used.
double RandomNumberGenerator::density(double x)
{
	FatalError("(RandomNumberGenerator::density) Subclass must implement");
	return 0;
}

// Return a probability density at k.
// Normally a subclass would implement this if it is used.
double RandomNumberGenerator::density(int k)
{
	return density(double(k)); // try to use double version
}

// Return the density at an adjacent point using an already
// computed density. This would typically be overridden if
// there is recursive formulation of the density.
// dir is either +1 or -1. The density at k+dir is returned.
double RandomNumberGenerator::adjDensity(int k, double denAtK, int dir)
{
	return density(k+dir);
}

// Return a mode value as a double for a continuous distribution.
// This must be implemented by a subclass if used.
// Normally a subclass would cache the value derived
// since it will not change for a given distribution.
double RandomNumberGenerator::mode()
{
	return imode();	// convert the integer value to a double
}

// Return a mode value as an integer for a discrete distribution.
// A subclass must implement this if it is used.
int RandomNumberGenerator::imode()
{
	FatalError("(RandomNumberGenerator::imode) Subclass must implement");
	return 0;
}

// Return the density value at the mode.
// Normally a subclass will cache this value if used.
double RandomNumberGenerator::modeDensity()
{
	return density(mode());
}



// ====================================================================
// RandomAlgorithm class body
// ====================================================================



// Constructor and destructor
RandomAlgorithm::RandomAlgorithm(NonuniformRandom* base) {
	_base = base;
}

RandomAlgorithm::~RandomAlgorithm() {}

// Return a random value as a double (continuous distribution)
double RandomAlgorithm::next()
{
	FatalError("(RandomAlgorithm::next) Subclass must implement");
	return 0;
}

// Return a random value as an integer (discreet distribution)
int RandomAlgorithm::inext()
{
	FatalError("(RandomAlgorithm::inext) Subclass must implement");
	return 0;
}



// ====================================================================
// InverseCDF class body
// ====================================================================



// Constructor and destructor
InverseCDF::InverseCDF(NonuniformRandom* base, int maxSize, int minValue)
: RandomAlgorithm(base)
{
	_maxSize = maxSize;
	_minValue = minValue;
}

InverseCDF::~InverseCDF() {}

// Do an inverse CDF lookup and return the index of the
// corresponding entry. If the table cannot be expanded
// further, return minValue-1 to indicate a lookup failure.
// The density function is supplied by the base object.
int InverseCDF::inext()
{
	int			k,kmax,m;
	double		u,cum;

	// First make sure _cdf is not empty
	if (_cdf.size() == 0) {
		_cdf.push_back(_base->density(_minValue)); // Add first entry
	}

	// Generate the uniform value for inverse cdf
	u = base()->runif();
	cum = _cdf.back();

	// See if we are adding new entries or else
	// searching the existing cdf values. This search
	// assumes that the table is relatively small.
	if (u<=cum) {

		m = _cdf.size()/2;	// Locate the middle of the table

		// See where we are with respect to the middle of the table.
		// For many distributions, much of the total density is near here. 

		if (m<_cdf.size() && u>=_cdf[m]) {

			// Search forwards starting at the middle
			for (k=m+1;k<_cdf.size()-1;k++) {
				if (u < _cdf[k]) {
					return k+_minValue;
				}
			}
			return k+_minValue;
		}
		else {
			// Search backwards starting at the middle
			kmax = m < _cdf.size() ? m : _cdf.size();
			for (k=m;k>0;k--) {
				if (u >= _cdf[k-1]) {
					return k+_minValue;
				}
			}
			return _minValue;
		}

	}
	else {
		// Expand the cdf table as necessary.
		for (k=_cdf.size();k<=_maxSize;k++) {
			cum=cum+_base->density(k);
			_cdf.push_back(cum);
			if (u<cum) {
				return k+_minValue;
			}
		}
	}

	// If we get here, there was a lookup failure
	// that is, the cdf table cannot be expanded further
	// and the random value exceeds the existing table.
	return _minValue-1;
}

// Return the size of the CDF table
int InverseCDF::cdfSize()
{
	return _cdf.size();
}

// Return the last entry in the CDF table
double InverseCDF::cdfMax()
{
	return _cdf[cdfSize()-1];
}



// ====================================================================
// RatioOfUniforms class body
// ====================================================================



// Constructor and destructor
RatioOfUniforms::RatioOfUniforms(NonuniformRandom* base)
: RandomAlgorithm(base)
{
	double pmu;	// probability density at the mode

	// Set initial values
	_mu = base->mode();
	pmu = base->modeDensity();
	_um = sqrt(pmu);
	_vmax = 1/_um;
	_vmin = -1/_um;

	// Set initial squeeze points.
	_squpl = _um/2;
	_squmn = _um/2;
	_sqvpl = 0;
	_sqvmn = 0;

	// Set precomputed squeeze ratios
	_sqr1pl = 0;
	_sqr2pl = 0;
	_sqr1mn = 0;
	_sqr2mn = 0;	
}

RatioOfUniforms::~RatioOfUniforms() {}

// Generate and return a random value for a continuous distribution.
// performance. The basic algorithm is SROUC from Leydold with some
// improvements specific to this implementation.
double RatioOfUniforms::next()
{
	double	u,v,r,x,den;
	double	tsqu,tsqv;

	// Sample u,v space for a value to return
	for(;;) {

		// Generate the random point in U,V space
		u = _um*base()->runif();
		v = _vmin+(_vmax-_vmin)*base()->runif();

		// Determine the offset from the mode value to get a
		// candidate value for the generated value.
		r = v/u;
		x = r+_mu;

		// Check u,v being inside the squeeze region if squeeze is used.
		if (v>=0) {
			if (_sqr1pl>r && _sqr2pl>v/(_um-u)) {
				return x;
			}
		}
		else {
			if (_sqr1mn<r && _sqr2mn<v/(_um-u)) {
				return x;
			}
		}

		// Since u,v is not in the squeeze region, get the density.
		// Quickly reject values of x that are out of range while
		// remembering the bounds for later.
		den=_base->density(x);
		if (den==0) {
			continue;
		}

		// Try to expand the squeeze region for next time.
		// If the expansion is successful, then readjust v bounds.
		// The rationale relies upon the squeeze region being 
		// made up of two triangles, one with v>=0 and the other v<0.

		tsqu = sqrt(den);
		tsqv = r*tsqu;
		if (r>0) {				
			if (tsqv>_sqvpl) {
				// Save the new squeeze triangle apex for v>=0
				_squpl = tsqu;
				_sqvpl = tsqv;

				// Recompute squeeze ratios
				_sqr1pl=_sqvpl/_squpl;
				_sqr2pl=_sqvpl/(_um-_squpl);

				// Since one squeeze region was increased in size,
				// the opposite v bound can be adjusted.
				// Note: this is the adjustment from SROUC when the
				// derivative of the density is not available. A better
				// adjustment can be done with some increase in complexity.
				_vmin = -1/_um + _sqvpl/2;
			}
		}
		else {
			if (tsqv<_sqvmn) {
				// Save the new squeeze triangle apex for v<0
				_squmn = tsqu;
				_sqvmn = tsqv;

				// Recompute squeeze ratios
				_sqr1mn=_sqvmn/_squmn;
				_sqr2mn=_sqvmn/(_um-_squmn);

				// Since one squeeze region was increased in size,
				// the opposite v bound can be adjusted.
				// Note: this is the adjustment from SROUC when the
				// derivative of the density is not available. A better
				// adjustment can be done with some increase in complexity.
				_vmax = 1/_um + _sqvmn/2;
			}
		}

		// Do the ROU accept/reject test for this u,v pair
		if (u*u<=den) {
			return x;
		}
	}
}

// Generate and return a random value for a discrete distribution.
// This code could be unified with next() but with some cost in
// performance. The basic algorithm is SROUD from Leydold with some
// improvements specific to this implementation.
int RatioOfUniforms::inext()
{
	double	u,v,r,den;
	double	tsqu,tsqv;
	double	aden,tadu,tadv,vtmp,slope;
	int		x,dx;

	// Sample u,v space for a value to return
	for(;;) {

		// Generate the random point in U,V space
		u = _um*base()->runif();
		v = _vmin+(_vmax-_vmin)*base()->runif();

		// Determine the integer offset from the mode value
		r = v/u;
		dx = int(floor(r));
		x = dx+int(_mu); // note: _mu must be integral

		// Check u,v being inside the squeeze region.
		if (v>=0) {
			if (_sqr1pl>r && _sqr2pl>v/(_um-u)) {
				return x;
			}
		}
		else {
			if (_sqr1mn<r && _sqr2mn<v/(_um-u)) {
				return x;
			}
		}

		// Since u,v is not in the squeeze region, get the density.
		// Quickly reject values that are out of range.
		den=_base->density(x);
		if (den==0) {
			continue;
		}

		// Try to expand the squeeze region for next time.
		// If the expansion is successful, then readjust v bounds.
		// The rationale relies upon the squeeze region being 
		// made up of two triangles, one with v>=0 and the other v<0.

		tsqu = sqrt(den);
		if (dx>0) {
			tsqv = dx*tsqu;	
			if (tsqv>_sqvpl) {
				// Save the new squeeze triangle apex for v>=0
				_squpl = tsqu;
				_sqvpl = tsqv;

				// Recompute squeeze ratios
				_sqr1pl=_sqvpl/_squpl;
				_sqr2pl=_sqvpl/(_um-_squpl);

				// Use the adjacent density to get the slope of
				// a tangent line passing through star points 
				// and see if it can improve the bounds on v.
				aden = _base->adjDensity(x,den,+1);
				tadu = sqrt(aden);
				tadv = (dx+2)*tadu;
				vtmp = (dx+1)*tsqu;
				slope = (tadv-vtmp)/(tadu-tsqu);
				if (slope>=0) 
					vtmp+=slope*(_um-tsqu);
				else
					vtmp-=slope*tsqv;
				if (vtmp<_vmax) {
					_vmax = vtmp;
				}
			}
		}
		else if (dx<-1) {
			tsqv = (dx+1)*tsqu;	// round to extreme case for floor()
			if (tsqv<_sqvmn) {
				// Save the new squeeze triangle apex for v<0
				_squmn = tsqu;
				_sqvmn = tsqv;

				// Recompute squeeze ratios
				_sqr1mn=_sqvmn/_squmn;
				_sqr2mn=_sqvmn/(_um-_squmn);

				// Use the adjacent density to get the slope of
				// a tangent line passing through star points 
				// and see if it can improve the bounds on v.
				// Note: this is an improvement over the adjustment
				// used in the published version of SROUD.
				aden = _base->adjDensity(x,den,-1);
				tadu = sqrt(aden);
				tadv = (dx-1)*tadu;
				vtmp = dx*tsqu;
				slope = (tadv-vtmp)/(tadu-tsqu);
				if (slope<=0) 
					vtmp+=slope*(_um-tsqu);
				else
					vtmp-=slope*tsqu;
				if (vtmp>_vmin) {
					_vmin = vtmp;
				}
			}
		}

		// Do the ROU accept/reject test for this u,v pair
		if (u*u<=den) {
			return x;
		}
	}
}



// ====================================================================
// UniformRandom class body
// ====================================================================



// Static values
unsigned int UniformRandom::_seedCount = 0;
UniformRandom* UniformRandom::_defaultGen = new MT19937_UniformRandom();

// Constructors and destructor
UniformRandom::UniformRandom() 
{
	_stride = 1;
	_offset = 0;
	_isSeeded = false;
}

UniformRandom::~UniformRandom() {}

// Set the default generator
void UniformRandom::defaultGen(UniformRandom* def)
{
	if (_defaultGen!=NULL) {
		delete _defaultGen;
	}
	_defaultGen = def;
}

// Set stride values
void UniformRandom::setStride(int stride,int offset)
{
	_stride = stride;
	_offset = offset;

	if (!_isSeeded ) {
		advance(offset);
	}
}

// Set the seed value based on an array of integers.
// Only some types of generators can set seeds by
// incorporating one integer at a time. Others would
// need to override this function.
void UniformRandom::setSeed(int nbr, unsigned int* seedValues)
{
	int k;

	// Clear any old seed
	resetSeed();
	_isSeeded=false;

	// Feed in new seed values one at a time
	// Note that addSeed does not imply commutativity
	for (k=0;k<nbr;k++) {
		addSeed(seedValues[k]);
	}

	// Set the flag the the seed was set
	// and advance to the current offset.
	_isSeeded = true;
	advance(_offset);
}

// Convenience function for setting the seed from one integer
void UniformRandom::setSeed(unsigned int seed) 
{
	const int		n=8;
	unsigned int	array[n];
	int				k;

	// Expand the array with a simple hash. This is useful
	// for breaking up obvious numerical patterns, especially
	// when successive seed values are used.
	array[0]=seed;
	for (k=1;k<=n;k++) {
		array[k%n]=array[k-1]<<23 ^ array[k-1]>>9;
		array[k%n] += 0xaaaaaaaaU;
	}

	// Pass values from the generated array
	setSeed(8,array);
}

// Set seed from a value in the interval (0,1)
void UniformRandom::setSeed(double seedValue)
{
	unsigned int seed = // convert to integer assuming 32 bits
		static_cast<unsigned int>(seedValue*0xffffffffU);
	setSeed(seed);
}

// Set the seed based on the current time.
// The selection of routines for getting time is
// somewhat limited by sticking to ANSI routines.
// This is not even remotely a good way to generate
// a large number of good seed values.
void UniformRandom::setSeed()
{
	time_t			t;
	clock_t			ct;

	const int		n=8;
	unsigned int	array[n];
	int				k;


	// Get times since they are the only entropy available
	time(&t);				// current time in seconds
	ct=clock();				// clock ticks for this task so far

	// Build an array with starting values
	array[0]=array[3]=++_seedCount;
	array[1]=array[4]=ct;
	array[2]=array[5]=t;
	array[6]=array[7]=0;

	// Expand the array with a simple hash. This is useful
	// for breaking up obvious numerical patterns when successive
	// runs do not differ by much in start time.
	for (k=1;k<=n;k++) {
		array[k%n] ^= array[k-1]<<23 ^ array[k-1]>>9;
		array[k%n] += 0xaaaaaaaaU;
	}
	
	// Pass values from the generated array
	setSeed(8,array);
}

// Reset the seed -- subclass must implement if used
void UniformRandom ::resetSeed()
{
	FatalError("(UniformRandom::resetSeed) Subclass must implement.");
}

// Add a portion of the seed
void UniformRandom::addSeed(unsigned int s)
{
	FatalError("(UniformRandom::addSeed) Subclass must implement.");
}

// Default implementation of advance.
// Some generators may implement this more efficiently
// than the brute force used here.
void UniformRandom::advance(int n) {
	int k;

	for (k=0;k<n;k++) {
		basicValue(); // discard the value
	}
}

// Set the seed of the default generator
void UniformRandom::setDefaultSeed(int seedCount, unsigned int* seedValues)
{
	defaultGen()->setSeed(seedCount,seedValues);
}

// Set the seed of the default generator
void UniformRandom::setDefaultSeed(unsigned int seedValue)
{
	defaultGen()->setSeed(seedValue);
}

// Set the seed of the default generator from date and time
void UniformRandom::setDefaultSeed()
{
	defaultGen()->setSeed();
}

// Return a value using default min=0 and max=1 values
double UniformRandom::value()
{
	return defaultGen()->next();
}

// Return a value with specified min and max values
double UniformRandom::value(double min, double max)
{
	return min+(max-min)*(defaultGen()->next());
}

// Return an integer value in the range of
// 0 to n-1 using the default uniform random generator
int UniformRandom::choice(int n)
{
	return int(n*value());
}

// Return an integer value in the range of
// min to max-1 inclusive using the default
// uniform random generator
int UniformRandom::choice(int min, int max)
{
	return int( min + floor((max-min)*value()) );
}

// Return a random value in the interval (0,1).
double UniformRandom::next()
{
	double u;

	// If generator is not seeded, do it now
	if (!isSeeded() ) {
		setSeed();
	}

	// Look for a value in the proper range
	do {
		// Get a value and check the range
		u = basicValue();

		// Check for a stride
		if (_stride > 1) {
			advance(_stride - 1);
		}

	} while (u<=0.0 || u>=1.0);

	// Return value found
	return u;
}

// Return a random value in (min,max)
double UniformRandom::next(double min, double max)
{
	return min+(max-min)*next();
}



// ====================================================================
// WH_UniformRandom class body
// ====================================================================



// Constructors and destructor
WH_UniformRandom::WH_UniformRandom()
{
	// Initialize prime mods
	_P[0]=30269;
	_P[1]=30307;
	_P[2]=30323;

	// Initialize multipliers
	_A[0]=171;
	_A[1]=172;
	_A[2]=170;

	// Initialize state (seed will override)
	resetSeed();
}

WH_UniformRandom::~WH_UniformRandom() {}

// Copy states to an external array
void WH_UniformRandom::getState(double* states)
{
	int k;

	for (k=0;k<3;k++) {
		states[k]=_X[k];
	}
}

// Set state values from an external array
void WH_UniformRandom::setState(double* states)
{
	int k;

	for (k=0;k<3;k++) {
		_X[k]=static_cast<unsigned int>(states[k]);
	}
}

// Get a random value
double WH_UniformRandom::basicValue()
{
	double u0,u1,u2;

	// Update the state for each generator
	_X[0]=(_A[0]*_X[0])%_P[0];
	_X[1]=(_A[1]*_X[1])%_P[1];
	_X[2]=(_A[2]*_X[2])%_P[2];
	
	// Get individual generator values
	u0=double(_X[0])/_P[0];
	u1=double(_X[1])/_P[1];
	u2=double(_X[2])/_P[2];

	// Combine modulo 1 while trying to somewhat minimize
	// loss of low order bits due to floating point rounding
	u0 = u0+u1<1.0 ? u0+u1 : (u0-0.5)+(u1-0.5);
	u0 = u0+u2<1.0 ? u0+u2 : (u0-0.5)+(u2-0.5);
	
	// Return random value. u==0 or u==1 is handled by caller.	
	return u0;
}

// Clear the seed to a known value
void WH_UniformRandom::resetSeed()
{
	_X[0]=0;
	_X[1]=0;
	_X[2]=0;
}

// Add a seed value to each of the state variables.
// This is basically just a hash with no particular logic.
void WH_UniformRandom::addSeed(unsigned int seed)
{
	const unsigned int two16 = 0x10000;
	int k;

	// Multiply prior values by a constant and add the seed value.
	for (k=0;k<3;k++) {

		_X[k]=(two16*_X[k])%_P[k];
		_X[k]=(two16*_X[k])%_P[k];
		_X[k]=(seed + _X[k])%_P[k];
	}

	// If any of the state values are zero, advance to
	// the next state in sequence until there are no zeros.
	while (_X[0]==0 || _X[1]==0 || _X[2]==0) {
		_X[0]=(_X[0]+1)%_P[0];
		_X[1]=(_X[1]+1)%_P[1];
		_X[2]=(_X[2]+1)%_P[2];
	}
}



// ====================================================================
// LF_UniformRandom class body
// ====================================================================


// Static values
const int LF_UniformRandom::_L = 521;
const int LF_UniformRandom::_K = 32;

// Constructor and destructor
LF_UniformRandom::LF_UniformRandom(UniformRandom* igen)
{
	// Set initial values
	_offset = -1;							// special value for first time
	_initGen = igen;						// save default initialization generator
}

LF_UniformRandom::~LF_UniformRandom()
{
	delete _initGen;
}

// Set the initialization generator
void LF_UniformRandom::initGen(UniformRandom* gen)
{
	// Dispose of the old generator and set a new one
	delete _initGen;
	_initGen = gen;

	// Force filling X with new values on next access
	_offset = -1;
}

// Return the next random value using a (521,32)
// additive lagged Finonacci rule.
double LF_UniformRandom::basicValue()
{
	// Parameters for the irreducible polynomial
	int			i;


	// Advance the current offset into the X array
	// Sequence is X is such that X[n], X[n+1], ...
	// correspond with more recent to least recent.
	if (--_offset<0) {

		// First time (or refill), need to fill X with values
		if (_offset < -1) {
			fillX();
		}

		// Wrap around the table
		_offset = _L-1;
	}

	// Locate entry for lagged Finonacci generator
	if ((i=_offset+_K)>=_L) {
		i -= _L;
	}

	// Apply the ALF rule of adding modulo 1.0.
	_X[_offset] += _X[_offset]+_X[i]<1 ? _X[i] : _X[i]-1.0;

	return _X[_offset];
}

// Fill the state vector X with values
void LF_UniformRandom::fillX()
{
	int i,j;
	double u;

	// Fill the state vector
	for (i=0;i<_L;i++) {
		_X[i] = _initGen->next();
	}

	// Do a random permutation using the initialization generator.
	// This should reduce any obvious correlations between adjacent
	// values, as for example can occur in linear conguent methods.
	for (i=0;i<_L;i++) {
		j = int(_L*_initGen->next());
		u = _X[i];
		_X[i]=_X[j];
		_X[j]=u;
	}
}

// Seed the initialization generator
void LF_UniformRandom::setSeed(int seedCount, unsigned int* seedValues)
{
	_initGen->setSeed(seedCount, seedValues);
	_isSeeded = true;
}

// Seed the initialization generator
void LF_UniformRandom::setSeed(unsigned int seedValue)
{
	_initGen->setSeed(seedValue);
	_isSeeded = true;
}

// Seed the initialization generator
void LF_UniformRandom::setSeed()
{
	_initGen->setSeed();
	_isSeeded = true;
}



// ====================================================================
// MT19937_UniformRandom class body.
// This is only glue code. The bulk of the implementation is in
// the third party source file bnsf_math_3rd_party.cpp.
// ====================================================================



// Constructor and destructor
MT19937_UniformRandom::MT19937_UniformRandom() 
{
	// Force initialization of the state values on first use
	mti = 624+1;
}

MT19937_UniformRandom::~MT19937_UniformRandom() {}


// Set set from array values
void MT19937_UniformRandom::setSeed(int seedCount, unsigned int* seedValues)
{	
	// Must copy seed to long array to match type
	unsigned long* key = new unsigned long[seedCount];
	for (int k=0;k<seedCount;k++) key[k]=seedValues[k];
	
	// Invoke the init function in the package
	init_by_array(key,seedCount); 
	delete key;
	
	// Remember that seeding is done
	_isSeeded = true;
	advance(_offset);
}

// Set seed from an integer
void MT19937_UniformRandom::setSeed(unsigned int seedValue)
{	
	init_genrand(static_cast<unsigned long>(seedValue));
	_isSeeded = true;
	advance(_offset);
}




// ====================================================================
// NonuniformRandom class body
// ====================================================================



// Constructor and destructor
NonuniformRandom::NonuniformRandom(UniformRandom* unif) 
{ 
	_uniformRandom=unif; 
}

NonuniformRandom::~NonuniformRandom() {}

// Set the source of uniform random numbers
void NonuniformRandom::uniformRandom(UniformRandom* unif)
{
	// Make sure change is allowed
	if (_uniformRandom != unif && _uniformRandom!=NULL) {
		FatalError("(NonuniformRandom::uniformRandom) "
			"Cannot change source of uniform random numbers once set");
	}

	// Save the source of uniform random numbers
	_uniformRandom = unif;
}



// ====================================================================
// NormalRandom class body
// ====================================================================



// Constructors and destructor
NormalRandom::NormalRandom(double mean, double sdev, UniformRandom* unif)
: NonuniformRandom(unif)
{
	_mean = mean;
	_sdev = sdev;
	_nX = 0;
	_X[0] = _X[1] = 0;
}

NormalRandom::~NormalRandom() {}

// Return a normally distributed random value 
// This uses the Box-Muller method which generates values
// in N(0,1) from pairs of uniform random values.
double NormalRandom::next()
{
	double u1,u2,r,theta;

	if (_nX == 0) {

		// Get the pair of uniform values
		u1 = runif();
		u2 = runif();

		// Convert to a pair with normal density
		r=sqrt(-2.0*log(u1));
		theta=2*Pi*u2;
		_X[0]=r*sin(theta);
		_X[1]=r*cos(theta);
		_nX = 2;
	}

	// Return a scaled value using the cached N(0,1) value
	return _mean + _sdev * _X[--_nX];
}

// Return a normally distributed random value with the
// mean and standard deviation as supplied using the source
// of random numbers provided without causing any interactions
// among the streams of random numbers.
double NormalRandom::value(double mean, double sdev, UniformRandom* unif)
{
	// This uses the Box-Muller method to generate a single
	// normal value. Caching is not done to avoid interactions
	// among the streams of uniform random numbers.

	double u1,u2,r,theta;
	
	u1 = runif(unif);
	u2 = runif(unif);

	r=sqrt(-2.0*log(u1));
	theta=2*Pi*u2;

	return mean + r*sin(theta)*sdev;
}



// ====================================================================
// ExponentialRandom class body
// ====================================================================



// Constructors and destructor
ExponentialRandom::ExponentialRandom(double mean, UniformRandom* unif)
: NonuniformRandom(unif) 
{
	_mean = mean;
}

ExponentialRandom::~ExponentialRandom() {}


// Return an exponentially distributed random value 
// generated from a uniform random number.
double ExponentialRandom::next()
{
	// Get a uniform random value in (0,1)
	double u = runif();

	// Use the standard transform. In theory the structure
	// of floating point could distort the distribution slightly.
	// Using 1-u helps avoid generating excess large values.
	return -log(1-u)*_mean;
}

// Return a exponentially distributed random value using the
// mean parameter value and the uniform random source given.
double ExponentialRandom::value(double mean, UniformRandom* unif)
{
	// Get a uniform random value in (0,1)
	double u = runif(unif);

	// Use the standard transform. In theory the structure
	// of floating point could distort the distribution slightly.
	// Using 1-u helps avoid generating excess large values.
	return -log(1-u)*mean;
}



// ====================================================================
// PoissonRandom class body
// ====================================================================



PoissonRandom::PoissonRandom(
	double			mean,
	UniformRandom*	unif, 
	MethodType		useMethod) : NonuniformRandom(unif)
{
	// First a sanity check
	if (mean<0) {
		FatalError("(PoissonRandom::PoissonRandom) Mean cannot be negative");
	}

	// Save parameters and derived values
	_mean = mean;
	_logMean = log(_mean);
	_modeDensity = 0;
	_cdfLookup = NULL;
	_generator = NULL;

	// If default method is selected, choose which real method is to be used.
	_method = useMethod;
	if (useMethod == defaultMethod) {
		if (mean<=10) {
			_method = inverseCDF;
		}
		else if (mean<=100) {
			_method = ROUMethod;
		}
		else {
			_method = normalApprox;
		}
	}

	// Setup for the method selected
	switch (_method) {

	case inverseCDF:
		_cdfLookup = new InverseCDF(this,int(1+mean+6*sqrt(mean)));
		break;

	case ROUMethod:
		_generator = new RatioOfUniforms(this);
		break;

	case normalApprox:
 		_method = normalApprox;
		break;

	default:
		FatalError("(PoissonRandom::PoissonRandom) Invalid method param.");	
	}
}

PoissonRandom::~PoissonRandom() 
{
	// Delete allocated objects
	delete _cdfLookup;
	delete _generator;
}

// Generate a Poisson distributed random variable using either
// an inverse CDF lookup or a generalized ratio-of-uniforms method.
int PoissonRandom::inext()
{
	if (_method==inverseCDF) {
		// Generate using an inverse CDF lookup
		// Limit table size to twice the mean, which
		// should be more than enough to almost always
		// find and entry in the table.

		int k = _cdfLookup->inext();
		int ix;
		double expMean,expDen,c,u;

		if (k>=0) {
			return k;
		}
		else {
	
		// When a value beyond the table is generated, an
		// accept-reject scheme is used for the right
		// tail values with a discrete exponential distribution
		// as the covering distribution for the tail.

			k = _cdfLookup->cdfSize()-1;	
			expMean = 1/log(k/_mean);
			c = density(k);
			for (;;) {
				ix = int(ExponentialRandom::value(expMean, uniformRandom() ));
				expDen = exp(-ix/expMean);
				u = runif();
				if (c*u < density(k+ix)/expDen) {
					return k+ix;
				}
			}
		}
	}

	else if (_method==ROUMethod) {
		// Use a Ratio-Of-Uniforms Function to generate a number
		return _generator->inext();
	}

	else if (_method==normalApprox) {
		// Use a normal approximation -- of limited accuracy for small means.
		int k;
		NormalRandom rnorm(_mean,sqrt(_mean));
		do {
			k = int(rnorm.next()+0.5); // round the continuous value
		} while (k<0);
		return k;
	}

	else {
		FatalError("(PoissonRandom::inext) Invalid method found");
		return 0;
	}
}

// Return the mode value as an integer
int PoissonRandom::imode()
{
	// Round down the mean for the mode.
	return int(_mean);
}

// Return the density at the mode
double PoissonRandom::modeDensity()
{
	if (_modeDensity == 0.0) {
		// Compute density on first reference
		_modeDensity = density(imode() );
	}
	return _modeDensity;
}

// Compute the density at an integer value
double PoissonRandom::density(int k)
{
	double	den;
	int		n;

	// Compute by brute force for small k
	if (k<10) {
		if (k<0) {
			den = 0;
		}
		else {
			den = exp(- _mean);
			for (n=1;n<=k;n++) {
				den *= _mean / n;
			}
		}
	}
	else {
		// Compute density while avoiding overflow
		den = exp(-_mean + k*_logMean - logFactorial(k));
	}
	return den;
}

// Return a Poisson distributed random value given the mean
int PoissonRandom::ivalue(double mean, UniformRandom* unif)
{
	// For a small mean, use a direct method.
	// Otherwise, use a one-time instance.

	if (mean<16) {

		double	expMean = exp(-mean);
		double	x=runif(unif);
		int		ivalue=0;

		while (x>expMean) {
			ivalue++;
			x *= runif(unif);
		}
		return ivalue;
	}
	else {		
		PoissonRandom prand(mean,unif,ROUMethod);
		return prand.inext();
	}
}



// ====================================================================
// BinomialRandom class body
// ====================================================================



// Constructe a new instance with parameters provided
BinomialRandom::BinomialRandom(
	int				size, 
	double			prob,
   UniformRandom*	unif, 
   MethodType		useMethod) : NonuniformRandom(unif)
{
	// First some sanity checks
	if (size<0) {
		FatalError("(BinomialRandom::BinomialRandom) Size is negative");
	}
	if (prob<0 || prob>1) {
		FatalError("(BinomialRandom::BinomialRandom) Probability not in range");
	}

	// Initialize things
	_N = size;
	_prob = prob;
	_logP = log(prob);
	_logQ = log(1-prob);
	_generator = NULL;
	_mode = -1;
	_modeDensity = 0;

	// If default is specified, select the method to use
	_method = useMethod;
	if (useMethod == defaultMethod) {
		if (_N <= 50) {
			_method = inverseCDF;
		}
		else {
			_method = ROUMethod;
		}
	}

	// Set up for the method selected
	switch (_method) {

	case inverseCDF:
		_generator = new InverseCDF(this, _N);
		break;

	case ROUMethod:
		_generator = new RatioOfUniforms(this);
		break;

	case normalApprox:
 		_method = normalApprox;
		break;

	default:
		FatalError("(BinomialRandom::BinomialRandom) Invalid method param.");	
	}
}

// Destroy this instance
BinomialRandom::~BinomialRandom()
{
	// Clean up any allocated objects
	delete _generator;
}

// Compute the probability density at at value
double BinomialRandom::density(int k)
{
	double		den;
	int			i;

	// For out of range values, return 0
	if (k<0 || k>_N) {
		return 0;
	}

	// Check for extreme probability values
	if (_prob==0 || _prob ==1) {
		if ( (_prob==0 && k==0) || (_prob==1 && k==_N))
			return 1;
		else
			return 0;
	}

	// For small k, try to compute exactly via brute force.
	if (k<=10) {
		den=comb(_N,k);
		for (i=1;i<=k;i++) {
			den *= _prob;
		}
		for (i=1;i<=_N-k;i++) {
			den *= (1-_prob);
		}
		
	}
	else {
		// For larger values use logs to avoid overflows
		den = exp(logComb(_N,k)+k*_logP+(_N-k)*_logQ);
	}

	return den;
}

// Return the density at an adjacent point using an already
// computed density. dir is either +1 or -1. 
// The density at k+dir is returned.
double BinomialRandom::adjDensity(int k, double denAtK, int dir)
{
	// Look at cases with special care for extreme values
	switch (dir) {
	case +1:
		if (k<_N) {
			if (_prob<1)
				return denAtK*(_N-k)/(k+1)*_prob/(1-_prob);
			else if (k==_N-1)
				return 1;
			else
				return 0;
		}
		else {
			return 0;
		}

	case -1:
		if (k>0) {
			if (_prob>0)
				return denAtK*k/(_N-k+1)*(1-_prob)/_prob;
			else if (k==1)
				return 1;
			else
				return 0;
		}
		else {
			return 0;
		}

	default:
		FatalError("(BinomialRandom::adjDensity) Invalid dir value.");
		return 0;
	}
}

// Get a mode value and save it for later
int BinomialRandom::imode()
{
	double		mean,tden;

	mean = _N * _prob;
	_mode = int(mean);
	_modeDensity = density(_mode);

	// Unless mean is exactly an integer, need to look
	// for a mode at integer values above and below the mean.
	if (_mode != mean) {
		tden=adjDensity(_mode,_modeDensity,+1); // density at _mode+1
		if (tden>_modeDensity) {
			_mode++;
			_modeDensity = tden;
		}
	}
	return _mode;
}

// Get the density at the mode
double BinomialRandom::modeDensity()
{
	// First time, let imode get the density
	if (_mode<0) {
		imode();
	}
	return _modeDensity;
}

// Return a random value as an integer
int BinomialRandom::inext()
{
	int		k;

	// See which method to use
	switch( _method) {

	case ROUMethod:
	case inverseCDF:
		// When doing a look up in the CDF table. There is a small 
		// probability that round-off error can cause a cdf look up
		// failure, in which case we try again.
		do {
			k = _generator->inext();
		} while (k<0 || k>_N);
		break;
		
	default:
		// Use a normal approximation
		NormalRandom norm(_N*_prob,sqrt(_N*_prob*(1-_prob)));
		do {
			k = int(norm.next()+0.5); // round the continuous value
		} while (k<0 || k>_N);
		break;

	}
	return k;
}


// Return a binomial random value using parameters
int BinomialRandom::ivalue(int size, double prob, UniformRandom* unif)
{
	// For small size values use a direct method.
	// Otherwise, use a one-time instance.
	if (size<=16) {

		int k, n=0;

		for (k=0; k<size; k++) {
			if (runif(unif)<prob) n++;
		}
		return n;
	}
	else {
		BinomialRandom brand(size, prob, unif, ROUMethod);
		return brand.inext();
	}
}
