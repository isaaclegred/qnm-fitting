#pragma once
/*NR3BEG nr3testutils*/
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif

#ifndef _NR_RAN
#include "ran.h"
#define _NR_RAN
#endif

Ran ran(10101); // global ranno generator

void ranmat(MatDoub &a, Doub dadd=0.) { // fill matrix with ran
	Int m = a.nrows(), n=a.ncols();
	for (Int i=0;i<m;i++) for (Int j=0;j<n;j++) a[i][j] = 2.*ran.doub() - 1.;
	for (Int i=0;i<MIN(m,n);i++) a[i][i] += dadd;
}

void ranvec(VecDoub &a) { // fill vector with ran
	Int m = a.size();
	for (Int i=0;i<m;i++) a[i] = 2.*ran.doub() - 1.;
}

template <class T>
NRmatrix<T> matmul(const NRmatrix<T> &a, const NRmatrix<T> &b) {
	Int i,j,k,m=a.nrows(), n=b.ncols(), p=a.ncols();
	if (p != b.nrows()) throw("impossible matrix multiply");
	Doub sum;
	NRmatrix<T> c(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		sum = 0.;
		for (k=0;k<p;k++) sum += a[i][k]*b[k][j];
		c[i][j] = sum;
	}
	return c;
}

template <class T>
NRmatrix<T> transpose(const NRmatrix<T> &a) {
	Int i,j,m=a.nrows(), n=a.ncols();
	NRmatrix<T> c(n,m);
	for (i=0;i<m;i++) for (j=0;j<n;j++) { c[j][i] = a[i][j]; }
	return c;
}

template <class T>
NRvector<T> matmul(const NRmatrix<T> &a, const NRvector<T> &b) {
	Int i,k,m=a.nrows();
	Doub sum;
	NRvector<T> c(m);
	for (i=0;i<m;i++) {
		sum = 0.;
		for (k=0;k<m;k++) sum += a[i][k]*b[k];
		c[i] = sum;
	}
	return c;
}

template <class T>
T maxel(const NRmatrix<T> &a) {
	Int i,j,m=a.nrows(), n=a.ncols();
	Doub max = 0.;
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		if (abs(a[i][j]) > max) max = abs(a[i][j]);
	}
	return T(max);
}

template <class T>
T maxel(const NRvector<T> &a) {
	Int i,m=a.size();
	Doub max = 0.;
	for (i=0;i<m;i++) {
		if (abs(a[i]) > max) max = abs(a[i]);
	}
	return T(max);
}

template <class T>
NRmatrix<T> ident(int n, T dum) {
	NRmatrix<T> c;
	c.assign(n,n,(T)0);
	for (Int i=0;i<n;i++) c[i][i] = (T)1;
	return c;
}

template <class T>
NRmatrix<T> matsub(const NRmatrix<T> &a, const NRmatrix<T> &b) {
	Int i,j,m=a.nrows(), n=a.ncols();
	if (a.nrows() != b.nrows() || a.ncols() != b.ncols()) throw("bad matsub");
	NRmatrix<T> c(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		c[i][j] = a[i][j]-b[i][j];
	}
	return c;
}

template <class T>
NRvector<T> vecsub(const NRvector<T> &a, const NRvector<T> &b) {
	Int i,m=a.size();
	if (a.size() != b.size()) throw("bad vecsub");
	NRvector<T> c(m);
	for (i=0;i<m;i++) {
		c[i] = a[i]-b[i];
	}
	return c;
}

template <class T>
NRvector<T> vecadd(const NRvector<T> &a, const NRvector<T> &b) {
	Int i,m=a.size();
	if (a.size() != b.size()) throw("bad vecadd");
	NRvector<T> c(m);
	for (i=0;i<m;i++) {
		c[i] = a[i]+b[i];
	}
	return c;
}

/*NR3END nr3testutils*/
