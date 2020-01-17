#pragma once
#include <utility>
#include "nr3devel.h"
#include "svd.h"
#include "marquardt.hpp"
#include "vectorhelpers.hpp"
/*NR3BEG sep_marquardt*/
// This class is now templated on a functor `T`, which must overload the function `void operator()`
// with arguments `(const Doub t, VecDoub_I & params, Doub & f_eval, VecDoub_O &jac)`, and with
// `(const Doub  t, VecDoub_I & params , VecDoub_O & f_vals, MatDoub_O & jac)` because `Marquardt`
// requires the first signature, and Sep_marquardt requires the second.  This is not ideal, and
// hopefully can be improved in the future.  Note that the respective `Marquardt` signature is
// not actually called in the file, so the body of the function can be anything, including {}
template <typename T>
struct Sep_marquardt : Marquardt<T> {
/*290001*/
  Int p, pfit;
  /*590002 1 2 2 2*/
  T funcs;
  VecBool ic;
  VecDoub c;
  using Marquardt<T>::n; using Marquardt<T>::m; using Marquardt<T>::mfit;
  using Marquardt<T>::x; using Marquardt<T>::y; using Marquardt<T>::sig;
  using Marquardt<T>::tol;
  using Marquardt<T>::ia; using Marquardt<T>::a;
  using Marquardt<T>::alpha; using Marquardt<T>::covar; using Marquardt<T>::tempmat;
  using Marquardt<T>::chisq;
  using Marquardt<T>::iter;
  using Marquardt<T>::cosalpha;
  using Marquardt<T>::done;

  /*390003*/

  // Alternative version of the function, `T` is a functor whose `operator()` has the same signature
  // as `funks` above
  Sep_marquardt(VecDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig, VecDoub_I &aa,
                T funks, const Int pp,
                const Doub TOL=1.e-3) : Marquardt<T>(xx,yy,ssig,aa,funks,TOL), p(pp),
                                        funcs(funks), ic(p,true), c(p) {
    /*390004*/
    covar.resize(p + m, p + m);
  }

  void holdc(const Int i, const Doub val) {ic[i]=false; c[i]=val;}
  void freec(const Int i) {ic[i]=true;}
  /*390005*/

  void mrqcof(VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta) {
    /*390006*/
    const Doub SVDTOL=1.e-12;
    /*590007 4 5 5 5*/
    Int i,j,k,l;
    Doub tmp,sum,dy;
    pfit=0;
    /*590008 4 5 5 5*/
    for (j=0;j<p;j++) if (ic[j]) pfit++;
    VecDoub b(n),u(p),cfit(pfit);
		MatDoub au(n,pfit),duda(p,m);
		for (i=0;i<n;i++) {
                  /*590009 4 5 5 5*/
                  funcs(x[i],a,u,duda);
                  /*590010 4 5 5 5*/
                  tmp=1.0/sig[i];
                  //std::cout << "sig["<< i << "] is" << sig[i] << "\n";
                  sum=0.0;
                  for (l=0,j=0;j<p;j++) {
                    if (ic[j])
                      /*590011 4 5 5 5*/
                      au[i][l++]=u[j]*tmp;
                    else
                      /*590012 4 5 5 5*/
                      sum += c[j]*u[j];
                  }
                  b[i]=(y[i]-sum)*tmp;
                  //std::cout << "b[" << i << "] is"  << b[i] << "\n";
		}
		SVD svd(au);
                /*590013 4 5 5 5*/
		Doub thresh=SVDTOL*svd.w[0];
		svd.solve(b,cfit,thresh);
                /*590014 4 5 5 5*/
		for (l=0,j=0;j<p;j++){
                  if (ic[j]) c[j]=cfit[l++];}
		for (j=0;j<mfit;j++) {
                  /*590015 4 5 5 5*/
                  for (k=j;k<mfit;k++)
                    {alpha[j][k]=0.0;}
                  beta[j]=0.0;
		}
		chisq=0.0;
		VecDoub ac(mfit);
                /*590020 4 5 5 5*/
		MatDoub uac(pfit,mfit,0.0);
                /*590021 4 5 5 5*/
		if (done) {
                  /*590016 4 6 6 6*/
                  tempmat.resize(pfit+mfit,pfit+mfit);
                  for (i=0;i<mfit+pfit;i++)
                    for (j=i;j<mfit+pfit;j++)
                      tempmat[i][j]=0.0;
		}
		for (i=0;i<n;i++) {
                  /*590017 4 5 5 5*/
                  funcs(x[i],a,u,duda);
                  tmp=1.0/sig[i];
                  sum=0.0;
                  for (l=0;l<p;l++)
                    sum += c[l]*u[l];
                  dy=(y[i]-sum)*tmp;
                  /*590018 4 5 5 5*/
                  chisq += dy*dy;
                  /*590019 4 5 5 5*/
                  for (l=0,k=0;k<m;k++) {
                    if (ia[k]) {
                      sum=0.0;
                      /*590022 4 5 5 5*/
                      for (j=0;j<p;j++) {
                        sum += c[j]*duda[j][k];
                      }
                      sum *= tmp;
                      ac[l]=sum;
                      beta[l] += dy*sum;
                      /*590023 4 5 5 5*/
                      for (j=0;j<pfit;j++)
                        /*590024 4 5 5 5*/
                        uac[j][l] += svd.u[i][j]*sum;
                      l++;
                    }
                  }
                  for (j=0;j<mfit;j++)
                    /*590025 4 5 5 5*/
                    for (k=j;k<mfit;k++)
                      alpha[j][k] += ac[k]*ac[j];
                  if (done) {
                    /*590026 4 5 5 5*/
                    for (j=0;j<pfit;j++)
                      for (k=j;k<pfit;k++)
                        tempmat[j][k] += au[i][j]*au[i][k];
                    for (j=0;j<pfit;j++)
                      for (k=0;k<mfit;k++)
                        tempmat[j][pfit+k] += au[i][j]*ac[k];
                  }
		}
		if (done) {
                  /*590027 4 5 5 5*/
                  for (j=0;j<mfit;j++)
                    for (k=j;k<mfit;k++)
                      tempmat[pfit+j][pfit+k]=alpha[j][k];
		}
		for (j=0;j<mfit;j++) {
                  /*590028 4 5 5 5*/
                  for (k=j;k<mfit;k++) {
                    sum=0.0;
                    for (l=0;l<pfit;l++)
                      sum += uac[l][j]*uac[l][k];
                    alpha[j][k] -= sum;
                  }
		}
  }

  void compute_covar() {
    /*390029*/
    VecDoub tempvec(mfit);
    /*590030 4 5 5 5*/
    mrqcof(a,alpha,tempvec);
    /*590031 4 5 5 5*/
    Cholesky chol(tempmat);
    /*590032 4 5 5 5*/
    if (chol.sing) throw("singular covariance matrix in marquardt");
    chol.inverse(tempmat);
    for (Int i=0;i<m+p;i++)
      /*590033 4 5 5 5*/
      for (Int j=0;j<m+p;j++)
        covar[i][j]=0.0;
    for (Int k=0,i=0;i<p;i++)
      if (ic[i]) {
        for (Int l=0,j=0;j<p;j++)
          if (ic[j]) covar[i][j]=tempmat[k][l++];
        for (Int l=0,j=0;j<m;j++)
          if (ia[j]) covar[i][p+j]=tempmat[k][pfit+l++];
        k++;
      }
    for (Int k=0,i=0;i<m;i++)
      if (ia[i]) {
        for (Int l=0,j=0;j<p;j++)
          if (ic[j]) covar[p+i][j]=tempmat[pfit+k][l++];
        for (Int l=0,j=0;j<m;j++)
          if (ia[j]) covar[p+i][p+j]=tempmat[pfit+k][pfit+l++];
        k++;
      }
  }
};
/*NR3END sep_marquardt*/
