#pragma once
#include "nr3devel.h"
#include "cholesky.h"
/*NR3BEG marquardt*/
template<typename T>
struct Marquardt {
/*240001*/
  Int n, m, mfit;
  /*540020 1 2 2 2*/
  VecDoub_I &x,&y,&sig;
  T funcs;
  Doub tol;
  VecBool ia;
  VecDoub a;
  /*540012 1 2 2 2*/
  MatDoub alpha, covar, tempmat;
  Doub chisq;
  Int iter;
  Doub cosalpha;
  Bool done;
  Marquardt(VecDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig, VecDoub_I &aa,
            T funks, const Doub TOL=1.e-3) : n(xx.size()), m(aa.size()), x(xx), y(yy), sig(ssig),
                         funcs(funks), tol(TOL), ia(m,true), a(aa), covar(m,m), cosalpha(1.0) {}
  /*340013*/

  void hold(const Int i, const Doub val) {ia[i]=false; a[i]=val;}
  void free(const Int i) {ia[i]=true;}
  /*340014*/

  void fit() {
    /*340002*/
    const Int ITMAX=200;
    const Doub MINLAMBDA=100.0*numeric_limits<Doub>::epsilon(),
      TOLX=1.0e-8,TOLCHISQ=1.0e-9, LAMBDAFACT=1.0e-3;
    /*440021*/
    Doub ochisq,alphamax,lambda,normda,norma,dl,df,nu=2.0;
    Bool sing;
    done=false;
    mfit=0;
    /*540022 4 6 5 5*/
    for (Int j=0;j<m;j++) if (ia[j]) mfit++;
    VecDoub beta(mfit),betanew(mfit),da(mfit),anew(a);
    MatDoub alphanew(mfit,mfit);
    alpha.resize(mfit,mfit);
    mrqcof(a,alpha,beta);
    /*540003 4 5 5 5*/
    ochisq=chisq;
    alphamax=0.0;
    for (Int i=0;i<mfit;i++){
      /*540023 4 5 5 5*/
      if (alpha[i][i] > alphamax) alphamax=alpha[i][i];}
    lambda=LAMBDAFACT*alphamax;
    for (Int its=0;its<ITMAX;its++) {
      iter=its;
      alphamax=0.0;
      for (Int i=0;i<mfit;i++)
        if (alpha[i][i] > alphamax) alphamax=alpha[i][i];
      lambda=MAX(lambda,MINLAMBDA*alphamax);
      /*540024 6 7 7 7*/
      sing=true;
      Cholesky *chol;
      while (sing) {
        MatDoub aa(alpha);
        /*540025 4 5 5 5*/
        for (Int i=0;i<mfit;i++)
          aa[i][i] += lambda;
        chol=new Cholesky(aa);
        /*540026 4 5 5 5*/
        sing=chol->sing;
        if (sing) {
          /*540027 4 5 5 5*/
          lambda *=10.0;
          delete chol;
        }
      }
      chol->solve(beta,da);
      /*540005 4 5 5 5*/
      delete chol;
      normda=norma=0.0;
      /*540028 4 5 5 5*/
      for (Int j=0,i=0;i<m;i++) {
        if (ia[i]) {
          norma += SQR(a[i]);
          normda += SQR(da[j++]);
        }
      }
      normda=sqrt(normda);
      norma=sqrt(norma);
      if (normda <= TOLX*(norma+TOLX)) {
        done=true;
        break;
      }
      for (Int i=0,j=0;i<m;i++) {
        /*540029 4 5 5 5*/
        if (ia[i])
          anew[i]=a[i]+da[j++];
      }
      dl=0.0;
      /*540030 4 5 5 5*/
      for (Int i=0;i<mfit;i++) {
        dl += da[i]*(lambda*da[i]+beta[i]);
      }
      mrqcof(anew,alphanew,betanew);
      /*540031 4 5 5 5*/
      df=ochisq-chisq;
      if (abs(chisq) < TOLCHISQ)
        done=true;
      cosalpha=sqrt(abs(dl/chisq));
      if (!done && cosalpha < tol)
        /*540032 4 5 5 5*/
        done=true;
      if (done) {
        /*540033 4 5 5 5*/
        if (df > 0.0) {
          for (Int j=0;j<mfit;j++)
            for (Int k=j;k<mfit;k++) alpha[j][k]=alphanew[j][k];
          for (Int j=0;j<m;j++) a[j]=anew[j];
        } else
          chisq=ochisq;
        break;
      }
      if (dl > 0.0 && df > 0.0) {
        /*540034 4 5 5 5*/
        ochisq=chisq;
        for (Int j=0;j<mfit;j++) {
          for (Int k=j;k<mfit;k++) alpha[j][k]=alphanew[j][k];
          beta[j]=betanew[j];
        }
        for (Int j=0;j<m;j++) a[j]=anew[j];
        lambda *= MAX(1.0/3.0,1.0-pow(2.0*df/dl-1.0,3));
        /*54003510111111*/
        nu=2.0;
      } else {
        /*540009 4 5 5 5*/
        lambda *= nu;
        nu *= 2.0;
      }
    }
    if (iter >= ITMAX) throw("too many iterations in marquardt");
    compute_covar();
    /*540036 4 5 5 5*/
  }

  virtual void mrqcof(VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta) {
    /*339901*/
    Int i,j,k,l,mm;
    Doub ymod,wt,sig2i,dy;
    VecDoub dyda(m);
    for (j=0;j<mfit;j++) {
      /*539902 2 3 3 3*/
      for (k=j;k<mfit;k++) alpha[j][k]=0.0;
      beta[j]=0.0;
    }
    chisq=0.;
    for (i=0;i<n;i++) {
      /*539903 2 3 3 3*/
      funcs(x[i],a,ymod,dyda);
      sig2i=1.0/(sig[i]*sig[i]);
      dy=y[i]-ymod;
      for (j=0,l=0;l<m;l++) {
        if (ia[l]) {
          wt=dyda[l]*sig2i;
          for (k=j,mm=l;mm<m;mm++)
            if (ia[mm]) alpha[j][k++] += wt*dyda[mm];
          beta[j++] += dy*wt;
        }
      }
      chisq += dy*dy*sig2i;
      /*539904 2 3 3 3*/
    }
  }

  virtual void compute_covar() {
    /*332801*/
    Cholesky chol(alpha);
    /*532802 2 3 3 3*/
    if (chol.sing) throw("singular covariance matrix in marquardt");
    chol.inverse(tempmat);
    /*532803 2 3 3 3*/
    for (Int i=0;i<m;i++)
      /*532804 2 3 3 3*/
      for (Int j=0;j<m;j++)
        covar[i][j]=0.0;
    for (Int k=0,i=0;i<m;i++)
      if (ia[i]) {
        for (Int l=0,j=0;j<m;j++)
          if (ia[j]) covar[i][j]=tempmat[k][l++];
        k++;
      }
  }

};
/*NR3END marquardt*/
