#include "nr3devel.h"
#include "nr3testutils.h"
#include <vector>

#include "test_Fitsep.hpp"



int main(int argc,char *argv[])
{
  struct fgauss{
    void operator()(const Doub x, VecDoub_I &a, VecDoub &u, MatDoub_O &duda){
      return fgauss_sep(x, a, u, duda);}
    void operator()(const Doub t, VecDoub_I & params, Doub & f_eval, VecDoub_O &jac){}
  };
  bool debug=false;
  Int i,j,N=100,MA=4,mfit=MA,P=MA/2,pfit=P;
  VecDoub a{2.0,3.0,5.0,3.0};
  VecDoub c{5.0,2.0};
  //	Doub gguess[]={4.5,2.2,2.8,2.5,4.9,2.8};
  VecDoub guess{2.2,2.8,4.9,2.8};
  Doub SPREAD=0.01;
  VecDoub x(N),y(N),sig(N);
  bool localflag, globalflag=false;

  if ((argc > 1) && (argv[1] == "debug")) debug=true;

  // Test Sep_Marquardt
  cout << "Testing Sep_Marquardt" << endl;

  Normaldev ndev(0.0,1.0,17);
  // First try a sum of two Gaussians
  for (i=0;i<N;i++) {
    x[i]=0.1*(i+1);
    y[i]=0.0;
    for (j=0;j<MA;j+=2)
      y[i] += c[j/2]*exp(-SQR((x[i]-a[j])/a[j+1]));
    y[i] *= (1.0+SPREAD*ndev.dev());
    sig[i]=SPREAD*y[i];
  }
  Sep_marquardt<fgauss> myfit(x,y,sig,guess,fgauss(),P);
  //Sep_marquardt myfit(x,y,sig,guess,fgauss_sep,P,1.e-6);
  myfit.fit();

  cout << setw(18) << "chi-squared:" << setw(13) << myfit.chisq << endl;
  cout << fixed << setprecision(6);
  for (i=0;i<P;i++) cout << setw(9) << myfit.c[i];
  for (i=0;i<MA;i++) cout << setw(9) << myfit.a[i];
  cout << endl;
  cout << "Uncertainties:" << endl;
  for (i=0;i<P;i++) cout << setw(9) << sqrt(myfit.covar[i][i]);
  for (i=0;i<MA;i++) cout << setw(9) << sqrt(myfit.covar[P+i][P+i]);
  cout << endl;
  cout << "Expected results:" << endl;
  for (i=0;i<P;i++) cout << setw(9) << c[i];
  for (i=0;i<MA;i++) cout << setw(9) << a[i];
  cout << endl << endl;

  for (j=0;j<P;j++) {
    localflag = fabs(myfit.c[j]-c[j]) > 2.0*sqrt(myfit.covar[j][j]);
    globalflag = globalflag || localflag;
    if (localflag) {
      cout << "*** Sep_Marquardt: Linear fitted parameters not within estimated uncertainty" << endl;
      if (debug) throw("");
    }
  }

  for (j=0;j<MA;j++) {
    localflag = fabs(myfit.a[j]-a[j]) > 2.0*sqrt(myfit.covar[P+j][P+j]);
    globalflag = globalflag || localflag;
    if (localflag) {
      cout << "*** Sep_Marquardt: Nonlinear fitted parameters not within estimated uncertainty" << endl;
      if (debug) throw("");
    }
  }

  // Test the hold() method on linear parameter 0 and nonlinear parameter 2
  myfit.holdc(0,5.0);
  myfit.hold(2,5.0);
  myfit.fit();

  cout << setw(18) << "chi-squared:" << setw(13) << myfit.chisq << endl;
  for (i=0;i<P;i++) cout << setw(9) << myfit.c[i];
  for (i=0;i<MA;i++) cout << setw(9) << myfit.a[i];
  cout << endl;
  cout << "Uncertainties:" << endl;
  for (i=0;i<P;i++) cout << setw(9) << sqrt(myfit.covar[i][i]);
  for (i=0;i<MA;i++) cout << setw(9) << sqrt(myfit.covar[P+i][P+i]);
  cout << endl;
  cout << "Expected results:" << endl;
  for (i=0;i<P;i++) cout << setw(9) << c[i];
  for (i=0;i<MA;i++) cout << setw(9) << a[i];
  cout << endl << endl;

  localflag = (myfit.c[0] != c[0]);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Sep_Marquardt: A held parameter does not have its assigned value" << endl;
    if (debug) throw("");
  }

  localflag = (myfit.a[2] != a[2]);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Sep_Marquardt: A held parameter does not have its assigned value" << endl;
    if (debug) throw("");
  }

  localflag = (myfit.covar[0][0] != 0.0);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Sep_Marquardt: A held parameter does not have uncertainty=0.0" << endl;
    if (debug) throw("");
  }

  localflag = (myfit.covar[P+2][P+2] != 0.0);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Sep_Marquardt: A held parameter does not have uncertainty=0.0" << endl;
    if (debug) throw("");
  }

  for (j=0;j<P;j++) {
    localflag = fabs(myfit.c[j]-c[j]) > 2.0*sqrt(myfit.covar[j][j]);
    globalflag = globalflag || localflag;
    if (localflag) {
      cout << "*** Sep_Marquardt: Fitted parameters (with 2 parameters held) not within estimated uncertainty" << endl;
      if (debug) throw("");
    }
  }
  for (j=0;j<MA;j++) {
    localflag = fabs(myfit.a[j]-a[j]) > 2.0*sqrt(myfit.covar[P+j][P+j]);
    globalflag = globalflag || localflag;
    if (localflag) {
      cout << "*** Sep_Marquardt: Fitted parameters (with 2 parameters held) not within estimated uncertainty" << endl;
      if (debug) throw("");
    }
  }

  localflag=false;
  for (i=0;i<P+MA;i++) {
    for (j=0;j<P+MA;j++) {
      if (i==0 || i==P+2 || j==0 || j==P+2)
        localflag = localflag || myfit.covar[i][j] != 0.0;
      else
        localflag = localflag || myfit.covar[i][j] == 0.0;
    }
  }
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Sep_Marquardt: Covariance matrix with 2 held parameters has incorrect pattern" << endl;
    if (debug) throw("");
  }

  // Test the free() method
  myfit.freec(0);
  myfit.fit();

  //	cout << setw(18) << "chi-squared:" << setw(13) << myfit.chisq << endl;
  //	for (i=0;i<MA;i++) cout << setw(9) << myfit.a[i];
  //	cout << endl;
  //	cout << "Uncertainties:" << endl;
  //	for (i=0;i<MA;i++) cout << setw(9) << sqrt(myfit.covar[i][i]);
  //	cout << endl;
  //	cout << "Expected results:" << endl;
  //	for (i=0;i<MA;i++) cout << setw(9) << a[i];
  //	cout << endl << endl;

  localflag = (myfit.c[0] == c[0]);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Fitlin: A freed parameter still has its assigned value" << endl;
    if (debug) throw("");
  }

  localflag = (myfit.a[2] != a[2]);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Fitlin: A held parameter does not have its assigned value" << endl;
    if (debug) throw("");
  }

  localflag = (myfit.covar[0][0] == 0.0);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Fitlin: A freed parameter still has uncertainty=0.0" << endl;
    if (debug) throw("");
  }

  localflag = (myfit.covar[P+2][P+2] != 0.0);
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Fitlin: A held parameter does not have uncertainty=0.0" << endl;
    if (debug) throw("");
  }

  for (j=0;j<P;j++) {
    localflag = fabs(myfit.c[j]-c[j]) > 2.0*sqrt(myfit.covar[j][j]);
    globalflag = globalflag || localflag;
    if (localflag) {
      cout << "*** Fitlin: Fitted parameter (with 1 parameters held) not within estimated uncertainty" << endl;
      if (debug) throw("");
    }
  }
  for (j=0;j<MA;j++) {
    localflag = fabs(myfit.a[j]-a[j]) > 2.0*sqrt(myfit.covar[P+j][P+j]);
    globalflag = globalflag || localflag;
    if (localflag) {
      cout << "*** Fitlin: Fitted parameter (with 1 parameters held) not within estimated uncertainty" << endl;
      if (debug) throw("");
    }
  }

  localflag=false;
  for (i=0;i<P+MA;i++) {
    for (j=0;j<P+MA;j++) {
      if (i==P+2 || j==P+2)
        localflag = localflag || myfit.covar[i][j] != 0.0;
      else
        localflag = localflag || myfit.covar[i][j] == 0.0;
    }
  }
  globalflag = globalflag || localflag;
  if (localflag) {
    cout << "*** Fitlin: Covariance matrix with 1 held parameters has incorrect pattern" << endl;
    if (debug) throw("");
  }

  if (globalflag) cout << "Failed" << endl << endl;
  else cout << "Passed" << endl << endl;
}
