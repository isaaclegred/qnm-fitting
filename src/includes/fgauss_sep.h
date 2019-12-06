#pragma once
/*NR3BEG fgauss_sep*/
void fgauss_sep(const Doub x, VecDoub_I &a, VecDoub &u, MatDoub_O &duda) {
/*290034*/
	Int i,j,na=a.size(),m=u.size();
	Doub ex,arg;
	for (j=0;j<m;j++)
		for (i=0;i<na;i++)
			duda[j][i]=0.0;
	for (i=0;i<na-1;i+=2) {
		arg=(x-a[i])/a[i+1];
		ex=exp(-SQR(arg));
		u[i/2]=ex;
		duda[i/2][i]=2.0*ex*arg/a[i+1];
		duda[i/2][i+1]=2.0*ex*SQR(arg)/a[i+1];
	}
}
/*NR3END fgauss_sep*/
