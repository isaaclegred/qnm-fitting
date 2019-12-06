#pragma once
/*NR3BEG svd*/
struct SVD {
/*247616*/
	Int m,n;
	MatDoub u,v;
/*547617 4 5 5 5*/
	VecDoub w;
/*547618 4 5 5 5*/
	Doub eps, tsh;
	SVD(MatDoub_I &a) : m(a.nrows()), n(a.ncols()), u(a), v(n,n), w(n) {
/*347619*/
		eps = numeric_limits<Doub>::epsilon();
		decompose();
		reorder();
		tsh = 0.5*sqrt(m+n+1.)*w[0]*eps;
/*547631 4 5 5 5*/
	}

	void solve(VecDoub_I &b, VecDoub_O &x, Doub thresh);
	void solve(MatDoub_I &b, MatDoub_O &x, Doub thresh);
/*347620*/

	Int rank(Doub thresh);
/*547621 4 5 5 5*/
	Int nullity(Doub thresh);
	MatDoub range(Doub thresh);
	MatDoub nullspace(Doub thresh);

	Doub inv_condition() {
/*547632 410 5 5*/
		return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
	}

	void decompose();
/*547622 4 5 5 5*/
	void reorder();
	Doub pythag(const Doub a, const Doub b);
};
/*NR3END svd*/
/*NR3BEG svd_solve*/
void SVD::solve(VecDoub_I &b, VecDoub_O &x, Doub thresh = -1.) {
/*247502*/
	Int i,j,jj;
	Doub s;
	if (b.size() != m || x.size() != n) throw("SVD::solve bad sizes");
	VecDoub tmp(n);
	tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
	for (j=0;j<n;j++) {
/*547503 1 2 2 2*/
		s=0.0;
		if (w[j] > tsh) {
/*547504 1 2 2 2*/
			for (i=0;i<m;i++) s += u[i][j]*b[i];
			s /= w[j];
/*547505 1 2 2 2*/
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) {
/*547506 1 2 2 2*/
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
}

void SVD::solve(MatDoub_I &b, MatDoub_O &x, Doub thresh = -1.)
/*247627*/
{
	int i,j,p=b.ncols();
	if (b.nrows() != m || x.nrows() != n || x.ncols() != p)
		throw("SVD::solve bad sizes");
	VecDoub xx(n),bcol(m);
	for (j=0;j<p;j++) {
/*519008 1 2 2 2*/
		for (i=0;i<m;i++) bcol[i] = b[i][j];
		solve(bcol,xx,thresh);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}
/*NR3END svd_solve*/
/*NR3BEG svd_rank*/
Int SVD::rank(Doub thresh = -1.) {
/*247623*/
	Int j,nr=0;
	tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
	for (j=0;j<n;j++) if (w[j] > tsh) nr++;
	return nr;
}

Int SVD::nullity(Doub thresh = -1.) {
/*247624*/
	Int j,nn=0;
	tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
	for (j=0;j<n;j++) if (w[j] <= tsh) nn++;
	return nn;
}

MatDoub SVD::range(Doub thresh = -1.){
/*247625*/
	Int i,j,nr=0;
	MatDoub rnge(m,rank(thresh));
	for (j=0;j<n;j++) {
		if (w[j] > tsh) {
			for (i=0;i<m;i++) rnge[i][nr] = u[i][j];
			nr++;
		}
	}
	return rnge;
}

MatDoub SVD::nullspace(Doub thresh = -1.){
/*247626*/
	Int j,jj,nn=0;
	MatDoub nullsp(n,nullity(thresh));
	for (j=0;j<n;j++) {
		if (w[j] <= tsh) {
			for (jj=0;jj<n;jj++) nullsp[jj][nn] = v[jj][j];
			nn++;
		}
	}
	return nullsp;
}
/*NR3END svd_rank*/
/*NR3BEG svd_decompose*/
void SVD::decompose() {
/*247602*/
	bool flag;
	Int i,its,j,jj,k,l,nm;
	Doub anorm,c,f,g,h,s,scale,x,y,z;
	VecDoub rv1(n);
	g = scale = anorm = 0.0;
/*547603 2 3 3 3*/
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += abs(u[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					u[k][i] /= scale;
					s += u[k][i]*u[k][i];
				}
				f=u[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
					f=s/h;
					for (k=i;k<m;k++) u[k][j] += f*u[k][i];
				}
				for (k=i;k<m;k++) u[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += abs(u[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					u[i][k] /= scale;
					s += u[i][k]*u[i][k];
				}
				f=u[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=u[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += u[j][k]*u[i][k];
					for (k=l-1;k<n;k++) u[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) u[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(abs(w[i])+abs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
/*547604 2 3 3 3*/
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
/*547605 2 3 3 3*/
					v[j][i]=(u[i][j]/u[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
/*547606 2 3 3 3*/
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) u[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
				f=(s/u[i][i])*g;
				for (k=i;k<m;k++) u[k][j] += f*u[k][i];
			}
			for (j=i;j<m;j++) u[j][i] *= g;
		} else for (j=i;j<m;j++) u[j][i]=0.0;
		++u[i][i];
	}
	for (k=n-1;k>=0;k--) {
/*547607 2 3 3 3*/
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
/*547608 2 3 3 3*/
				nm=l-1;
				if (l == 0 || abs(rv1[l]) <= eps*anorm) {
					flag=false;
					break;
				}
				if (abs(w[nm]) <= eps*anorm) break;
			}
			if (flag) {
				c=0.0;
/*547610 2 3 3 3*/
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (abs(f) <= eps*anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=u[j][nm];
						z=u[j][i];
						u[j][nm]=y*c+z*s;
						u[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
/*547611 2 3 3 3*/
				if (z < 0.0) {
/*547612 2 3 3 3*/
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) throw("no convergence in 30 svdcmp iterations");
			x=w[l];
/*547613 2 3 3 3*/
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
/*547614 2 3 3 3*/
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
/*547615 2 3 3 3*/
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=u[jj][j];
					z=u[jj][i];
					u[jj][j]=y*c+z*s;
					u[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}

void SVD::reorder() {
/*247628*/
	Int i,j,k,s,inc=1;
	Doub sw;
	VecDoub su(m), sv(n);
	do { inc *= 3; inc++; } while (inc <= n);
/*547629 6 6 6 6*/
	do {
		inc /= 3;
		for (i=inc;i<n;i++) {
			sw = w[i];
			for (k=0;k<m;k++) su[k] = u[k][i];
			for (k=0;k<n;k++) sv[k] = v[k][i];
			j = i;
			while (w[j-inc] < sw) {
				w[j] = w[j-inc];
				for (k=0;k<m;k++) u[k][j] = u[k][j-inc];
				for (k=0;k<n;k++) v[k][j] = v[k][j-inc];
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for (k=0;k<m;k++) u[k][j] = su[k];
			for (k=0;k<n;k++) v[k][j] = sv[k];

		}
	} while (inc > 1);
	for (k=0;k<n;k++) {
/*547630 6 6 6 6*/
		s=0;
		for (i=0;i<m;i++) if (u[i][k] < 0.) s++;
		for (j=0;j<n;j++) if (v[j][k] < 0.) s++;
		if (s > (m+n)/2) {
			for (i=0;i<m;i++) u[i][k] = -u[i][k];
			for (j=0;j<n;j++) v[j][k] = -v[j][k];
		}
	}
}

Doub SVD::pythag(const Doub a, const Doub b) {
/*241101*/
	Doub absa=abs(a), absb=abs(b);
	return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}
/*NR3END svd_decompose*/
