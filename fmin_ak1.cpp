// fmin_ak1.cpp -- This program finds the maximum of the following function:
//
// asech[sin(ET4)] - 0.81732634658 * e
//
// The main minimum-finding routine is a translation of the FORTRAN routine
// FMIN.F, which is a slightly modifed version of the ALGOL 60 procedure
// LOCALMIN written by Richard Brent.
//
// The equation chosen for use in this sample program includes a
// transcendental function for which a numerical root-finer is required to
// solve. The root-finding routine used is a translation of the FORTRAN
// routine FZERO.F, written by L. F. Shampine (SNLA) and H. A. Watts (SNLA),
// based upon a method by T. J. Dekker.
//
// References:
//
// Brent, Richard
// "Algorithms for Minimization Without Derivatives"
// Prentice - Hall, Inc.
// 1973
//
// Shampine, L.F. (SNLA) and H.A.Watts (SNLA)
// "FZERO, A Root-solving Code"   Report SC - TM - 70 - 631
// Sandia Laboratories
// September, 1970
//
// Dekker, T.J.
// "Finding a Zero by Means of Successive Linear Interpolation"
// "Constructive Aspects of the Fundamental Theorem of Algebra"
// edited by B.Dejon and P.Henrici
// Wiley - Interscience
// 1969
//
// To distinguish this version of fmin from other translations,
// an '_ak1' suffix has been appended to its name.
//
// 2 August 2017
//
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
//

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cctype>

using namespace std;

#define PI 3.14159265358979323846
#define HALF_PI 0.5 * PI
#define SQ_GOLDEN_INV 0.38196601125010515179541	//Squared inverse of Golden Ratio

struct ZTYPE {
	double b, c, re, ae;
	short iflag, kount;
	double mAnom;
};

double kep_f(double meanAnom, double eccAnom, double ecc);
double f(double x);
double sign(double y);
void Compute_ET4(ZTYPE *zd, double ecc);
double fmin(ZTYPE *zd, double tol);

double kep_f(double meanAnom, double eccAnom, double ecc){
	return eccAnom - ecc*sin(eccAnom) - meanAnom;
}

double f(double x){	//The function to be minimized

	// **************** VARIABLES FOR KEPLER'S EQUATION *************

	ZTYPE KepSolvVar;		// Variable for Kepler Solver Routine

	int seg_ind;			// Section indicator, to indicate which segment of the piecewise bounds to use.
	double E_QuarterT, sin_ET4, temp1;

	// KepSolvVar.mAnom = 0.5 * PI;
	KepSolvVar.mAnom = HALF_PI;
	KepSolvVar.ae = KepSolvVar.re = DBL_EPSILON;

	//KepSolvVar.b = KepSolvVar.mAnom + 0.5 * x;
	//KepSolvVar.c = KepSolvVar.mAnom + x;
	//
	// Use new limits
	//KepSolvVar.b = KepSolvVar.mAnom + x * (1.0 - mDot * x);
	//KepSolvVar.c = KepSolvVar.mAnom + x * cos(DOTTIE *x);
	//
	// *******************************************************
	// Use the latest piecewise linear bounds

	seg_ind = static_cast<int> (10.0 * x);

	switch (seg_ind){
	case 0: if (x < 0.05){
				KepSolvVar.b = 1.0 - 0.0249325256* x;
				KepSolvVar.c = KepSolvVar.b + 0.0006;
			}
			else {
				KepSolvVar.b = 1.0024532 - 0.07399894* x;
				KepSolvVar.c = KepSolvVar.b + 0.00060;
			}
			break;
	case 1:	KepSolvVar.b = 1.0092853 - 0.1423205 * x;
			KepSolvVar.c = KepSolvVar.b + 0.00216;
			break;
	case 2:	KepSolvVar.b = 1.0246502 - 0.21914456 * x;
			KepSolvVar.c = KepSolvVar.b + 0.001670;
			break;
	case 3:	KepSolvVar.b = 1.0414315 - 0.2750822 * x;
			KepSolvVar.c = KepSolvVar.b + 0.00115;
			break;
	case 4:	KepSolvVar.b = 1.0555247 - 0.3103152 * x;
			KepSolvVar.c = KepSolvVar.b + 0.00065;
			break;
	case 5:	KepSolvVar.b = 1.064431 - 0.3281283 * x;
			KepSolvVar.c = KepSolvVar.b + 0.00029;
			break;
	case 6:	KepSolvVar.b = 1.067247 - 0.3328443 * x;
			KepSolvVar.c = KepSolvVar.b + 0.00005;
			break;
	case 7:	KepSolvVar.c = 1.064235 - 0.3285209 * x;
			KepSolvVar.b = KepSolvVar.c - 0.0001;
			break;
	case 8:	KepSolvVar.c = 1.056150 - 0.3184143 * x;
			KepSolvVar.b = KepSolvVar.c - 0.0002;
			break;
	case 9:
	default:KepSolvVar.c = 1.043999 - 0.3049128 * x;
			KepSolvVar.b = KepSolvVar.c - 0.0002;
			break;
	} // End switch

	KepSolvVar.b = KepSolvVar.mAnom + x * KepSolvVar.b;
	KepSolvVar.c = KepSolvVar.mAnom + x * KepSolvVar.c;

	Compute_ET4(&KepSolvVar, x);

	E_QuarterT = KepSolvVar.b;
	sin_ET4 = sin(E_QuarterT);

	temp1 = acosh(1/sin_ET4) - 0.81732634657933280881822 * x;

	return -temp1;
}

double sign(double y)						//If y < 0, return -1, else +1
{
	return ((y < 0) ? -1 : 1);
}

void Compute_ET4(ZTYPE *zd, double ecc)
{/*	This function determines the zero of the function specified in f, above. The method
 employed is an efficient combination of bisection and secant rules.

 Parameters
 b				lower bound of interval in which zero to be found;
 returns as zero to function
 c				upper bound of interval in which zero to be found
 re			relative error
 ae			absolute error
 iflag			status code
 =1	The zero is within the requested tolerance, the interval has collapsed
 to the requested tolerance, the function changes sign over the interval,
 and the function decreased in magnitude as the interval collapsed.
 =2	A zero has been found, but the interval has not collapsed to the
 requested tolerance.
 =3	Possibly near a singular point. The interval has collapsed to the
 requested tolerance and the function changes sign over the interval,
 but the magnitude of the function increased as the interval collapsed.
 =4	The function does not change sign over the specified interval, which
 has collapsed to the requested tolerance. Possibly near a minimum of
 the function or a zero of even multiplicity.
 =5	More than MAXIT function evaluations used.
 kount			number of function calls

 Authors:	Shampine, L.F., SNLA
 Watts, H.A., SNLA */

	static const short MAXIT = 100;
	double RW = 2.0 * DBL_EPSILON;  //Machine epsilon for type double
	double a, acbs, acmb, cmb, fa, fb, fc, fx, fz, p, q, t, tol, z;
	short ic = 0;

	//Initialize and do some checks

	zd->kount = 0;

	if (RW < zd->re) RW = zd->re;

	z = zd->c;
	t = zd->b;

	if (z == t){  // Check if b and c are equal
		cout << "The interval endpoints that were input are equal." << endl;
		cout << "Please try again with a non-zero interval." << endl;
		cout << "Routine aborted." << endl;
		zd->iflag = 4;
		return;
	}

	fb = kep_f(zd->mAnom, t, ecc);
	zd->kount = 1;
	if (fabs(fb) < DBL_EPSILON){ // Zero at b
		zd->iflag = 2;
		return;
	}

	z = t + 0.5 * (z - t);

	fc = fz = kep_f(zd->mAnom, z, ecc);
	zd->kount = 2;

	if (sign(fz) == sign(fb)) {
		//cout << "ecc = " << ecc << ".\n";
		//cout << "sign fb == sign fz.\n";
		t = zd->c;
		fc = kep_f(zd->mAnom, t, ecc);
		zd->kount = 3;

		if (fabs(fc) < DBL_EPSILON){ // Zero at c
			zd->b = t;
			zd->iflag = 2;
			return;
		}

		if (sign(fz) != sign(fc)) {
			zd->b = z;
			fb = fz;
		}
		else {
			cout << "ecc = " << ecc << ".\n";
			cout << "The function sign does not seem to change over the interval." << endl;
			cout << "Please try again with endpoints such that the sign of the function changes over the interval." << endl;
			cout << "Routine aborted." << endl;
			zd->iflag = 4;
			return;
		}
	}
	else zd->c = z;

	a = zd->c;
	fa = fc;
	acbs = fabs(a - zd->b);

	fx = fabs(fb);
	if (fx < fabs(fc)) fx = fabs(fc);

	do	{
		// Arrange so fabs(f(b)) LE fabs(f(c))
		if (fabs(fc) < fabs(fb))		//Interchange if necessary
		{
			a = zd->b;
			fa = fb;
			zd->b = zd->c;
			fb = fc;
			zd->c = a;
			fc = fa;
		}

		cmb = 0.5 * (zd->c - zd->b);
		acmb = fabs(cmb);
		tol = RW * fabs(zd->b) + zd->ae;

		//Test-stopping criteria and function count

		if (acmb <= tol) break;

		if (fb == 0){
			zd->iflag = 2;
			return;
		}

		if (zd->kount >= MAXIT) {
			zd->iflag = 5;
			return;
		}

		/*Calculate new iterate implicitly as b + p/q, where p is arranged to be
		>= 0. This implicit form is used to prevent overflow.*/

		p = (zd->b - a)*fb;
		q = fa - fb;

		if (p < 0) {
			p = -p;
			q = -q;
		}

		/*Update a and check for satisfactory reduction in the size of the bracketing
		interval. If not, perform bisection.*/

		a = zd->b;
		fa = fb;
		++ic;

		if ((ic >= 4) && (8 * acmb >= acbs))
			zd->b = 0.5 * (zd->c + zd->b);				//Use bisection
		else
		{
			if (ic >= 4) { ic = 0; acbs = acmb; }
			if (p <= tol*fabs(q))			//Test for too small a change
				zd->b += tol*sign(cmb);
			else							//Root between b and (b + c)/2
			{
				if (p < cmb*q) zd->b += p / q;	//Use secant rule
				else zd->b = 0.5 * (zd->c + zd->b);
			}
		} // End else !((ic >= 4) && (8 * acmb >= acbs))

		//Have now computed new iterate,b.
		fb = kep_f(zd->mAnom, zd->b, ecc);
		zd->kount += 1;

		//Decide if next step interpolation or extrapolation.

		if (sign(fb) == sign(fc)) { zd->c = a; fc = fa; }
	} while (zd->kount < MAXIT);				//End while loop

	if (sign(fb) == sign(fc)) {
		zd->iflag = 4;
	} // end if (sign(fb) == sign(fc))
	else {// else (sign(fb) != sign(fc))

		if (fabs(fb) > fx)  zd->iflag = 3;
		else  zd->iflag = 1;

	} // end else (sign(fb) != sign(fc))

	return;
}											//End Compute_ET4

double fmin(ZTYPE *zd, double tol){
 // This function determines an approximation to the point where f attains a
 // minimum on the interval (ax, bx). The method used is a combination of
 // golden section search and successive parabolic interpolation.
 //
 // Parameters
 // ax		left endpoint of initial interval = zd->b
 // bx		right endpoint of initial interval = zd->c
 // tol		desired length of the interval of uncertainty of the final result

	double a, b, c, d, e, eps, xm, p, q, r, tol1, tol2, u, v, w, fu, fv, fw, fx, x;
	static const int FMINIT = 500;	//Maximum number of iterations

	//c is the squared inverse of the golden ratio
	//c = (3 - sqrt(5)) / 2;
	c = SQ_GOLDEN_INV;

	eps = sqrt(DBL_EPSILON);		//Machine epsilon for type double

	//Initialization
	a = zd->b;
	b = zd->c;
	e = 0;
	w = x = v = a + c*(b - a);
	fv = fw = fx = f(x);

	//Start main loop
	for (short i = 0; i < FMINIT; ++i) {

		xm = 0.5 * (a + b);
		tol1 = eps*fabs(x) + tol / 3;
		tol2 = tol1 * 2;

		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) break; // Conditions met; leave main loop

		//Is Golden Section necessary?
		if (fabs(e) <= tol1) {
			if (x >= xm) e = a - x;	//Golden-Section step
			else e = b - x;
			d = c*e;
		}
		else {						//Else fit parabola
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2 * (q - r);
			if (q > 0) p = -p;
			q = fabs(q);
			r = e;
			e = d;

			//Is parabola acceptable?
			if ((fabs(p) >= fabs(0.5*q*r)) || (p <= (q*(a - x))) || (p >= q*(b - x))) {
				if (x >= xm) e = a - x;		//Golden-Section step
				else e = b - x;
				d = c*e;
			}
			else {
				d = p / q;				//Parabolic interpolation step
				u = x + d;
				//f must not be evaluated too close to ax or bx
				if ((u - a) < tol2) d = tol1*sign(xm - x);
				if ((b - u) < tol2) d = tol1*sign(xm - x);
			}		//End else, compund OR statement, p >= q*r/2

		}		//End else fabs(e) <= tol1

		if (fabs(d) >= tol1) u = x + d;
		else u = x + tol1*sign(d);
		fu = f(u);

		//Update a, b, v, w, and x
		if (fu <= fx) {
			if (u >= x) a = x;
			else b = x;
			v = w;
			fv = fw;
			w = x;
			fw = fx;
			x = u;
			fx = fu;
			continue;
		}	//End if fu <= fx

		if (u < x) a = u;
		else b = u;

		if ((fu <= fw) || (w == x)) {
			v = w;
			fv = fw;
			w = u;
			fw = fu;
			continue;
		}

		if ((fu <= fv) || (v == x) || (v == w)) {
			v = u;
			fv = fu;
		}
	}						//End main (for) loop
	return x;
}							//End fmin		

int main()
{
	char rflag;			//Readiness flag

	cout << "                       fmin_ak1 (2 August 2017)" << endl;
	cout << "======================================================================" << endl;
	cout << "\nThis program finds where the following function takes a maximum:" << endl << endl;
	cout << "f = asech[sin(ET4)] - 0.81732634658 * e" << endl << endl;
	cout << "The results are calculated to double precision--" << DBL_DIG << " decimal places." << endl << endl;
	cout << "\nEverything ready? If yes, press y." << endl;
	cout << "Otherwise, press any other key." << endl;
	cin >> rflag;
	if (toupper(rflag) == 'Y') {

		ZTYPE dumvar;	// Variable for fmin
		double res;		// Result

		cout.precision(DBL_DIG);

		dumvar.b = 0.45;
		dumvar.c = 0.55;
		//dumvar.ae = dumvar.re = DBL_EPSILON; // Not used in fmin

		//cout << "A maximum for the function will be sought on the interval" << endl;
		//cout << dumvar.b << " to " << dumvar.c << endl << endl;

		res = fmin(&dumvar, DBL_EPSILON);

		cout << "The point at which f is a maximum is e = " << res << endl;
		cout << "The value of f at this point is f = " << -f(res) << endl;

	}	//End if ready
	else cout << "Not ready. Try again when ready with information.\n";

	cout << "\nEnter any key to continue." << endl;
	cin >> rflag;
	return 0;
}	//End main program
