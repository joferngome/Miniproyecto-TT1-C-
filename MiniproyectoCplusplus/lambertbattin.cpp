#include <cstring>
#include <iostream>
#include <math.h>
#include "seebatt.h"
#include "seebattk.h"
#include "lambertbattin.h"
#include "vector.h"
using namespace std;


/*-----------------------------------------------------------------------------
*
*                          LAMBERBATTIN
*
*   this subroutine solves Lambert's problem using Battins method. The method
*   is developed in Battin (1987). It uses contiNued fractions to speed the
*   solution and has several parameters that are defined differently than
*   the traditional Gaussian technique.
*
* Inputs:         Description                    Range/Units
*   ro          - IJK Position vector 1          m
*   r           - IJK Position vector 2          m
*   dm          - direction of motion            'pro','retro'
*   Dtsec       - Time between ro and r          s
*
* OutPuts:
*   vo          - IJK Velocity vector            m/s
*   v           - IJK Velocity vector            m/s
*
* Reference:
* Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
* , New York; 3rd edition(2007).
* 
* Last modified:   2015/08/12   M. Mahooti
* 
*------------------------------------------------------------------------------*/



void LAMBERTBATTIN(double ro[],double r[],char *dm, double Dtsec, double vo[], double v[])

{




	double small = 0.000001;
	double mu = 3.986004418e14;   // m3/s2
	double y1 = 0;
	double magr = norm(r);
	double magro = norm(ro);
	double CosDeltaNu= dot(ro,r)/(magro*magr);
	double rcrossr[3];
	int n;
	cross(rcrossr,n,ro,r);
	double magrcrossr = norm(rcrossr);
	
	double SinDeltaNu;
	if(strcmp(dm,"pro") == 0)
	{
	    SinDeltaNu = magrcrossr/(magro*magr);
	}else
	{
	    SinDeltaNu = -magrcrossr/(magro*magr);
	}
	
	double DNu = atan2(SinDeltaNu,CosDeltaNu);
	
	//the angle needs to be positive to work for the long way
	if(DNu < 0.0)
	{
	    DNu = 2.0*PI_+DNu;
	}


	

	double RoR   = magr/magro;
	double eps   = RoR - 1.0;
	double tan2w = 0.25*eps*eps/(sqrt(RoR) + RoR*(2.0+sqrt(RoR)));
	double rp    = sqrt( magro*magr )*(pow(cos(DNu*0.25),2) + tan2w);

	double L;
	if (DNu < PI_)
	{
		L = (pow(sin(DNu*0.25),2) +tan2w)/(pow(sin(DNu*0.25),2) +tan2w +cos(DNu*0.5));
	}else
	{
		L = (pow(cos(DNu*0.25),2) +tan2w -cos(DNu*0.5))/(pow(cos(DNu*0.25),2) +tan2w);
	}
	
	double m    = mu*Dtsec*Dtsec / (8.0*rp*rp*rp);
	double x    = 10.0;
	double xn   = L;
	double chord= sqrt(magro*magro + magr*magr - 2.0*magro*magr*cos(DNu));
	double s    = (magro + magr + chord)*0.5;
	double lim1 = sqrt(m/L);
	
	int Loops= 1;
	double tempx, Denom, h1, h2, b, u, k2, y;
	
	while (1)
	{
	    x    = xn;    
	    tempx= seebatt(x);
	    Denom= 1.0/((1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x)));
	    h1   = pow(L+x,2) * (1.0 + 3.0*x + tempx)*Denom;
	    h2   = m*( x - L + tempx )*Denom;
	  
	    //----------------------- Evaluate CUBIC ------------------
	    b = 0.25*27.0*h2 / pow(1.0+h1,3);
	    
	    if (b < -1.0) //reset the initial condition
	    {
		xn = 1.0 - 2.0*L;
	    }else
	    {
		if (y1 > lim1)
		{
		    xn = xn * (lim1/y1);
		}else
		{
		    u = 0.5*b / (1.0 + sqrt(1.0 + b));            
		    k2 = seebattk(u);
		    y = (( 1.0+h1)/3.0)*(2.0 + sqrt(1.0+b)/(1.0+2.0*u*k2*k2));
		    xn= sqrt(pow((1.0-L)*0.5,2) + m/(y*y)) - (1.0+L)*0.5;
		}
	    }
	    
	    Loops = Loops + 1;
	    y1=sqrt(m/((L+x)*(1.0+x)));
	    
	    if((fabs(xn-x) < small) && (Loops > 30))
	    {
		break;
	    }
	}


	double a = mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );


	// ------------------ Find Eccentric anomalies -----------------
	// ------------------------ Hyperbolic -------------------------


	
	double arg1, arg2, AlpH, BetH, DH, F, GDot, G, Sinv, Cosv, AlpE, BetE, DE, am, ae, be ,tm;
	
	if (a < -small)
	{
	    arg1 = sqrt( s / ( -2.0*a ) );
	    arg2 = sqrt( ( s-chord ) / ( -2.0*a ) );
	    //------- Evaluate f and g functions --------
	    AlpH = 2.0 * asinh( arg1 );
	    BetH = 2.0 * asinh( arg2 );
	    DH   = AlpH - BetH;
	    F    = 1.0 - (a/magro)*(1.0 - cosh(DH) );
	    GDot = 1.0 - (a/magr) *(1.0 - cosh(DH) );
	    G    = Dtsec - sqrt(-a*a*a/mu)*(sinh(DH)-DH);
	}else
	{

	    //------------------------ Elliptical ---------------------
	    if (a > small)
	    {
		arg1 = sqrt( s / (2.0*a));
		arg2 = sqrt((s-chord) / (2.0*a));
		Sinv = arg2;
		Cosv = sqrt( 1.0 - (magro+magr-chord)/(4.0*a));
		//BetE = 2.0*acos(Cosv);
		BetE = 2.0*asin(Sinv);
		if ( DNu > PI_)
		{
		    BetE= -BetE;
		}
		Cosv= sqrt( 1.0 - s/(2.0*a));
		Sinv= arg1;
		am  = s*0.5;
		ae  = PI_;
		be  = 2.0*asin(sqrt((s-chord)/s));
		tm  = sqrt(am*am*am/mu)*(ae - (be-sin(be)));
		if (Dtsec > tm)
		{
		    AlpE= 2.0*PI_-2.0*asin( Sinv );
		}else
		{
		    AlpE= 2.0*asin( Sinv );
		}
		DE  = AlpE - BetE;
		F   = 1.0 - (a/magro)*(1.0 - cos(DE) );
		GDot= 1.0 - (a/magr)* (1.0 - cos(DE) );
		G   = Dtsec - sqrt(a*a*a/mu)*(DE - sin(DE));
	    }else
	    {
		// --------------------- Parabolic ---------------------
		arg1 = 0.0;
		arg2 = 0.0;
		cout<<"Error: a parabolic orbit"<<endl;
	    }


	}
	
	for(int i=0;i<3;i++)
	{
	    vo[i]= ( r[i] - F*ro[i])/G;
	    v[i] = ( GDot*r[i] - ro[i])/G;
	}
}
