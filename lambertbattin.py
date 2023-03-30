#-------------------------------------------------------------------------------
# Name:        LAMBERT BATTIN
# Author:      Jorge Fernández 
# Created:     27/03/2022
# Copyright:   (c) Jorge

#-------------------------------------------------------------------------------

import numpy as np
import math
import doctest 

import unittest



# ------- recursion algorithms needed by the lambertbattin routine
def seebatt( v ):
    """
        >>> seebat(3.76982021737918)
        8.03007718341174
        """
    #-------------------------  implementation   -------------------------
    c=np.zeros(21)
    
    
    c[0] =    0.2
    c[1] =    9.0 /  35.0
    c[2] =   16.0 /  63.0
    c[3] =   25.0 /  99.0
    c[4] =   36.0 / 143.0
    c[5] =   49.0 / 195.0
    c[6] =   64.0 / 255.0
    c[7] =   81.0 / 323.0
    c[8] =  100.0 / 399.0
    c[9] =  121.0 / 483.0
    c[10]=  144.0 / 575.0
    c[11]=  169.0 / 675.0
    c[12]=  196.0 / 783.0
    c[13]=  225.0 / 899.0
    c[14]=  256.0 /1023.0
    c[15]=  289.0 /1155.0
    c[16]=  324.0 /1295.0
    c[17]=  361.0 /1443.0
    c[18]=  400.0 /1599.0
    c[19]=  441.0 /1763.0
    c[20]=  484.0 /1935.0
    sqrtopv=np.sqrt(1.0 + v)
    eta    = v / (1.0 + sqrtopv)**2


    # ------------------- process forwards ----------------------
    delold = 1.0
    termold= c[0]   #eta
    sum1   = termold
    i= 0
    while ((i < 20) and (np.fabs(termold) > 0.00000001)):
        dele  = 1.0 / ( 1.0 + c[i+1]*eta*delold)
        term = termold * (dele - 1.0)
        sum1 = sum1 + term
        i    = i + 1;
        delold = dele
        termold= term

    return (1.0/ ((1.0/(8.0*(1.0+sqrtopv))) * ( 3.0 + sum1 / ( 1.0+eta*sum1 ))))

def seebattk(v):
    """
        >>> seebattk(5)
        0.191489361702128
        """
        

	#-----------------------------  Locals  ------------------------------
	#--------------------  Implementation   ----------------------

    d=np.zeros(21)
    d[0] =     1.0 /    3.0
    d[1] =     4.0 /   27.0
    d[2] =     8.0 /   27.0
    d[3] =     2.0 /    9.0
    d[4] =    22.0 /   81.0
    d[5] =   208.0 /  891.0
    d[6] =   340.0 / 1287.0
    d[7] =   418.0 / 1755.0
    d[8] =   598.00 / 2295.0
    d[9] =   700.0 / 2907.0
    d[10]=   928.0 / 3591.0
    d[11]=  1054.0 / 4347.0
    d[12]=  1330.0 / 5175.0
    d[13]=  1480.0 / 6075.0
    d[14]=  1804.0 / 7047.0
    d[15]=  1978.0 / 8091.0
    d[16]=  2350.0 / 9207.0
    d[17]=  2548.0 /10395.0
    d[18]=  2968.0 /11655.0
    d[19]=  3190.0 /12987.0
    d[20]=  3658.0 /14391.0

	#----------------- Process Forwards ------------------------
    sum1   = d[0]
    delold = 1.0
    termold= d[0]
    i = 0
    while(1):
        dele  = 1.0 / ( 1.0 + d[i+1]*v*delold )
        term = termold * ( dele - 1.0 )
        sum1 = sum1 + term
        i    = i + 1
        delold = dele
        termold= term
        if ((i<20) or (np.fabs(termold) > 0.000001 )):
            break
    return sum1

"""-----------------------------------------------------------------------------
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
*------------------------------------------------------------------------------"""

def LAMBERTBATTIN(ro,r,dm,Dtsec):

    small = 0.000001
    mu = 3.986004418e14   #m3/s2
    y1 = 0
    magr = np.linalg.norm(r)
    magro = np.linalg.norm(ro)
    CosDeltaNu= np.dot(ro,r)/(magro*magr)
    rcrossr=np.cross(ro,r)
    magrcrossr = np.linalg.norm(rcrossr)

    if(dm=="pro"):
        SinDeltaNu = magrcrossr/(magro*magr)
    else:
        SinDeltaNu = -magrcrossr/(magro*magr)

    DNu = math.atan2(SinDeltaNu,CosDeltaNu)

    #the angle needs to be positive to work for the long way
    if(DNu < 0.0):
        DNu = 2.0*np.pi+DNu

    RoR   = magr/magro
    eps   = RoR - 1.0
    tan2w = 0.25*eps*eps/(np.sqrt(RoR) + RoR*(2.0+np.sqrt(RoR)))
    rp    = np.sqrt( magro*magr )*((np.cos(DNu*0.25))**2 + tan2w)

    if (DNu < np.pi):
        L = ((np.sin(DNu*0.25)**2) +tan2w)/((np.sin(DNu*0.25)**2) +tan2w +np.cos(DNu*0.5))
    else:
        L = ((np.cos(DNu*0.25)**2) +tan2w -np.cos(DNu*0.5))/((np.cos(DNu*0.25)**2) +tan2w)

    m    = mu*Dtsec*Dtsec / (8.0*rp*rp*rp)
    x    = 10.0
    xn   = L;
    chord= np.sqrt(magro*magro + magr*magr - 2.0*magro*magr*np.cos(DNu))
    s    = (magro + magr + chord)*0.5
    lim1 = np.sqrt(m/L)

    Loops= 1

    while (1):
        x    = xn
        tempx= seebatt(x)
        Denom= 1.0/((1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x)))
        h1   = ((L+x)**2) * (1.0 + 3.0*x + tempx)*Denom
        h2   = m*( x - L + tempx )*Denom

        #----------------------- Evaluate CUBIC ------------------
        b = 0.25*27.0*h2 / ((1.0+h1)**3)

        if (b < -1.0): #reset the initial condition
            xn = 1.0 - 2.0*L
        else:
            if (y1 > lim1):
                xn = xn * (lim1/y1)
            else:
                u = 0.5*b / (1.0 + np.sqrt(1.0 + b))
                k2 = seebattk(u)
                y = (( 1.0+h1)/3.0)*(2.0 + np.sqrt(1.0+b)/(1.0+2.0*u*k2*k2))
                xn= np.sqrt(((1.0-L)*0.5)**2 + m/(y*y)) - (1.0+L)*0.5

        Loops = Loops + 1
        y1=np.sqrt(m/((L+x)*(1.0+x)))

        if((np.fabs(xn-x) < small) and (Loops > 30)):
            break

    a = mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );

	# ------------------ Find Eccentric anomalies -----------------
	#------------------------ Hyperbolic -------------------------
    if (a < -small):
        arg1 = np.sqrt( s / ( -2.0*a ) )
        arg2 = np.sqrt( ( s-chord ) / ( -2.0*a ) )
        #------- Evaluate f and g functions --------
        AlpH = 2.0 * np.asinh( arg1 )
        BetH = 2.0 * np.asinh( arg2 )
        DH   = AlpH - BetH
        F    = 1.0 - (a/magro)*(1.0 - np.cosh(DH) )
        GDot = 1.0 - (a/magr) *(1.0 - np.cosh(DH) )
        G    = Dtsec - np.sqrt(-a*a*a/mu)*(np.sinh(DH)-DH)
    else:
        if(a>small):
            #------------------------ Elliptical ---------------------
            arg1 = np.sqrt( s / (2.0*a))
            arg2 = np.sqrt((s-chord) / (2.0*a))
            Sinv = arg2
            Cosv = np.sqrt( 1.0 - (magro+magr-chord)/(4.0*a));
    		#BetE = 2.0*acos(Cosv)
            BetE = 2.0*math.asin(Sinv)
            if ( DNu > np.pi):
                BetE= -BetE
            Cosv= np.sqrt( 1.0 - s/(2.0*a))
            Sinv= arg1
            am  = s*0.5
            ae  = np.pi
            be  = 2.0*math.asin(np.sqrt((s-chord)/s))
            tm  = np.sqrt(am*am*am/mu)*(ae - (be-np.sin(be)))
            if(Dtsec > tm):
                AlpE= 2.0*np.pi-2.0*math.asin( Sinv )
            else:
                AlpE= 2.0*math.asin( Sinv )
            DE  = AlpE - BetE
            F   = 1.0 - (a/magro)*(1.0 - np.cos(DE) )
            GDot= 1.0 - (a/magr)* (1.0 - np.cos(DE) )
            G   = Dtsec - np.sqrt(a*a*a/mu)*(DE - np.sin(DE))
        else:
            #--------------------- Parabolic ---------------------
            arg1 = 0.0
            arg2 = 0.0
            print("Error: a parabolic orbit")

    vo=np.zeros(3)
    v=np.zeros(3)
    for i in range(2):
        vo[i]= ( r[i] - F*ro[i])/G
        v[i] = ( GDot*r[i] - ro[i])/G

    return vo, v



doctest.testmod()




