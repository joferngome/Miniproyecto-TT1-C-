
#include <math.h>
#include "seebattk.h"

//------- recursion algorithms needed by the lambertbattin routine -------


double seebattk(double v)
{


	//-----------------------------  Locals  ------------------------------
	//--------------------  Implementation   ----------------------



	double d[21];
	d[0] =     1.0 /    3.0;
	d[1] =     4.0 /   27.0;
	d[2] =     8.0 /   27.0;
	d[3] =     2.0 /    9.0;
	d[4] =    22.0 /   81.0;
	d[5] =   208.0 /  891.0;
	d[6] =   340.0 / 1287.0;
	d[7] =   418.0 / 1755.0;
	d[8] =   598.00 / 2295.0;
	d[9] =   700.0 / 2907.0;
	d[10]=   928.0 / 3591.0;
	d[11]=  1054.0 / 4347.0;
	d[12]=  1330.0 / 5175.0;
	d[13]=  1480.0 / 6075.0;
	d[14]=  1804.0 / 7047.0;
	d[15]=  1978.0 / 8091.0;
	d[16]=  2350.0 / 9207.0;
	d[17]=  2548.0 /10395.0;
	d[18]=  2968.0 /11655.0;
	d[19]=  3190.0 /12987.0;
	d[20]=  3658.0 /14391.0;

	//----------------- Process Forwards ------------------------


	double sum1   = d[0];
	double delold = 1.0; 
	double termold= d[0];
	double term, del;
	int i = 0;
	while(1){
	    del  = 1.0 / ( 1.0 + d[i+1]*v*delold );
	    term = termold * ( del - 1.0 );
	    sum1 = sum1 + term;
	    i    = i + 1;
	    delold = del;
	    termold= term;
	    if ((i<20) || (fabs(termold) > 0.000001 )){
		break;
	    }
	}


	return sum1;

}
