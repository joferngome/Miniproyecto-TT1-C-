#include <iostream>
#include <math.h>
#include <stdio.h>
#include "seebatt.h"
#include "seebattk.h"
#include "lambertbattin.h"
using namespace std;

int test_run = 0;

#define FAIL() printf("\n fallo en %s() linea %d\n",__func__,__LINE__)
#define _assert(test) do { if(!(test)) { FAIL(); return 1;}} while(0)
#define _verify(test) do {int r=test(); test_run++;if(r) return r;} while(0)

int seebatt_01()
{
	_assert(seebatt(0.0)-5.0 <= pow(10,-14));
	return 0;
}

	
int seebattk_01()
{
	_assert(seebattk(0.0)-0.333333333333333<= pow(10,-14));
	return 0;
}


int lambertbattin_01()
{
	double r1[3]={20.0e6, 20.0e6, 0};
	double r2[3]={-20.0e6, 10.0e6, 0};
	double tof = 1.0 * 86400;
	char dm[]= "retro";
	double v1[3], v2[3];
	double lim = pow(10,-10);
	
	LAMBERTBATTIN(r1, r2, dm, tof, v1, v2);
	
	_assert(v1[0]-4144.30717367666 <= lim);
	_assert(v1[1]+1571.15318557576 <= lim);
	_assert(v1[2]-0 <= lim);
	
	_assert(v2[0]-3223.39508300487 <= lim);
	_assert(v2[1]-4103.76281774998 <= lim);
	_assert(v2[2]-0 <= lim);
	
	return 0;
}



int all_tests()
{
	_verify(seebatt_01);

	
	_verify(seebattk_01);

	
	_verify(lambertbattin_01);

	
	return 0;
}

int main()
{

}
