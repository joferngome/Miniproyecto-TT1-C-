#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include "vector.h"
#include "seebatt.h"
#include "seebattk.h"
#include "lambertbattin.h"
using namespace std;

int main(){

	double v[] = {1,1,1};
	double v1[] = {1,-2,3};
	double v2[] = {1,2,3};
	int n;
	
	//norm
	try{
		cout<<"Norma = "<<norm(v)<<endl;
		cout<<"Norma = "<<norm(v,0)<<endl;
	}catch(const char *msg)
	{
		cerr<<msg<<endl;
	}
	
	//dot
	try{
		cout<<"Producto escalar = "<<dot(v1,v2)<<endl;
		cout<<"Producto escalar = "<<dot(v1,v2,0,6)<<endl;
	}catch(const char *msg)
	{
		cerr<<msg<<endl;
	}
	
	//cross
	try{

		cout<<"Producto vectorial:"<<endl;
		cross(v,n,v1,v2);

		cout<<"   v[0] = "<<v[0]<<endl;
		cout<<"   v[1] = "<<v[1]<<endl;
		cout<<"   v[2] = "<<v[2]<<endl;
		cout<<"   dimension = "<<n<<endl;

	}catch(const char *msg)
	{
		cerr<<msg<<endl;
	}
	cout<<endl;

	//seebatt
	cout<<"SEEBATT"<<endl;
	if(seebatt(0.0)-5.0 <= pow(10,-14))
	{
		cout<<"funciona 1"<<endl;
	}else
	{
		cout<<"fall1"<<endl;
	}





	//seebattk
	cout<<endl<<"SEEBATTK"<<endl;
	if(seebattk(0.0)-0.333333333333333 <= pow(10,-14))
	{
		cout<<"seebatt(0.0) es correcto"<<endl;
	}else
	{
		cout<<"Falla seebatt(0.0)"<<endl;
	}


	
	//Test Lambert Battin
	double r1[3]={20.0e6, 20.0e6, 0};
	double r2[3]={-20.0e6, 10.0e6, 0};
	double tof = 1.0 * 86400;
	char dm[]= "retro";
	LAMBERTBATTIN(r1, r2, dm, tof, v1, v2);
	
	cout<<endl<<"LAMBERT BATTIN"<<endl;
	cout<<"["<< v1[0]<<","<<v1[1]<<","<<v1[2]<<"]"<<endl;
	cout<<"["<< v2[0]<<","<<v2[1]<<","<<v2[2]<<"]"<<endl;

	
	return 0;
}

