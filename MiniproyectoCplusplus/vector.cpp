#include "vector.h"

//Method that calculates the norm of a given vector v.
//v vector (double)
//n dimension  (int)
//Out: return the norm of the given vector v



double norm(double v[], int n)
{
	double sum = 0;
	if(n<=0)
	{
		throw "Vector vacío";
	}else{
		for(int i=0;i<n;i++)
		{
			sum+=v[i]*v[i];
		}
	}
	return sqrt(sum);
}

//Method that calculates the dot product of two given vectors.
//In: v1 , v2 vectors (double)
//In: n1 , n2 dimensions (int)
//Out: returns the dot product the given vectors v1,v2.
double dot(double v1[], double v2[], int n1, int n2)
{
	double sum = 0;
	if(n1<=0 || n2<=0 || n1!=n2)
	{
		throw "Vacíos o dimensión distinta";
	}else{
		for(int i=0;i<n1;i++)
		{
			sum+=v1[i]*v2[i];
		}
	}
	return sum;
}

//Method that calculates the cross product of 2 given vectors
//In: v1, v2 vectors (double)
//In: n1 , n2 dimensiones (int)
//Out: returns the cross product of the given vectors v1,v2
void cross(double v[], int &n, double v1[], double v2[], int n1, int n2)
{
	double sum = 0;
	if(n1<=0 || n2<=0 || n1!=n2)
	{
		throw "Vacíos o dimensión distinta";
	}else{
		v[0]=-v1[2]*v2[1] + v1[1]*v2[2];
		v[1]=v1[2]*v2[0] - v1[0]*v2[2];
		v[2] = -v1[1]*v2[0] + v1[0]*v2[1];
		n=n1;
	} 
}
