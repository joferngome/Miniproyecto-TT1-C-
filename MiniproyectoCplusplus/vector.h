#include <math.h>

#ifndef _NORMA_
#define _NORMA_

double norm(double v[], int n = 3);
#endif

#ifndef _DOT_
#define _DOT_

double dot(double v1[], double v2[], int n1 = 3, int n2 = 3);

#endif

#ifndef _CROSS_
#define _CROSS_

void cross(double v[], int &n, double v1[], double v2[], int n1 = 3, int n2 = 3);

#endif
