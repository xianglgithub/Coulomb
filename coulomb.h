#ifndef COULOMB_H
#define COULOMB_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>   
#include <iostream>  
#include <fstream>

struct my_params 
{
  int num;
  std::vector<double> mass;//=
  std::vector<double> charge;//= 
}; 

int rk8(size_t dims, my_params *para, double t, double t1, int iternum, double hstart, double epsabs, double epsrel, double y[]);

double potential(my_params *para, const double d[]);
  
#endif

