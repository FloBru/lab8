// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
#include "stdlib.h"
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
        const int dim = 2;
	double dx = 1,x=0;
	const double L = 20;
        double y0[dim];
	double yn[dim],y_temp,m;
	double eps=0.0001;
	double theta_beg, theta_end,theta_mit, b[4];
        double k1[dim], k2[dim], k3[dim], k4[dim];
 //       out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	for(double p0 = 0.1; p0 < 5; p0 +=0.1){
	  y0[0] = p0; y0[1] = 0; x =0;
	do
	{
		x += dx;
		
		RKstep(yn, y0, x, dx,k1,k2,k3,k4);
		
		if((yn[1]<0) && (y0[1]>0)) // steigungsänderung von + auf - ?
		   break;
		 
		
		y_temp=y0[1];
		for(int i=0; i<dim; i++) {
		  y0[i] = yn[i];
		  
		  
		}
		
		//s out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
		
	}while(x<=L);
	  m = 1;
	  theta_beg=0;
	  theta_end=1;
	  while(abs(m)>eps)
	{
	  theta_mit=(theta_end + theta_beg)/2;
	  b[0]=theta_mit-3*pow(theta_mit,2)/2+2*pow(theta_mit,3)/3;
	  b[1] = pow(theta_mit,2)-2.*pow(theta_mit,3)/3;
	  b[2] = b[1];
	  b[3]=-pow(theta_mit,2)/2+2.*pow(theta_mit,3)/3;
	  m = y0[1] + (b[0]*k1[1] +b[1]*k2[1] + b[2]*k3[1] + b[3]*k4[1])*dx;
	
	  if(m<0)
	    theta_end=theta_mit;
	  else
	    theta_beg=theta_mit;
	  
	  
	}
	out <<  p0 << "\t" << x+theta_mit*dx -dx <<"\t"<<x<< endl;
        }
	
	out.close();
	
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;
	

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
	double y[2] = { y0[0], y0[1] };

	y0[0] = y[1];;
	y0[1] = -y[0]/sqrt(1+pow(y[0],2));
	
}
