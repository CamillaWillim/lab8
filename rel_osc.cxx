// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,double* k1,double* k2,double* k3,double* k4);
void bisection(double* const yn, const double* const y0, const double dx,double* k1,double* k2,double* k3,double* k4,double& theta);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 2;
	double dx = 0.1,x=0;
	const double L = 20;
	//double P0=0.1;
	int N = 100;
  double y0[dim];
	
	double yn[dim];
	double k1[dim],k2[dim],k3[dim],k4[dim];
	double theta;
	

  //out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;
  
  for(double p0=0.1;p0<=5;p0+=0.1){
    x=0;
    y0[0] = p0;
    y0[1] = 0;//-p0/sqrt(1+pow(p0,2));
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx,k1,k2,k3,k4);
		if(y0[1]>0 && yn[1]<0)
		  break;
    for(int i=0; i<dim; i++) y0[i] = yn[i];
		//out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	}
	bisection(yn,y0,dx,k1,k2,k3,k4,theta);
    //out << x-dx+theta*dx<< "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;

  out << x-dx+theta*dx<< "\t" << p0 <<endl;
}
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx,double* k1,double* k2,double* k3,double* k4)
{
	const int dim = 2;
	//double k1[dim], k2[dim], k3[dim], k4[dim];

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

  y0[0] = y[1];
  y0[1] = -y[0]/sqrt(1+pow(y[0],2));
}

void bisection(double* const yn, const double* const y0, const double dx,double* k1,double* k2,double* k3,double* k4,double& theta){
  double l=0; // linkes theta
  double r=1; //rechtes theta
  double b1,b2,b3,b4;
  
  
  while(abs(yn[1])>1e-8){
    
    theta = (l+r)/2.0;
    b1 = theta - 3*pow(theta,2)/2.0 + 2*pow(theta,3)/3.0;
    b2 = pow(theta,2) - 2*pow(theta,3)/3.0;
    b3 = b2;
    b4 = -pow(theta,2)/2.0 + 2*pow(theta,3)/3.0;
    
    
  yn[1] = y0[1] + dx*(b1*k1[1] + b2*k2[1] + b3*k3[1] + b4*k4[1]);
  
  if(yn[1]>0)
    l=theta;
  else
    r = theta;
  }
  
  
}
  
  

