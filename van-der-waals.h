
double HAM = 1e-10;
  
#define a_VdW(eta,i) (-HAM*(pow(2.*eta[i],-3)-pow(2.*eta[i-1],-3))/Delta )
//#define a_VdW(eta,i) (-HAM*(pow(2.*eta[i+1],-3)-pow(2.*eta[i-1],-3))/(2.*Delta))
//#define a_VdW(eta,i) (3*HAM*pow(2.*eta[i],-4)*2.*(eta[i]-eta[i-1])/(Delta))
//#define a_VdW(eta,i) (3*HAM*pow(2.*eta[i],-4)*(eta[i+1]-eta[i-1])/(Delta))

#include "./hydro-tension.h"
