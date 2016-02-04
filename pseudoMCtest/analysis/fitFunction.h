#include <iostream>
#include <TMath.h>

double fitZpmass(double* v, double* p){

  double x = v[0];
  return p[0]*TMath::Exp(p[1]*x + p[2]/x);

}

double divFunc(double* v, double* p){

  double x = v[0];
  return (p[0]*TMath::Exp(p[1]*x+p[2]/x))/(p[3]*TMath::Exp(p[4]*x+p[5]/x));

}

double ErfExp(double x, double c, double offset, double width){
  
  if( width < 1e-2 ) width = 1e-2;
  if( c == 0 ) c = -1e-7;
  return TMath::Exp(c*x)*(1.+TMath::Erf((x-offset)/width))/2.;

}

double integral_ErfExp(const double c, const double offset, double width, const double xmin, const double xmax){

  double width_tmp = width; 
  double minTerm   = 0;
  double maxTerm   = 0;

  if( width < 1e-2 ) width = 1e-2;

  if( c == 0 ){ 
    
    double delta = -1e-7;

    minTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
	       TMath::Erf((2*xmin-delta*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(delta*xmin) * 
	       TMath::Erf((xmin-offset)/width_tmp) - 
	       TMath::Exp(delta*xmin))/-2/delta;

    maxTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
	       TMath::Erf((2*xmax-delta*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(delta*xmax) * 
	       TMath::Erf((xmax-offset)/width_tmp) - 
	       TMath::Exp(delta*xmax))/-2/delta;
        
  }

  else{

    minTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
	       TMath::Erf((2*xmin-c*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(c*xmin) * 
	       TMath::Erf((xmin-offset)/width_tmp) - 
	       TMath::Exp(c*xmin))/-2/c;

    maxTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
	       TMath::Erf((2*xmax-c*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(c*xmax) * 
	       TMath::Erf((xmax-offset)/width_tmp) - 
	       TMath::Exp(c*xmax))/-2/c;
  }

  return (maxTerm-minTerm);

}  

double fitPRmass(double* v, double* p){

  double x = v[0];
  Double_t width_tmp = p[3];
  Double_t binwidth  = p[4];

  if( p[3] < 1e-2 ) width_tmp = 1e-2;

  return p[0]*ErfExp(x,p[1],p[2],width_tmp)/integral_ErfExp(p[1],p[2],width_tmp,40,240)*binwidth;

}

double hollow_fitPRmass(double* v, double* p){

  double x = v[0];
  double width_tmp = p[3];
  double binwidth  = p[4];

  if( p[3] < 1e-2 ) width_tmp = 1e-2;

  return p[0]*ErfExp(x,p[1],p[2],width_tmp)/(integral_ErfExp(p[1],p[2],width_tmp,40,65)+integral_ErfExp(p[1],p[2],width_tmp,145,240))*binwidth;

}
