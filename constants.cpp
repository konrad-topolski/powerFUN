#include "constants.h"


int nthreads=omp_get_max_threads();
int NGRID=512;
int NYQUIST=int(NGRID/2);
double L=500; //in Mpc/h
double H0=70.4;	 //in (km/s)/Mpc
double H1=70.4;	 //in (km/s)/Mpc
double H2=70.4;	 //in (km/s)/Mpc
double redshift=0;
double redshift1=0;
double redshift2=0;
double omega_matter=0.272;
double omega_lambda=0.728;

std::string filetype="BinaryFile";

//in this code, the minimal LCDM model without radiation  is assumed - only Omega_Matter, Omega_Lambda
// therefore H(a) = H0 * sqrt( Omega_Matter * a^{-3} + Omega_Lambda    )
