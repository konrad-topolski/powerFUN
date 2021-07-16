#include <boost/program_options.hpp>
#include <omp.h>
#include <filesystem>
#include <stdexcept>
#include <cstdlib>
#include <string> 

#define LOC(i,j,k)  (((i)*NGRID*NGRID + (j)*NGRID + (k)))

extern int nthreads;
extern int NGRID;
extern int NYQUIST;
extern double L; 
extern double H0;
extern double H1;
extern double H2;
extern double redshift;
extern double redshift1;
extern double redshift2;
extern double omega_matter;
extern double omega_lambda;

extern std::string filetype;
















