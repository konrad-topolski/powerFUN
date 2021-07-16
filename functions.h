#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <array>
#include <vector>


/*==============================================================================================
Disclaimer: for the calculations of power spectra / auto / cross - correlation functions, prefarably take the
fields averaged by volume
==============================================================================================*/



/*=================================================================================================================================
a special type of parallelised reduction designed to work with the std::vector entries 
===================================================================================================================================*/
#pragma omp declare reduction(vec_double_plus : \
		    std::vector<double> : \
                    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
                    
 #pragma omp declare reduction(vec_int_plus : \
		    std::vector<int> : \
                    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))      
                    
           
/*==================================================================================================================================*/
/*==================================================================================================================================*/  
/*==================================================================================================================================*/



double density_average(fftw_complex* density);
fftw_complex * dens_contrast(fftw_complex* density, double dens_avg);
double win_func(double x);



                 
