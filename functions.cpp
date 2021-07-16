#include "functions.h"
#include "constants.h"


double density_average(fftw_complex* density){
double result=0;

#pragma omp parallel num_threads(nthreads)
{

#pragma omp single //even though this says pragma omp single, it is executed within a parallel section
{
std::cout << "Threads working for average density calculation: " << omp_get_num_threads() << std::endl;
}

#pragma omp parallel for reduction(+:result)
for(int i=0; i<NGRID*NGRID*NGRID;i++){
result+=density[i][0];}
}

result/=NGRID*NGRID*NGRID;
return result;
}

/*=====================================================================================
The omp_get_num_threads() function returns the number of threads 
that are currently in the team executing the parallel region from which it is called.
The omp_get_num_threads() therefore gets the number of threads 
of the parallel region within which it is enclosed, regardless of omp single, 
which is there just to print it once, not nthreads-number-of-times. 
======================================================================================*/

fftw_complex* dens_contrast(fftw_complex* density,double dens_avg){
#pragma omp parallel for num_threads(nthreads)
	for(int i=0; i<NGRID*NGRID*NGRID;i++){
		density[i][0]=density[i][0]/dens_avg -1.;
	}
return density;
}

double win_func_sinc(double x){		//the definition of the used window function
return sin(x)/x;
}






