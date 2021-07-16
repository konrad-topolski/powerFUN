/* This simple program serves to calculate the PDFs (probability distribution functions), 
   CDFs (cumulative distribution functions), averages, mean standard deviations, kurtoses 
   (potentially higher moments too?) from  a binary file containing a given field 
   Features:

   (*) PDF of a field (both with and without additional deconvolution through the top-hat window function in Fourier space):
       supported modes of operation chosen:
		a) log(x/x_avg)  //dedicated to density contrast plus one, which is, as a matter of fact, just rho/rho_avg = delta+1//LOG_X_OVER_XAVG
		b) log(x)    //dedicated to fields that we know to be non-negative (e.x. norm of velocity |v|) 	//LOG_X
		c) log(x+|x_min|) // LOG_X_PLUS_XMIN dedicated to fields that exhibit most informative behaviour when plotted in logspace, but attain negative values (e.x. velocity divergence) 
		d) x         //dedicated to fields whose PDF is best represented in linear space (e.x. velocity average, vorticity average) //LINEAR_X
 
   (**) CDF of a field (cumulative distribution function)
   (***) Negentropy of a random field 
   (****) Central moments of the distrubutions above (mean, deviation, skewedness, kurtosis)

   
   
   Compilation: 
   USE THE PrgEnv-gnu (with the GNU C++ compiler version 8.X)

   Usage: 

   ./statistics_calc binary_file output_pathname MODE NGRID BOXSIZE PDF_BINS NBINS RSTART REND FILTER_TYPE H0 REDSHIFT

   Names of available modes of operation:
	- LOG_X
	- LOG_X_PLUS_XMIN
	- LOG_X_OVER_XAVG
	- LINEAR_X

   Filter types:
	- TopHat
	- Gaussian
	- None

   Example:
./statistics_calc /home/topolski/dtfe_files/512x500Mpc/2/512x500Mpc_R02_W7_LCDM_z1.130.DTFE512.a_den  /home/topolski/CCG_CODES/statistics_calc/testing LOG_X_OVER_XAVG 512 500 1000 5 1 5 Gaussian 70.4 1.130
MOST RECENT UPDATE:	19.05.2021
AUTHOR:			Konrad Topolski
@:			k.topolski2@student.uw.edu.pl

Implementation of smoothing filters inspired by code written by:
AUTHOR:			Wojciech Hellwing
@: 			hellwing@cft.edu.pl
 																       */

/*===================================================================================================================================================
PREPROCESSOR COMMANDS 
===================================================================================================================================================*/

/*================================================
 Necessary definitions. Do not comment it out.
===================================================*/
#define LOC(ix,iy,iz)  (((ix)*NGRID*NGRID + (iy)*NGRID + (iz))    )

/*====================================================
 Enable or disable some functionalities of the program 
=======================================================*/
#define smoothing_filters //enable the calculation of fields with smoothing in Fourier space (by a top-hat function)
#define moments_calc //enable the calculation of central moments
#define negentropy_calc //enable the calculation of negentropy of a d - the measure of its total deviation from gaussian random distribution

/*====================================================
 Technical aspects of running the program 
=======================================================*/

//#define FIXED_ABS  //enable fixed max and min values of the PDF bins 

#define ARTIFICIAL_POSITIVE  //after performing the filtering (or working on a dominantly positive field), erase all the nonpositive values from the field in LOG_X and LOG_X_OVER_XAVG modes

//#define NEGENTROPY_LOG //instead of calculating the negentropy with the referential Gaussian function taking moments from distribution X, take the moments from f(X);  works only if NOT working in LINEAR_X mode   

/*===================================================================================================================================================
HEADERS
===================================================================================================================================================*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <numeric> 
#include <omp.h>

#ifdef smoothing_filters
#include <fftw3.h>
#endif 


#ifdef negentropy_calc
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_deriv.h>
#endif	
//#include <execution> NOT SUPPORTED BY GNU C++ 8.3, ONLY AS OF 9.1 PROBABLY




#pragma omp declare reduction(vec_int_plus : \
		    std::vector<int> : \
                    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))      
                    
/*===================================================================================================================================================
FUNCTION DECLARATIONS
===================================================================================================================================================*/
std::array<double,4> moments(const std::vector<float>& field); //return the mean, deviation, skewedness and kurtosis

void get_PDF_CDF(float rbin, std::string filtertype, 
		 std::vector<float> field, std::vector<float>& PDF_x, 
		 std::vector<int>& PDF_counts, std::vector<float>& PDF, std::vector<float>& CDF);

#ifdef smoothing_filters
void top_hat_f(int NGRID, float fsize, fftw_complex* Fourierfield);
void gaussian_f(int NGRID, float fsize, fftw_complex* Fourierfield);
#endif 

#ifdef negentropy_calc
float negentropy(std::vector<float> field);
static double pdfspline(double x, void * p);
static double integrand_f(double x, void* p);
static double gaussian(double x, void * p);
struct GAUSS_PARAMS { double MEAN; double SIGMA;};
#endif
/*===================================================================================================================================================
GLOBAL VARIABLES
===================================================================================================================================================*/
std::string binary_file, output, mode, filtertype;
int PDF_BINS=1000; //number of pdf bins 
int NGRID;
int NBINS=1000;  //number of 
int NYQUIST;
long int FFTNORM;
double BOXSIZE, RSTART, REND;
float REDSHIFT, H0;
int nthreads=omp_get_max_threads();


#ifdef negentropy_calc
gsl_interp_accel *acc_pk;
gsl_spline *pkl_spline;
float NEGENTROPY;
#endif

float ABS_MIN=-pow(10,30);
float ABS_MAX=pow(10,30);
float ABS_EPSILON=pow(10,-3);

/*===================================================================================================================================================
THE MAIN FUNCTION 
===================================================================================================================================================*/

int main(int argc, char* argv[]){
 
if(argc<11){
	throw std::invalid_argument("Not enough arguments supplied!");
        }


   binary_file = argv[1];
   output= argv[2];
   mode = argv[3];
   NGRID    = atoi(argv[4]);
   BOXSIZE = atof(argv[5]);
   PDF_BINS=atoi(argv[6]);
   NBINS = atoi(argv[7]);
   RSTART = atof(argv[8]);
   REND = atof(argv[9]);
   filtertype = argv[10];
   
   if(argc==12){H0=atof(argv[11]);}			//an optional argument
   else{H0=70.4;}

   if(argc==13){REDSHIFT=atof(argv[12]);}                     //an optional argument
   else{REDSHIFT=0.0;}
		
   NYQUIST = NGRID/2;
   FFTNORM = pow(NGRID,3);


if(!std::filesystem::exists(binary_file)){
	throw std::invalid_argument("Invalid name of the input file! Check the name!");
	}
	
if( mode.compare("LOG_X") != 0 & mode.compare("LINEAR_X") != 0 & mode.compare("LOG_X_OVER_XAVG") != 0 & mode.compare("LOG_X_PLUS_XMIN")!=0 ){
	throw std::invalid_argument("Unrecognised mode of operation chosen!");
	}

if(RSTART>REND || RSTART <=0 || REND <=0){
	throw std::range_error("The range of smoothing radius is incorrect!");
	}

if(filtertype.compare("Gaussian") != 0 & filtertype.compare("TopHat") != 0 & filtertype.compare("None") != 0){
	throw std::invalid_argument("Unrecognised filter type! It must be one of the following: \n TopHat \n Gaussian \n None"); 
	}

std::cout << "The binary file name and location is: " << binary_file << std::endl;
std::cout << "The output name and location: " << output << std::endl;
std::cout << "Grid size: " << NGRID << std::endl;
std::cout << "Box size: " << BOXSIZE << std::endl;
std::cout << "The calculations will yield the PDFs and CDFs in the following mode of operation:" <<  mode << " [LINEAR_X = x, LOG_X = log(x), LOG_X_PLUS_XMIN = log(x+|x_min|), LOG_X_OVER_XAVG = log(x/xavg)] "<< std::endl;
std::cout << "The number of bins for the PDFs and CDFs is: " << PDF_BINS << std::endl;
std::cout << "The smoothing radius range is: [" << RSTART << "," << REND << "]" << " Mpc/h" << std::endl;
std::cout << "Overall, the calculations will proceed with: " << NBINS << " smoothing kernels of different radii" << std::endl;
std::cout << "The chosen filter is: " << filtertype << std::endl;

if(argc==12){std::cout<< "The value of Hubble constant is: " << H0 << std::endl;}
if(argc==13){std::cout<< "The value of redshift is: " << REDSHIFT << std::endl;}
                    

#ifdef FIXED_ABS
std::cout << "During compilation time, the incorporation of fixed absolute and minimum values for any field and mode of operation was chosen" << std::endl;
std::cout << "ABS_MIN=" << ABS_MIN << ", ABS_MAX=" << ABS_MAX ", ABS_EPSILON=" << ABS_EPSILON << std::endl;
#endif

#ifdef smoothing_filters
std::cout << "The operation of smoothing the fields has been chosen during compilation time and the program will perform FFTs if you've chosen appropriate filters" << std::endl;
if(filtertype.compare("None")==0){
        std::cout<< "However, the filter chosen is: " << filtertype <<". Is this intentional?" <<std::endl;
        }
#endif

#ifdef ARTIFICIAL_POSITIVE
std::cout << "The option of artificial positivity is enabled. After the Fourier transform and in the LOG_X or LOG_X_OVER_XAVG mode, if the field attains some negative values, these will be removed" << std::endl;
#endif

std::vector<float> field(FFTNORM);

std::ifstream fin(binary_file,std::ios::binary | std::ios::in);
fin.read(reinterpret_cast<char*>(field.data()), field.size()*sizeof(float));
fin.close();


std::cout << "The data loaded is at least: " << (FFTNORM*sizeof(float))/pow(1024,3) << "GB in size" << std::endl;

/*========================================================================================
        There will be PDF_BINS bins of the histogram
        There will be PDF_BINS+1 values corresponding to each bins upper and lower value,
        i.e. for each bin there are two values: x_{i}, x_{i+1}, that
        describe the range of its contents
========================================================================================*/


std::vector<int> PDF_counts(PDF_BINS+1,0); //unnormalized probability distribution function(just counts)
std::vector<float> PDF(PDF_BINS+1,0); //probability distributution function, i.e. the above, but normalized to 1 in the sum/integral 
std::vector<float> CDF(PDF_BINS+1,0); //cumulative distribution function 
std::vector<float> PDF_x(PDF_BINS+1,0); //the values of the field investigated (to be represented on the x axis), and given in the first column of any file

#ifdef moments_calc
std::array<double,4> MOMENTS_RAW;
std::array<double,4> MOMENTS;
#endif
 
double rbin=RSTART;
double rstep=(REND-RSTART)/(NBINS-1);


//If the filtertype was chosen to be either Gaussian or TopHat, a suitable range of smoothing kernels of different radii will be covered;
// otherwise, the output will only be one realization (no smoothing)

if( (filtertype.compare("Gaussian") == 0)   || (filtertype.compare("TopHat")==0)    ){

for (int counter=0; counter < NBINS; counter++){

FFTNORM=pow(NGRID,3);

std::cout << "Currently, the smoothing radius of index: " << counter << " is R=" << rbin << " Mpc/h" << std::endl;  
std::cout << "Zeroing the PDF_x, PDF_counts, PDF, CDF arrays from previous iteration..." << std::endl;

std::fill(PDF.begin(), PDF.end(), 0);
std::fill(PDF_x.begin(), PDF_x.end(), 0);
std::fill(PDF_counts.begin(), PDF_counts.end(), 0);
std::fill(CDF.begin(), CDF.end(), 0);

std::cout << "Calculating the PDF..." << std::endl;

get_PDF_CDF(rbin, filtertype, field, PDF_x, PDF_counts, PDF, CDF);

std::cout << "Writing the data  to textfile..." << std::endl;


#ifdef moments_calc
MOMENTS_RAW=moments(field);
MOMENTS=moments(PDF_x);
#endif


std::ofstream outfile;
std::stringstream ss;
ss << std::fixed << std::setprecision(2) << rbin;
std::string rbin_cut = ss.str();
outfile.open(output+"_"+mode+"_"+filtertype+"_RBIN_"+rbin_cut+"_stats", std::ios::trunc | std::ios::out);

if(outfile.is_open()){



	outfile << "# PDF Data generated for a simulation of size: " << BOXSIZE << "Mpc/h" << " on a grid of " << NGRID << " points" << std::endl;
	outfile << "# Smoothing kernel was of type: " << filtertype << " with smoothing radius: " << rbin << " Mpc/h" << std::endl;
	outfile << "# Overall there are " << PDF_BINS << " bins of PDF in this file and the field values have a step dictated by mode " << mode << std::endl;

	#ifdef moments_calc
	outfile << "# Central moments of the RAW input x field: " << std::endl;
	outfile << "# Mean: " << MOMENTS_RAW[0] << " " << "Deviation: " << MOMENTS_RAW[1] << " " << "Skewedness: " << MOMENTS_RAW[2] << " " << "Kurtosis: " << MOMENTS_RAW[3]  << std::endl;
	outfile << MOMENTS_RAW[0] << " " << MOMENTS_RAW[1] <<  " " << MOMENTS_RAW[2]  << " " << MOMENTS_RAW[3]  << std::endl;
	outfile << "# Central moments of the input x field as expressed in " << mode << " mode of operation, so for f(x) in one of the forms: " << std::endl;
	outfile << "# [LINEAR_X = x, LOG_X = log(x), LOG_X_PLUS_XMIN = log(x+|x_min|), LOG_X_OVER_XAVG = log(x/x_avg)] " << std::endl; 
	outfile << "# Mean: " << MOMENTS[0] << " " << "Deviation: " << MOMENTS[1] << " " << "Skewedness: " << MOMENTS[2] << " " << "Kurtosis: " << MOMENTS[3]  << std::endl;
	outfile << MOMENTS[0] << " " << MOMENTS[1] <<  " " << MOMENTS[2]  << " " << MOMENTS[3]  << std::endl;
	#endif

	#ifdef negentropy_calc 
	outfile << "# The value of negentropy calculated from the obtained PDF (calculating integral of f(x)*logf(x) from x_min to x_max and G(x)logG(x) from mean-5*sigma to mean+5*sigma) is: " << NEGENTROPY  << "." << std::endl;
	outfile << NEGENTROPY << std::endl;	
	#endif



	outfile << "#" << "X (field value) PDF_COUNTS PDF CDF \n";                                                  //PDF=probability distribution (density) function, CDF=cumulative distribution function
	for(int i=0;i<=PDF_BINS;i++){
			outfile 
			<< PDF_x[i]
			<< " "
			<< PDF_counts[i]
			<< " "
			<< PDF[i]
			<< " "
			<< CDF[i]
			<< "\n";
		}
}
	
outfile.close();


std::cout << "The normalization check for realization with this smoothing radius yields: " << std::accumulate(PDF.begin(), PDF.end(), double(0))  << std::endl;

if(counter<NBINS-1){
			std::cout << "Increasing the radius... " << std::endl;
				}
rbin+=rstep;
std::cout << std::endl;
}


}
					//THE CASE FOR NO SMOOTHING
else{

std::cout << "Calculating the PDF..." << std::endl;
get_PDF_CDF(rbin, filtertype, field, PDF_x, PDF_counts, PDF, CDF);

std::cout << "Writing the data  to textfile..." << std::endl;
std::ofstream outfile;
outfile.open(output+"_"+mode+"_stats", std::ios::trunc | std::ios::out);

if(outfile.is_open()){

outfile << "# PDF Data generated for a simulation of size: " << BOXSIZE << "Mpc/h" << " on a grid of " << NGRID << " points" << std::endl;
outfile << "# No smoothing was chosen for the data " << std::endl;
outfile << "# Overall there are " << PDF_BINS << " bins of PDF in this file and the field values have a step dictated by mode " << mode << std::endl;

#ifdef moments_calc
std::array<double,4> MOMENTS_RAW=moments(field);
std::array<double,4> MOMENTS=moments(PDF_x);
outfile << "# Central moments of the RAW input x field: " << std::endl;
outfile << "# Mean: " << MOMENTS_RAW[0] << " " << "Deviation: " << MOMENTS_RAW[1] << " " << "Skewedness: " << MOMENTS_RAW[2] << " " << "Kurtosis: " << MOMENTS_RAW[3]  << std::endl;
outfile << MOMENTS_RAW[0] << " " << MOMENTS_RAW[1] <<  " " << MOMENTS_RAW[2]  << " " << MOMENTS_RAW[3]  << std::endl;
outfile << "# Central moments of the input x field as expressed in " << mode << " mode of operation, so for f(x) in one of the forms: " << std::endl;
outfile << "# [LINEAR_X = x, LOG_X = log(x), LOG_X_PLUS_XMIN = log(x+|x_min|), LOG_X_OVER_XAVG = log(x/x_avg)] " << std::endl;
outfile << MOMENTS[0] << " " << MOMENTS[1] <<  " " << MOMENTS[2]  << " " << MOMENTS[3]  << std::endl;
#endif

#ifdef negentropy_calc
outfile << "# The value of negentropy calculated from the obtained PDF (calculating integral of f(x)*logf(x) from x_min to x_max and G(x)logG(x) from mean-5*sigma to mean+5*sigma) is: " <<  NEGENTROPY  << "." << std::endl;
outfile << NEGENTROPY << std::endl;
#endif


outfile << "#" << "X (field value) PDF_COUNTS[X] PDF[X] CDF[X] \n";                                                //PDF=probability distribution (density) function, CDF=cumulative distribution function

for(int i=0;i<=PDF_BINS;i++){
		outfile
		<< PDF_x[i]
		<< " "
		<< PDF_counts[i]
		<< " "
		<< PDF[i]
		<< " "
		<< CDF[i]
		<< "\n";
		}
}

outfile.close();

std::cout << "The normalization check for the PDF without any smoothing yields: " << std::accumulate(PDF.begin(), PDF.end(), double(0)) << std::endl;
}


return 0;
}

/*===================================================================================================================================================
FUNCTION DEFINITIONS
===================================================================================================================================================*/


/*===================================================================================================================================================
The primary function of the program: calculating the counts, PDF and CDF of a given field 
===================================================================================================================================================*/

void get_PDF_CDF(float rbin, std::string filtertype, std::vector<float> field, 
		std::vector<float>& PDF_x, std::vector<int>& PDF_counts, std::vector<float>& PDF, std::vector<float>& CDF){

double deltax, x_max, x_min, x_avg;
float fsize;
int pdf_index;
int count_zeros=0;
int count_negative=0;

x_min=field[std::distance(field.begin(), std::min_element(field.begin(), field.end()))];
x_max=field[std::distance(field.begin(), std::max_element(field.begin(), field.end()))];
std::cout << "Before any filtering, the x_min and x_max values are equal: " <<  std::endl;
std::cout << "x_min=" << x_min << std::endl;
std::cout << "x_max=" << x_max << std::endl;
std::cout << "There are currently " << nthreads << " threads working in parellel" << std::endl;

/*===========================================================
Here we go to Fourier space for applying smoothing filters
=============================================================*/
#ifdef smoothing_filters
if( (filtertype.compare("Gaussian") == 0)   || (filtertype.compare("TopHat")==0)    ){

fftw_init_threads();

fftw_complex * field_to_fourier;
field_to_fourier = (fftw_complex*) fftw_malloc(FFTNORM*sizeof(fftw_complex));

#pragma omp parallel for num_threads(nthreads) 
for(int i=0; i<FFTNORM; i++){
	field_to_fourier[i][0]=field[i];
	field_to_fourier[i][1]=0.0;
	}

fftw_complex * Fourierfield;
Fourierfield = (fftw_complex*) fftw_malloc(FFTNORM*sizeof(fftw_complex));

fsize = rbin*(double)NGRID/(double)BOXSIZE;

std::cout << "The FFTW multi-threaded calculations are " << (fftw_init_threads() ? "working" : "NOT working")<< std::endl;
std::cout << "Creating plans for forward and backward Fast Fourier transform..." << std::endl;

fftw_plan_with_nthreads(nthreads);
fftw_plan plan_forward;
plan_forward = fftw_plan_dft_3d(NGRID, NGRID, NGRID, field_to_fourier, Fourierfield,FFTW_FORWARD, FFTW_ESTIMATE);


std::cout << "Executing the forward Fourier transform..." << std::endl;
fftw_plan_with_nthreads(1);
fftw_execute(plan_forward);
fftw_destroy_plan(plan_forward);

/*========================
Applying the filter here
========================*/


std::cout << "Applying the " << filtertype << " filter... " << std::endl;


if(filtertype.compare("Gaussian")==0){
        gaussian_f(NGRID, fsize, Fourierfield);
        }

else if(filtertype.compare("TopHat")==0){
        top_hat_f(NGRID, fsize, Fourierfield);
        }


std::cout << "Executing the backward Fourier transform..." << std::endl;

fftw_plan_with_nthreads(nthreads);
fftw_plan plan_backward;
plan_backward = fftw_plan_dft_3d(NGRID,NGRID,NGRID,Fourierfield, field_to_fourier,FFTW_BACKWARD, FFTW_ESTIMATE);

fftw_plan_with_nthreads(1);
fftw_execute(plan_backward);
fftw_destroy_plan(plan_backward);


#pragma omp parallel for num_threads(nthreads)
for(int i =0 ; i<FFTNORM; i++){
        field[i]=field_to_fourier[i][0];
        }

fftw_free(field_to_fourier);
fftw_free(Fourierfield);
fftw_cleanup_threads();


std::cout << "Rescaling the field... " << std::endl;
//RESCALE THE FIELDS (in 3D FFT transform, a factor of NGRID^3 is gained upon returning to the real space)

std::for_each(field.begin(), field.end(), [&](float& x){ x/=double(FFTNORM);} );



#ifdef ARTIFICIAL_POSITIVE
if(mode.compare("LOG_X")==0 | mode.compare("LOG_X_OVER_XAVG")==0){ 
std::cout << "Removing potential negative values... " << std::endl;
auto it = remove_if(field.begin(), field.end(),  [](float x) {return x <= 0.0000000000; });
field.erase(it, field.end());
count_negative=FFTNORM-field.size();
std::cout << count_negative << " non-positive values have been detected" << std::endl;
FFTNORM-=count_negative;
}
#endif

}
#endif //endif for smoothing filters...
/*===========================================================
Getting the x_min, x_max of the fields 
=============================================================*/
std::cout << "Getting the x_min and x_max..." << std::endl;
x_min=field[std::distance(field.begin(), std::min_element(field.begin(), field.end()))];
x_max=field[std::distance(field.begin(), std::max_element(field.begin(), field.end()))];

if(x_min <0 & ( mode.compare("LOG_X")==0 || mode.compare("LOG_X_OVER_XAVG")==0 )  ){
	std::cout << "The minimal field value is negative and you chose the log(x) or log(x/x_avg) presentation - it is not possible!" << std::endl;
        throw std::logic_error("The minimal field value is negative and you chose the log(x) or log(x/x_avg) presentation - it is not possible!");
	}
if(x_min==0){
        std::cout << "Careful, x_min is zero!" << " x_min = " << x_min << std::endl;
        }
else if(x_min<0){
        std::cout << "Careful, x_min is negative! It equals " << x_min <<  std::endl;
        }
else{
 	std::cout << "x_min is equal to: " << x_min << std::endl;
	}
std::cout << "Meanwhile, x_max is equal to: " << x_max << std::endl;

if(x_min>0 & mode.compare("LOG_X_PLUS_XMIN")==0 ){
		std::cout << "x_min is positive, yet " << mode << "has been chosen - is it intentional?" << std::endl;
		}

/*===========================================================
Here we perform field-specific calculations, for example,
obtatining the density contrast from density, rescaling
the velocity divergence or otherwise any other operation
=============================================================*/
std::cout << "Performing field-specific operations... " << std::endl;

if(mode.compare("LOG_X_OVER_XAVG")==0){
	std::cout << "Calculating the quantity log(x/x_xavg)" << std::endl;
	x_avg=std::accumulate(field.begin(), field.end(), double(0));
	x_avg/=double(FFTNORM); //now x_avg is an average over all the positive values of the field! 
	std::cout << "The average value of the field is equal to: " << x_avg << std::endl;
	std::for_each(field.begin(), field.end(), [&](float& x){ x=x/x_avg;} ); //obtaining the x/x_avg, for example ---> delta+1 :) 
	}

else if(mode.compare("LOG_X_PLUS_XMIN")==0){

	if(x_min<0){x_min*=-1;}
	std::cout << "Taking the modulus of x_min gives: " << x_min << std::endl;
	std::cout << "Calculating the PDF for quantity log(x+|x_xmin|)" << std::endl;
	std::for_each(field.begin(), field.end(), [&](float& x){x=x+x_min;}); //obtaining x+|x_min|
	#ifdef ARTIFICIAL_POSITIVE
        std::cout << "Removing the zero values from the field..." << std::endl;
        auto it = remove_if(field.begin(), field.end(),  [](float x) {return x == 0.0; });
        field.erase(it, field.end());
        count_zeros=FFTNORM-field.size();
        std::cout << count_zeros << " zeroes have been detected" << std::endl;
        FFTNORM-=count_zeros;
        #endif  
	}


std::cout << "Getting the x_min and x_max once again, after field-specific calculations..." << std::endl;
x_min=field[std::distance(field.begin(), std::min_element(field.begin(), field.end()))];
x_max=field[std::distance(field.begin(), std::max_element(field.begin(), field.end()))];
std::cout << "x_min=" << x_min << std::endl;
std::cout << "x_max=" << x_max << std::endl;

/*==========================================================
Calculate the negentropy here...
===========================================================*/
#ifdef negentropy_calc
std::cout << "Calculating the negentropy... " << std::endl;
NEGENTROPY = negentropy(field);
std::cout << "Negentropy equal to: " << NEGENTROPY << std::endl;
#endif
/*===========================================================
Here we calculate the delta x in the desired representation 
of the field 
if we consider ABS_MIN and ABS_MAX, 
or we take these values from the field itself
=============================================================*/

if(mode.compare("LOG_X")==0){

		deltax=log10(x_max/x_min)/double(PDF_BINS);

		#ifdef FIXED_ABS
			deltax=log10(ABS_MAX/ABS_MIN)/double(PDF_BINS);
		#endif
	}

else if(mode.compare("LINEAR_X")==0){
		deltax=(x_max-x_min)/double(PDF_BINS);
		#ifdef FIXED_ABS
                        deltax=(ABS_MAX-ABS_MIN)/double(PDF_BINS);
                #endif
	}

else if(mode.compare("LOG_X_PLUS_XMIN")==0){
		
		deltax=log10(x_max/x_min)/double(PDF_BINS);

		#ifdef FIXED_ABS
                        deltax=log10(ABS_MAX/ABS_EPSILON))/double(PDF_BINS);
                #endif

	}

else if(mode.compare("LOG_X_OVER_XAVG")==0){
		deltax=log10(x_max/x_min)/double(PDF_BINS);
		
		#ifdef FIXED_ABS
                        deltax=log10(ABS_MAX/ABS_MIN)/double(PDF_BINS);
                #endif
	}

else{
	std::cout << "Unrecognised mode of operation! This should have been prevented at the start of the program..." << std::endl;
	}

std::cout << "With the number of given bins, mode of operation:" << mode << " and ABS_MIN, ABS_MAX values, the deltax is equal to: " << deltax << std::endl;

/*==========================================
Getting the PDF_counts filled up 
===========================================*/
std::cout << "Filling up the PDF array... " << std::endl;

if(mode.compare("LOG_X")==0){

                #ifndef FIXED_ABS

		#pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
		for(int i=0 ; i < FFTNORM ; i++){
	        	pdf_index=0;
			if(field[i]>=x_min & field[i]>0){
				pdf_index= log10(field[i]/x_min)/deltax ;
	                        PDF_counts[pdf_index]+=1;
				}
			}
		#endif

		#ifdef FIXED_ABS
		#pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
                for(int i=0 ; i < FFTNORM ; i++){
                        if(field[i]>=x_min & field[i]>0){
		                pdf_index= log10(field[i]/ABS_EPSILON)/deltax ;
	                        PDF_counts[pdf_index]+=1;
				}
                        }
		#endif
        }

else if(mode.compare("LINEAR_X")==0){

                #ifndef FIXED_ABS
	       	#pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
		for(int i=0 ; i < FFTNORM ; i++){
                        pdf_index= (field[i]-x_min)/deltax ;
                        PDF_counts[pdf_index]+=1;
                        }
		#endif

                #ifdef FIXED_ABS
		#pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
                for(int i=0 ; i < FFTNORM ; i++){
                        pdf_index= (field[i]-ABS_MIN)/deltax ;
                        PDF_counts[pdf_index]+=1;
                        }

                #endif
        }

else if(mode.compare("LOG_X_OVER_XAVG")==0){

                #ifndef FIXED_ABS
		if(x_min==0){
                        std::cout << "x_min turns out to be 0... changing the x_min to ABS_EPSILON" << std::endl;
                        x_min=ABS_EPSILON;}
		

		#pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
                for(int i=0 ; i < FFTNORM ; i++){
                        pdf_index= log10(field[i]/x_min)/deltax ;
                        PDF_counts[pdf_index]+=1;
                }
		#endif

		#ifdef FIXED_ABS
                #pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
                for(int i=0 ; i < FFTNORM ; i++){
                        pdf_index= log10(field[i]/ABS_EPSILON)/deltax ;
                        PDF_counts[pdf_index]+=1;
                        }

                #endif
}

else if(mode.compare("LOG_X_PLUS_XMIN")==0){

                #ifndef FIXED_ABS
                #pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
                for(int i=0 ; i < FFTNORM ; i++){
			pdf_index=0;
			if(field[i]>0){
			pdf_index=log10( field[i]/x_min)/ deltax;
			}
                        PDF_counts[pdf_index]+=1;
                }
		#endif

                #ifdef FIXED_ABS		//the minimum, being ABS_EPSILON is set by hand...
                #pragma omp parallel for num_threads(nthreads) private(pdf_index)  reduction(vec_int_plus: PDF_counts)
                for(int i=0 ; i < FFTNORM ; i++){
			pdf_index=0;
			if(field[i]>=ABS_EPSILON){
		                        pdf_index= log10( field[i]/ABS_EPSILON )/ deltax;			
				}
                        PDF_counts[pdf_index]+=1;
                        }

                #endif
}


else{
        std::cout << "Unrecognised mode of operation! This should have been prevented at the start of the program..." << std::endl;
        }



/*=====================================================
Filling the PDF_x array 
=====================================================*/
if(mode.compare("LOG_X")==0){

		for(int i = 0; i<=PDF_BINS; i++){
			PDF_x[i]=x_min*pow(10,i*deltax);			
		}

		#ifdef FIXED_ABS
		for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=ABS_EPSILON*pow(10,i*deltax);
                }
		#endif
}

else if(mode.compare("LINEAR_X")==0){

                for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=x_min+i*deltax;
                }

                #ifdef FIXED_ABS
                for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=ABS_EPSILON+i*deltax;
                }
                #endif
}

else if(mode.compare("LOG_X_OVER_XAVG")==0){

                for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=x_min*pow(10,i*deltax);
                }

                #ifdef FIXED_ABS
                for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=ABS_EPSILON*pow(10,i*deltax);
                }
                #endif
}

else if(mode.compare("LOG_X_PLUS_XMIN")==0){
		
		for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=x_min*pow(10,i*deltax);
			}

                #ifdef FIXED_ABS
                for(int i = 0; i<=PDF_BINS; i++){
                        PDF_x[i]=ABS_EPSILON*pow(10,i*deltax);
                }
                #endif
}

else{
        std::cout << "Unrecognised mode of operation! This should have been prevented at the start of the program..." << std::endl;
	}


/*=====================================================
Normalizing the PDF_counts array, thus obtaining 
PDF and CDF arrays
=====================================================*/

std::copy_n(PDF_counts.begin(), PDF_BINS+1, PDF.begin());
std::for_each(PDF.begin(), PDF.end(), [&](float& x){ x/=double(FFTNORM); }); //normalizing the probability distribution function

float temp=0;
for(int i=0; i<=PDF_BINS; i++){
        temp+=PDF[i];
        CDF[i]=temp;
	}

std::cout << "Succesfully ended the function get_PDF_CDF" << std::endl;
}



/*===================================================================================================================================================
Implementation of smoothing kernels in Fourier space: 
- Gaussian
- Top Hat 
===================================================================================================================================================*/
#ifdef smoothing_filters

/*==============================
GAUSSIAN FILTER FOR K-SPACE
==============================*/

void gaussian_f(int NGRID, float fsize, fftw_complex *field)
{

   int ni,nf,icx,icy,icz,isx,isy,isz;
   float irz,iry,irx, zk, zk2, rks, tmp =0.0, cfilter=0.0;
   rks = fsize*2.0*M_PI/(float)NGRID;

   ni = 1;
   nf = NGRID/2;

#pragma omp parallel for collapse(3) private(isz,irz,isy,iry,isx,irx,zk2,tmp) num_threads(nthreads)
  for(icz=ni; icz<=nf; icz++)
   {
    for(icy=ni; icy<=nf; icy++)
     {
      for(icx=ni; icx<=nf; icx++)
       {

        isy = NGRID - icy +1;
        iry = (float)icy - 0.5;
        isz = NGRID - icz +1;
        irz = (float)icz - 0.5;
        isx = NGRID - icx +1;
        irx = (float)icx - 0.5;
        zk2 = (irx*irx + iry*iry + irz*irz) * rks * rks;
        tmp = exp(-0.5*zk2);

        field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
        field[LOC(icx-1,isy-1,icz-1)][0] *= tmp; field[LOC(icx-1,isy-1,icz-1)][1] *= tmp;
        field[LOC(icx-1,icy-1,isz-1)][0] *= tmp; field[LOC(icx-1,icy-1,isz-1)][1] *= tmp;
        field[LOC(icx-1,isy-1,isz-1)][0] *= tmp; field[LOC(icx-1,isy-1,isz-1)][1] *= tmp;
        field[LOC(isx-1,icy-1,icz-1)][0] *= tmp; field[LOC(isx-1,icy-1,icz-1)][1] *= tmp;
        field[LOC(isx-1,isy-1,icz-1)][0] *= tmp; field[LOC(isx-1,isy-1,icz-1)][1] *= tmp;
        field[LOC(isx-1,icy-1,isz-1)][0] *= tmp; field[LOC(isx-1,icy-1,isz-1)][1] *= tmp;
        field[LOC(isx-1,isy-1,isz-1)][0] *= tmp; field[LOC(isx-1,isy-1,isz-1)][1] *= tmp;
       }
     }
   }
}



/*==============================
TOP HAT FILTER FOR K-SPACE
==============================*/

void top_hat_f(int NGRID, float fsize, fftw_complex *field)
{
   int ni,nf,icx,icy,icz,isx,isy,isz;
   float irz,iry,irx, zk, rks, tmp, cfilter=0.0;

    field[LOC(0,0,0)][0] *= 1.0; field[LOC(0,0,0)][1] *= 1.0;
    rks = fsize*2.0*M_PI/(float)NGRID;
    ni = 2;
    nf = NGRID/2 +1;

// os x
   icz =1;
   icy =1;

#pragma omp parallel for private(isx,zk,tmp) num_threads(nthreads)
   for(icx=ni; icx<=nf; icx++)
    {
     isx = NGRID - icx +2;
     zk = ((float)icx-1) * rks;
     tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

     field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
     field[LOC(isx-1,icy-1,icz-1)][0] *= tmp; field[LOC(isx-1,icy-1,icz-1)][1] *= tmp;
    }

// os y
   icz =1;
   icx =1;

#pragma omp parallel for private(isy,zk,tmp) num_threads(nthreads)                                                                                             
   for(icy=ni; icy<=nf; icy++)
    {
     isy = NGRID - icy +2;
     zk = ((float)icy-1) * rks;
     tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

     field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
     field[LOC(icx-1,isy-1,icz-1)][0] *= tmp; field[LOC(icx-1,isy-1,icz-1)][1] *= tmp;
    }

// os z
   icx =1;
   icy =1;

#pragma omp parallel for  private(isz,zk,tmp) num_threads(nthreads)                                                                                             
   for(icz=ni; icz<=nf; icz++)
    {
     isz = NGRID - icz +2;
     zk = ((float)icz-1) * rks;
     tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

     field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
     field[LOC(icx-1,icy-1,isz-1)][0] *= tmp; field[LOC(icx-1,icy-1,isz-1)][1] *= tmp;
    }

// k_z =0 plane
  icz = 1;

#pragma omp parallel for collapse(2) private(isy,iry,isx,irx,zk,tmp) num_threads(nthreads)                                                                                             
  for(icy = ni; icy <=nf; icy++)
   {
    for(icx=ni; icx<=nf; icx++)
     {
      isy = NGRID- icy +2;
      iry = pow((float)icy-1.0,2.0);
      isx = NGRID - icx +2;
      irx = pow((float)icx-1.0,2.0);
      zk = sqrt(irx+iry) * rks;
      tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

      field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
      field[LOC(icx-1,isy-1,icz-1)][0] *= tmp; field[LOC(icx-1,isy-1,icz-1)][1] *= tmp;
      field[LOC(isx-1,icy-1,icz-1)][0] *= tmp; field[LOC(isx-1,icy-1,icz-1)][1] *= tmp;
      field[LOC(isx-1,isy-1,icz-1)][0] *= tmp; field[LOC(isx-1,isy-1,icz-1)][1] *= tmp;
     }
   }

// k_y=0 plane
  icy = 1;

#pragma omp parallel for collapse(2) private(isx, irx, isz, irz, zk, tmp) num_threads(nthreads)                                                                                             
  for(icx = ni; icx <=nf; icx++)
   {
    for(icz=ni; icz<=nf; icz++)
     {

      isx = NGRID- icx +2;
      irx = pow((float)icx-1.0,2.0);
      isz = NGRID - icz +2;
      irz = pow((float)icz-1.0,2.0);
      zk = sqrt(irx+irz) * rks;
      tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

      field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
      field[LOC(icx-1,icy-1,isz-1)][0] *= tmp; field[LOC(icx-1,icy-1,isz-1)][1] *= tmp;
      field[LOC(isx-1,icy-1,icz-1)][0] *= tmp; field[LOC(isx-1,icy-1,icz-1)][1] *= tmp;
      field[LOC(isx-1,icy-1,isz-1)][0] *= tmp; field[LOC(isx-1,icy-1,isz-1)][1] *= tmp;
     }
   }

// k_x=0 plane
  icx = 1;

#pragma omp parallel for collapse(2) private(isy,iry,isz,irz,zk,tmp) num_threads(nthreads)                                                                                             
  for(icy = ni; icy <=nf; icy++)
   {
    for(icz=ni; icz<=nf; icz++)
     {

      isy = NGRID- icy +2;
      iry = pow((float)icy-1.0,2.0);
      isz = NGRID - icz +2;
      irz = pow((float)icz-1.0,2.0);
      zk = sqrt(iry+irz) * rks;
      tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

      field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
      field[LOC(icx-1,isy-1,icz-1)][0] *= tmp; field[LOC(icx-1,isy-1,icz-1)][1] *= tmp;
      field[LOC(icx-1,icy-1,isz-1)][0] *= tmp; field[LOC(icx-1,icy-1,isz-1)][1] *= tmp;
      field[LOC(icx-1,isy-1,isz-1)][0] *= tmp; field[LOC(icx-1,isy-1,isz-1)][1] *= tmp;
     }
   }

// the rest is volume-symetric

#pragma omp parallel for collapse(3) private(isz,irz,isy,iry,isx,irx,zk,tmp) num_threads(nthreads)                                                                                             
  for(icz=ni; icz<=nf; icz++)
   {
    for(icy=ni; icy<=nf; icy++)
     {
      for(icx=ni; icx<=nf; icx++)
       {

	isy = NGRID- icy +2;
        iry = pow((float)icy-1.0,2.0);
	isz = NGRID - icz +2;
    	irz = pow((float)icz-1.0,2.0);
        isx = NGRID - icx +2;
        irx = pow((float)icx -1,2.0);
        zk=sqrt(irx+iry+irz) * rks;
        tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);
  
        field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
        field[LOC(icx-1,isy-1,icz-1)][0] *= tmp; field[LOC(icx-1,isy-1,icz-1)][1] *= tmp;
        field[LOC(icx-1,icy-1,isz-1)][0] *= tmp; field[LOC(icx-1,icy-1,isz-1)][1] *= tmp;
        field[LOC(icx-1,isy-1,isz-1)][0] *= tmp; field[LOC(icx-1,isy-1,isz-1)][1] *= tmp;
        field[LOC(isx-1,icy-1,icz-1)][0] *= tmp; field[LOC(isx-1,icy-1,icz-1)][1] *= tmp;
        field[LOC(isx-1,isy-1,icz-1)][0] *= tmp; field[LOC(isx-1,isy-1,icz-1)][1] *= tmp;
        field[LOC(isx-1,icy-1,isz-1)][0] *= tmp; field[LOC(isx-1,icy-1,isz-1)][1] *= tmp;
        field[LOC(isx-1,isy-1,isz-1)][0] *= tmp; field[LOC(isx-1,isy-1,isz-1)][1] *= tmp;
       }
     }
    }

}
#endif


/*================================================================================================================================================
Additional functionalities:
        -calculation of moments
        -calculation of negentropy
===================================================================================================================================================*/


/*===================================================================================================================================================
Moments calculator: 
===================================================================================================================================================*/
std::array<double,4> moments(const std::vector<float>& field){

std::array<double,4> moments_array;

double mean=std::accumulate(field.begin(), field.end(), double(0));
mean/=field.size();

moments_array[0]=mean;
double temp;

for(int i=0; i<field.size(); i++){

	temp=pow(field[i]-mean,2);
	moments_array[1]+=temp;

	temp*=(field[i]-mean);
	moments_array[2]+=temp;

	temp*=(field[i]-mean);
	moments_array[3]+=temp;
	
	}

moments_array[1]/=field.size();
moments_array[1]=sqrt(moments_array[1]);

moments_array[2]/=field.size();
moments_array[2]/=pow(moments_array[1],3);

moments_array[3]/=field.size();
moments_array[3]/=pow(moments_array[1],4);


return moments_array;
}

/*===================================================================================================================================================
Negentropy calculator
===================================================================================================================================================*/

#ifdef negentropy_calc

#define PRECISION 1e-4
#define NUMPOINTS 10000
#define EPS 1e-9

float negentropy( std::vector<float> field){

double integral_f=0;
double integral_g=0;
double abserr;


std::array<double,4> MOMENTS=moments(field);
float mean, variance, deviation;

mean=MOMENTS[0];
deviation=MOMENTS[1];
variance=MOMENTS[1]*MOMENTS[1];

double x_min=field[std::distance(field.begin(), std::min_element(field.begin(), field.end()))];
double x_max=field[std::distance(field.begin(), std::max_element(field.begin(), field.end()))];
double deltax;

deltax=(x_max-x_min)/PDF_BINS;


double* pdf_val;
double* pdf_arg;
pdf_val = (double*) calloc((PDF_BINS+1),sizeof(double));
pdf_arg = (double*) calloc((PDF_BINS+1),sizeof(double));

int pdf_index;

for(int i=0 ; i < FFTNORM ; i++){
                        pdf_index= (field[i]-x_min)/deltax ;
                        pdf_val[pdf_index]+=1;
	}

for(int i=0; i<= PDF_BINS; i++){
	pdf_val[i]/=(double)FFTNORM;
	}


for(int i=0; i<=PDF_BINS; i++){
                        pdf_arg[i]=x_min+i*deltax;
      }


double integral_f_simple=0;
for(int i=0; i<PDF_BINS; i++){
	if(pdf_val[i]!=0){
		integral_f_simple+= pdf_val[i]*log10(pdf_val[i])*(pdf_arg[i+1]-pdf_arg[i]);
		}
}
std::cout << "Integral_f_simple is: " << integral_f_simple << std::endl;		


//gsl_set_error_handler_off();
acc_pk = gsl_interp_accel_alloc ();
pkl_spline = gsl_spline_alloc (gsl_interp_linear,PDF_BINS+1);
gsl_spline_init (pkl_spline, pdf_arg, pdf_val, PDF_BINS+1);
std::cout  << "Initialized splines" << std::endl;


gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUMPOINTS);

double a;
double b;
a=x_min;
b=x_max;



gsl_function f;
f.params = NULL;
f.function = &integrand_f;
std::cout  << "Integrating the f(x) function from " << a << " to " << b << "..." << std::endl;



if (gsl_integration_qag(&f, a, b, 0, 100*PRECISION, NUMPOINTS, GSL_INTEG_GAUSS61, w, &integral_f, &abserr) != GSL_SUCCESS){
			std::cout << "Something went wrong" << std::endl;
		}

std::cout << "Integral_f = " << integral_f << std::endl;

// testing on a binned gaussian
/*
std::cout << "TESTING: Testing on a binned Gaussian"<<std::endl;
a=mean-5*deviation;
b=mean+5*deviation;
deltax=(b-a)/PDF_BINS;

for(int i=0; i<=PDF_BINS; i++){
                        pdf_arg[i]=a+i*deltax;
      }


for(int i=0 ; i <=  PDF_BINS ; i++){
                        pdf_val[i]= gsl_ran_gaussian_pdf(pdf_arg[i]-mean, deviation) ;
        }

integral_f_simple=0;
for(int i=0; i<PDF_BINS; i++){
        if(pdf_val[i]!=0){
                integral_f_simple+= pdf_val[i]*log10(pdf_val[i])*(pdf_arg[i+1]-pdf_arg[i]);
                }
}



std::cout << "TESTING: Integral_f = " << integral_f_simple << std::endl;
*/
//end of testing 

/*
#ifdef NEGENTROPY_LOG
std::for_each(field.begin(),field.end(), [&](float& x){x=log10(x);});
MOMENTS=moments(field);
mean=MOMENTS[0];
deviation=MOMENTS[1];
variance=MOMENTS[1]*MOMENTS[1];
std::cout << "Parameters for the Gaussian are taken from log(X) distribution, mean = " << mean << " , deviation = " << deviation << "."  << std::endl;
#endif 
*/

a=mean-5*deviation;
b=mean+5*deviation;

std::cout  << "Integrating the g(x) function from mean-5*sigma to mean+5*sigma, i.e.: [" << a<< "," << b << "]." << std::endl;

gsl_function g;
struct GAUSS_PARAMS gauss_params;
gauss_params={mean,deviation};

g.params= &gauss_params;
g.function = &gaussian;

//DEBUGGING: save the corresponding gaussian to file 
/*
std::ofstream outg;
outg.open(output+"_gaussian_"+filtertype+"_"+mode, std::ios::trunc | std::ios::out);
double gdelta=(b-a)/100;
double xtemp;
for(int i=0; i<100; i++){
xtemp=a+i*gdelta;
outg << xtemp << " " << gsl_ran_gaussian_pdf(xtemp-mean, deviation) << "\n"; 
}
outg.close();*/
//END OF DEBUGGING

if (gsl_integration_qag(&g, a, b, 0, PRECISION, NUMPOINTS, GSL_INTEG_GAUSS21, w, &integral_g, &abserr) != GSL_SUCCESS){
                        std::cout << "Something went wrong" << std::endl;
                }

std::cout << "Integral_g = " << integral_g << std::endl;

std::cout  << "Freeing space of GSL routines..." << std::endl;

gsl_spline_free(pkl_spline);
gsl_interp_accel_free(acc_pk);
gsl_integration_workspace_free(w);
free(pdf_val);
free(pdf_arg);

std::cout  << "Done calculating negentropy!" << std::endl;


return integral_f-integral_g;
}


static double pdfspline(double x, void * p){
double result = gsl_spline_eval (pkl_spline, x, acc_pk);
return result;
}

static double integrand_f(double x, void* p){
if(pdfspline(x,NULL)==0){
	return 0;
	}
else
	return pdfspline(x,NULL)*log10(pdfspline(x,NULL));
}


static double gaussian(double x, void * p){
struct GAUSS_PARAMS * params = (struct GAUSS_PARAMS *)p;
double mean = (params->MEAN);
double sigma = (params->SIGMA);

double result = gsl_ran_gaussian_pdf(x-mean, sigma);
return result*log10(result);
}
#endif


