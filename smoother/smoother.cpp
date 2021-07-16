/* This simple program serves to smooth out a field using one of the two filters on a given scale, described by the variable RADIUS
   
   Compilation: 
   USE the PrgEnv-gnu (it provides the GNU C++ compiler version 8.X, necessary for the filesystem library ---> -lstdc++fs flag)

   g++ -std=c++17 smoother.cpp filters.cpp -I /lustre/tetyda/home/topolski/include -L /lustre/tetyda/home/topolski/lib -lfftw3_omp -lfftw3 -lm -fopenmp -lstdc++fs -o smoother

   Usage: 

   ./smoother BINARY_FILE OUTPUT NGRID BOXSIZE FilterConvolveName TypeConvolve RadiusConvolve FilterDeconvolveName TypeDeconvolve RadiusDeconvolve 

   Example:

   ./smoother /home/topolski/dtfe_files/512x500Mpc/2/512x500Mpc_R02_W7_LCDM_z1.130.DTFE512.a_den /home/topolski/CCG_CODES/smoother/testing 512 500 Gaussian Radial 3 NGP Cube 5

   Filter types:
	- TopHat
	- Gaussian


MOST RECENT UPDATE:	22.05.2021
AUTHOR:			Konrad Topolski
@:			k.topolski2@student.uw.edu.pl

Implementation of Gaussian, TopHat smoothing filters coutresy of:
Pawel Ciecielag
 
																       */

/*===================================================================================================================================================
HEADERS
===================================================================================================================================================*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <cmath>
#include <fftw3.h>
#include <numeric> 
#include <cstring>
#include "filters.h"

/*===================================================================================================================================================
GLOBAL VARIABLES
===================================================================================================================================================*/
std::string binary_file, output;
std::string FilterConvolveName = "None";
std::string TypeConvolve = "Cube";
std::string FilterDeconvolveName = "None";
std::string TypeDeconvolve = "Radial";
std::string RadiusConvolve = "5";
std::string RadiusDeconvolve = "2";
double radConv = std::stod(RadiusConvolve);
double radDeconv = std::stod(RadiusDeconvolve);
long int FFTNORM;

/*===================================================================================================================================================
THE MAIN FUNCTION 
===================================================================================================================================================*/

int main(int argc, char* argv[]){
 

if(argc<8){ //arguments for convolving are\ necessary, deconvolving - it's not 
	std::cout <<"Not enough arguments supplied!" << std::endl;
	std::exit(1);
        }


else{
   binary_file = argv[1];
   output= argv[2];
   NGRID    = atoi(argv[3]);
   BOXSIZE = atof(argv[4]);
   FFTNORM = (long int)pow(NGRID,3);
   FilterConvolveName = argv[5];
   TypeConvolve	 = argv[6];
   RadiusConvolve = argv[7];
   radConv = std::stod(RadiusConvolve);}


if(argc>8){
   FilterDeconvolveName = argv[8];
   TypeDeconvolve  = argv[9];
   RadiusDeconvolve = argv[10];
   radDeconv = std::stod(RadiusDeconvolve);
	}





if(!std::filesystem::exists(binary_file)){
	throw std::invalid_argument("Invalid name of the input file! Check the name!");
	}
	

if(radConv<0){
	throw std::range_error("The value of smoothing convolution radius is incorrect!");
	}
if(radDeconv<0){
        throw std::range_error("The value of smoothing convolution radius is incorrect!");
        }


if(FilterConvolveName.compare("Gaussian") != 0 & FilterConvolveName.compare("TopHat") != 0 & FilterConvolveName.compare("None")!=0 & FilterConvolveName.compare("NGP")!=0 & 
			FilterConvolveName.compare("CIC")!=0 & FilterConvolveName.compare("TSC")!=0){
	throw std::invalid_argument("Unrecognised convolving filter name! It must be one of the following: \n TopHat \n Gaussian \n NGP \n CIC \n TSC \n None"); 
	}

if(FilterDeconvolveName.compare("Gaussian") != 0 & FilterDeconvolveName.compare("TopHat") != 0 & FilterDeconvolveName.compare("None")!=0 & FilterDeconvolveName.compare("NGP")!=0 & 
                        FilterDeconvolveName.compare("CIC")!=0 & FilterDeconvolveName.compare("TSC")!=0){
        throw std::invalid_argument("Unrecognised deconvolving filter name! It must be one of the following: \n TopHat \n Gaussian \n NGP \n CIC \n TSC \n None"); 
        }


std::cout << "The binary file name and location is: " << binary_file << std::endl;
std::cout << "The output name and location: " << output << std::endl;
std::cout << "Grid size: " << NGRID << std::endl;
std::cout << "Box size: " << BOXSIZE << std::endl;
std::cout << "The smoothing convolution kernel is: " << FilterConvolveName << " of " << TypeConvolve << " type." << std::endl;
std::cout << "The smoothing convolution radius is: " << RadiusConvolve << " Mpc/h" << std::endl;

std::cout << "The smoothing convolution kernel is: " << FilterDeconvolveName << " of " << TypeDeconvolve << " type." << std::endl;
std::cout << "The smoothing deconvolution radius is: " << RadiusDeconvolve << " Mpc/h" << std::endl;



fftw_complex* field;
field = (fftw_complex*) fftw_malloc(FFTNORM * sizeof(fftw_complex));


std::ifstream fin(binary_file,std::ios::binary);
float dummy;
for(int i=0; i<FFTNORM; i++){
fin.read((char*)&dummy, sizeof(float));
field[i][0]=dummy;
field[i][1]=0.0;
}

fin.close();


std::cout << "The data loaded is at least: " << (FFTNORM*sizeof(float))/pow(1024,3) << "GB in size" << std::endl;

int count_negative;

std::cout << "There are currently " << nthreads << " threads working in parellel" << std::endl;



/*===========================================================
Here we go to Fourier space for applying smoothing filters
=============================================================*/


fftw_init_threads();

fftw_complex * Fourierfield;
Fourierfield = (fftw_complex*) fftw_malloc(FFTNORM*sizeof(fftw_complex));


std::cout << "The FFTW multi-threaded calculations are " << (fftw_init_threads() ? "working" : "NOT working")<< std::endl;
std::cout << "The number of threads used for FFT is: " << nthreads << std::endl;
std::cout << "Creating plans for forward and backward Fast Fourier transform..." << std::endl;

fftw_plan_with_nthreads(nthreads);
fftw_plan plan_forward;
plan_forward = fftw_plan_dft_3d(NGRID, NGRID, NGRID, field, Fourierfield,FFTW_FORWARD, FFTW_ESTIMATE);


std::cout << "Executing the forward Fourier transform..." << std::endl;
fftw_plan_with_nthreads(1);
fftw_execute(plan_forward);
fftw_destroy_plan(plan_forward);

/*========================
Applying the filter here
========================*/


std::cout << "Applying the filters..."  << std::endl;

FilterInterface(field, FilterConvolveName,FilterDeconvolveName, TypeConvolve, TypeDeconvolve, radConv, radDeconv);

std::cout << "Executing the backward Fourier transform..." << std::endl;

fftw_plan_with_nthreads(nthreads);
fftw_plan plan_backward;
plan_backward = fftw_plan_dft_3d(NGRID,NGRID,NGRID,Fourierfield, field,FFTW_BACKWARD, FFTW_ESTIMATE);

fftw_plan_with_nthreads(1);
fftw_execute(plan_backward);
fftw_destroy_plan(plan_backward);

std::cout << "Rescaling the field... " << std::endl;
//RESCALE THE FIELDS (in 3D FFT transform, a factor of NGRID^3 is gained upon returning to the real space)

std::for_each(field, field+FFTNORM, [&](fftw_complex& x){ x[0]/=float(FFTNORM);} );


std::ofstream outfile;

std::cout << "Writing the output to file..." << std::endl;

std::string OUTPUTNAME;
OUTPUTNAME=output;
if(FilterConvolveName.compare("None")!=0){
OUTPUTNAME= OUTPUTNAME+"_CONV_"+FilterConvolveName+"_"+TypeConvolve+"_R"+RadiusConvolve;
}
if(FilterDeconvolveName.compare("None")!=0){
OUTPUTNAME= OUTPUTNAME+"_DECONV_"+FilterDeconvolveName+"_"+TypeDeconvolve+"_R"+RadiusDeconvolve;
}

std::ofstream fout(OUTPUTNAME,std::ios::out | std::ios::binary);

float* savefield;
savefield = (float*) fftw_malloc(FFTNORM * sizeof(float));


#pragma omp parallel for num_threads(nthreads)
for(int i =0 ; i<FFTNORM; i++){
	savefield[i]=field[i][0];
	}


fout.write((char*)&savefield[0], FFTNORM*sizeof(float));
fftw_free(savefield);

fout.close();


fftw_cleanup_threads();
fftw_free(Fourierfield);


fftw_free(field);

return 0;

}

