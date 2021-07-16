/* This simple program subsamples the data in a binary file, artificially creating a data suite with smaller resolution (e.x. 1024 ---> 512), 
   but with the same box length as the original
   
   Compilation: 
   USE THE PrgEnv-gnu (with the GNU C++ compiler version 8.X)

   g++ -std=c++17 subsampler.cpp -I /lustre/tetyda/home/topolski/include -L /lustre/tetyda/home/topolski/lib -lm -fopenmp -lstdc++fs -o subsampler

   Usage: 

   ./subsampler.exe binary_file output_pathname NGRID_IN NGRID_OUT NCOMP


   Example:
./subsampler /home/topolski/dtfe_files/512x500Mpc/2/512x500Mpc_R02_W7_LCDM_z1.130.DTFE512.a_den /home/topolski/CCG_CODES/subsampler/SAMPLE256 512 256 1


MOST RECENT UPDATE:	22.05.2021
AUTHOR:			Konrad Topolski
@:			k.topolski2@student.uw.edu.pl
 
																       */

#define PARALLEL // Use parallel or otherwise, serial calculations; NOTE: parallel calculations are visibly faster,
		 // but with large sets of data (e.x. velShear for NGRID=1024), there is too much memory consumption
		 // This is possibly associated with some memory overheads done by OpenMP in order to introduce the parallelism
		 // In conclusion - unless your data has >= 5 components and NGRID>=1024, 
		 // you should be fine using this program with OpenMP support.
/*===================================================================================================================================================
Functions for picking out the element location
===================================================================================================================================================*/
long int LOC(int x_i, int y_j, int z_k, int NGRID);
long int LOC_comp(long int LOC, int ncomp);
/*===================================================================================================================================================
HEADERS
===================================================================================================================================================*/
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <cmath>

#include <numeric> 

#ifdef PARALLEL
#include <omp.h>
#pragma omp declare reduction(vec_plus_float : \
		    std::vector<float> : \
                    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))      
#endif     

/*===================================================================================================================================================
GLOBAL VARIABLES
===================================================================================================================================================*/
std::string binary_file, output;
int NGRID_IN, NGRID_OUT;
long int FFTNORM_IN, FFTNORM_OUT;
int NCOMP;
#ifdef PARALLEL
int nthreads=omp_get_max_threads();
#endif
/*===================================================================================================================================================
THE MAIN FUNCTION 
===================================================================================================================================================*/

int main(int argc, char* argv[]){
 
if(argc<6){
	throw std::invalid_argument("Not enough arguments supplied!");
        }


   binary_file = argv[1];
   output= argv[2];
   NGRID_IN    = atoi(argv[3]);
   NGRID_OUT = atoi(argv[4]);
   NCOMP=atoi(argv[5]);
   FFTNORM_IN = pow(NGRID_IN,3);
   FFTNORM_OUT= pow(NGRID_OUT,3);


if(!std::filesystem::exists(binary_file)){
	throw std::invalid_argument("Invalid name of the input file! Check the name!");
	}
	
std::cout << "The binary file name and location is: " << binary_file << std::endl;
std::cout << "The output name and location: " << output << std::endl;
std::cout << "Grid size for the input: " << NGRID_IN << std::endl;
std::cout << "Grid size for the output: " << NGRID_OUT << std::endl;
std::cout << "Number of components for the subsampled field: " << NCOMP << std::endl;

int num_divisions=0;
int temp=NGRID_IN/NGRID_OUT;
for(int i=0; i<8; i++){
	if(temp==pow(2,i)){
		num_divisions=i;
	}
}

if(num_divisions==0){
	std::cout << "The program either could not find the number of divisions needed, check the input and output grid sizes" << std::endl;
	std::cout << "Or found that NGRID_IN and NGRID_OUT are equal, therefore no subsampling can happen" << std::endl;		
}

std::cout << "Based on the number of grid points for the input and the output, the subsampling (by averaging over neighbourhoods) will be performed: " << num_divisions << " times" << std::endl;


std::vector<float> field(NCOMP*FFTNORM_IN);


std::ifstream fin(binary_file,std::ios::binary);
fin.read(reinterpret_cast<char*>(field.data()), field.size()*sizeof(float));
fin.close();
std::cout << "The data loaded is at least: " << (field.size()*sizeof(float))/pow(1024,3) << "GB in size" << std::endl;

float x_max=field[std::distance(field.begin(), std::max_element(field.begin(), field.end()))];
float x_min=field[std::distance(field.begin(), std::min_element(field.begin(), field.end()))];
std::cout << "The maximum value of the original field is: " << x_max << ", and the minimum is: " << x_min << "." << std::endl; 


int num, i,j,k;
int loc_i, loc_j, loc_k;

int NGRID_TEMP=NGRID_IN;
int NGRID_PREV=NGRID_IN;
long int FFTNORM_TEMP, FFTNORM_PREV;

FFTNORM_TEMP=pow(NGRID_TEMP,3);
FFTNORM_PREV=pow(NGRID_PREV,3);

std::vector<float> fieldout(NCOMP*FFTNORM_TEMP/2,0);

for(num=0; num < num_divisions; ++num){

	NGRID_TEMP/=2;
	FFTNORM_TEMP=pow(NGRID_TEMP,3);
	FFTNORM_PREV=pow(NGRID_PREV,3);
        std::cout << "Iteration number " << num << ". NGRID_TEMP is equal to: " << NGRID_TEMP << " and FFTNORM_TEMP=" << FFTNORM_TEMP << std::endl;
	fieldout.resize(NCOMP*FFTNORM_TEMP, 0);

#ifdef PARALLEL
#pragma omp parallel for collapse(7) reduction(vec_plus_float:fieldout) num_threads(nthreads)	
#endif
for(i=0; i<NGRID_PREV; i+=2){
	for(j=0; j<NGRID_PREV; j+=2){
		for(k=0; k<NGRID_PREV; k+=2){
			for(loc_k=-1; loc_k<=1; loc_k++){
                          	for(loc_j=-1; loc_j<=1; loc_j++){
               		                for(loc_i=-1; loc_i<=1; loc_i++){
						for(int comp=0; comp< NCOMP; comp++){ 
 			                                fieldout[LOC_comp(LOC(i/2,j/2,k/2,NGRID_TEMP),comp)]+=field[LOC_comp(LOC((1+i+loc_i)%NGRID_PREV,(1+j+loc_j)%NGRID_PREV,(1+k+loc_k)%NGRID_PREV,NGRID_PREV) , comp)];
	
							}
						}
					
					}
			}
		}
	}
}

#ifdef PARALLEL
#pragma omp parallel for collapse(4) num_threads(nthreads)   
for(int comp=0; comp< NCOMP; comp++){
                for(i=0; i<NGRID_PREV; i+=2){
                        for(j=0; j<NGRID_PREV; j+=2){
                                for(k=0; k<NGRID_PREV; k+=2){
						fieldout[LOC_comp(LOC(i/2,j/2,k/2,NGRID_TEMP),comp)]/=27;
				}
			}
		}
}
#endif

#ifndef PARALLEL 
for(i=0; i<NGRID_PREV; i+=2){
        for(j=0; j<NGRID_PREV; j+=2){
                for(k=0; k<NGRID_PREV; k+=2){
                        for(loc_k=-1; loc_k<=1; loc_k++){
                                for(loc_j=-1; loc_j<=1; loc_j++){
                                        for(loc_i=-1; loc_i<=1; loc_i++){
                                                for(int comp=0; comp< NCOMP; comp++){ 
                                                        fieldout[LOC_comp(LOC(i/2,j/2,k/2,NGRID_TEMP),comp)]+=field[LOC_comp(LOC((1+i+loc_i)%NGRID_PREV,(1+j+loc_j)%NGRID_PREV,(1+k+loc_k)%NGRID_PREV,NGRID_PREV) , comp)];

                                                        }
                                                }

                                        }
                                                        fieldout[LOC_comp(LOC(i/2,j/2,k/2,NGRID_TEMP),comp)]/=27;
                                                         
                        }
                }
        }
}
#endif

        if(num==num_divisions-1){	
		std::ofstream fout(output,std::ios::out | std::ios::binary);
		fout.write((char*)&fieldout[0], fieldout.size()*sizeof(float));
		fout.close();
		x_max=fieldout[std::distance(fieldout.begin(), std::max_element(fieldout.begin(), fieldout.end()))];
		x_min=fieldout[std::distance(fieldout.begin(), std::min_element(fieldout.begin(), fieldout.end()))];
		std::cout << "The maximum value of the subsampled field is: " << x_max << ", and the minimum is: " << x_min << "." << std::endl;

		}

        else{
		NGRID_PREV=NGRID_TEMP;
		field=fieldout;
		}	
	

}

return 0;
}



long int LOC(int x_i, int y_j, int z_k, int NGRID){
return (x_i*NGRID*NGRID + y_j*NGRID + z_k);   //in index notation (i,j,k), where k is the fastest-changing index 
}

long int LOC_comp(long int LOC, int ncomp){
return ncomp+NCOMP*LOC;
}

