/* This simple program glues together the data in a DTFE output binary file, which is obtained via the options --partNo
   Additionally, it must take the indices from the separate binary file 
   
   Compilation: 
   USE THE PrgEnv-gnu (with the GNU C++ compiler version 8.X)

   g++ -std=c++17 glue_DTFE.cpp -I /lustre/tetyda/home/topolski/include -L /lustre/tetyda/home/topolski/lib -lm -fopenmp -lstdc++fs -o glue_DTFE

   Usage: 

   ./glue_DTFE DTFE_loc output_loc  DTFE_binary_file_basename NGRID NCOMP NFILES 


   Example:
./glue_DTFE /home/topolski/dtfe_files/1024x1000Mpc/5/ /home/topolski/dtfe_files/1024x1000Mpc/5/  1024x1000Mpc_R025W7_LCDM_z0.000.DTFE1024.a_vel 1024 3 


MOST RECENT UPDATE:	22.05.2021
AUTHOR:			Konrad Topolski
@:			k.topolski2@student.uw.edu.pl
 
																       */

/*===================================================================================================================================================
Functions for picking out the element location
===================================================================================================================================================*/
long int LOC(int x_i, int y_j, int z_k);
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
#include <omp.h>
#include <numeric> 


#pragma omp declare reduction(vec_int_plus : \
		    std::vector<int> : \
                    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))      
                    

/*===================================================================================================================================================
GLOBAL VARIABLES
===================================================================================================================================================*/
std::string basename;
std::string path;
std::string output;
std::string extension;
int NGRID;
int NCOMP;
int NFILES;
int nthreads=omp_get_max_threads();
/*===================================================================================================================================================
THE MAIN FUNCTION 
===================================================================================================================================================*/

int main(int argc, char* argv[]){
 
if(argc<8){
	std::cout << "\n\t Not enough arguments supplied. Usage:\n \t ./glue_DTFE input_location output_loc basename field_extension NGRID NCOMP NFILES \n";
        std::cout << "\t The paths should be represented as follows: "  << std::endl;
	std::cout << "\t FILES_LOC/partNoX_BASENAME, and the index-keeping files should be saved as:" << std::endl;
	std::cout << "\t FILES_LOC/indices_partNoX_BASENAME" << std::endl;
	std::cout << std::endl;
	std::exit(1);
        }

//std::fstream indexFile("indices_partNo"+std::to_string(userOptions.partNo)+"_"+userOptions.outputFilename, 

   path  = argv[1];
   output = argv[2];
   basename = argv[3];
   extension = argv[4];
   NGRID=atoi(argv[5]);
   NCOMP=atoi(argv[6]);
   NFILES=atoi(argv[7]);


std::string filename;
std::string indexfilename;
int nx0,nxf,ny0,nyf,nz0,nzf; //indices to be read from the header of a specifically-prepared file 
std::cout << "\n\n\n";	
std::cout << "The binary file basename: " << basename << std::endl;
std::cout << "The binary file location: " << path << std::endl;
std::cout << "The binary file output and path: " << output+"/"+basename+extension << std::endl;
std::cout << "The file extension of the glued field is: " << extension << std::endl;
std::cout << "Grid size for the data: " << NGRID<< std::endl;
std::cout << "Number of components for the field: " << NCOMP << std::endl;
std::cout << "Number of files to glue data from: " << NFILES << std::endl;
std::cout << "The data saved will be at least: " << (pow(NGRID,3)*NCOMP*sizeof(float))/pow(1024,3) << "GB in size" << std::endl;




std::vector<float> field(NCOMP*pow(NGRID,3));
std::cout << "Allocated the vector" << std::endl;




int num;
float dummy;
std::cout <<"\n";
for(int num=0; num<NFILES; num++){
		       	        filename=path+"/partNo"+std::to_string(num)+"_"+basename+extension;
				std::cout << "[" << num << "/" << NFILES  << "]" << std::endl;
				std::cout << "Reading from file" << filename << std::endl;
				if(!std::filesystem::exists(filename)){
				     std::cout <<"Invalid name of the input file - the file doesn't exist! Check the name!" << std::endl;
				      }

				indexfilename=path+"/indices_partNo"+std::to_string(num)+"_"+basename;
                                //std::cout << "TEST:" << path+"/indices_partNo"+std::to_string(num)+"_"+basename <<  std::endl;
				//indexfilename="/home/topolski/DTFE/indices_partNo6_1024x1000Mpc_R06_W7_LCDM_z10.000.DTFE1024";
				std::cout << "Reading indices from file: " << indexfilename << std::endl;
				if(!std::filesystem::exists(indexfilename)){
                                     std::cout <<"Invalid name of the index file - the file doesn't exist! Check the name!" << std::endl;
                                      }
				

				std::ifstream IDfile(indexfilename,std::ios::binary | std::ios::in);
				IDfile.read((char*)&nx0, sizeof(int));
                                IDfile.read((char*)&nxf, sizeof(int));
                                IDfile.read((char*)&ny0, sizeof(int));
                                IDfile.read((char*)&nyf, sizeof(int));
                                IDfile.read((char*)&nz0, sizeof(int));
                                IDfile.read((char*)&nzf, sizeof(int));
				IDfile.close();

				std::cout  << "Reading indices: " <<  "[" << nx0 << " " << nxf << "] [" << ny0 << " " << nyf << "] [" << nz0 << " " << nzf << "]"<< std::endl;
				std::ifstream fin(filename, std::ios::binary | std::ios::in);
				for(int i=nx0; i<=nxf; i++){
					for(int j=ny0; j<=nyf ; j++){
						for(int k=nz0; k<=nzf ; k++){
		                                        for(int comp=0; comp<NCOMP; comp++){
										fin.read((char*)&dummy,sizeof(float));
										field[LOC_comp(LOC(i,j,k),comp)]=dummy;
								}
							}
						}
				}

				fin.close();
				std::cout << "\n";
}

std::cout << "Finished reading from all the files. Writing to a single binary file..." << std::endl;



std::ofstream fout(output+"/"+basename+extension,std::ios::out | std::ios::binary);
fout.write((char*)&field[0], field.size()*sizeof(float));
fout.close();

std::cout  << "Done. " << std::endl;

return 0;
}



long int LOC(int x_i, int y_j, int z_k){
return (x_i*NGRID*NGRID + y_j*NGRID + z_k);
}

long int LOC_comp(long int LOC, int ncomp){
return ncomp+NCOMP*LOC;
}

