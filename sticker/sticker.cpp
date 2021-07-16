/* This simple program serves to load a binary file containing many field components and join them into a binary file containing:
   (*) an average of field components, e.x. (v1+v2+v3)/3
   (*) the norm of the fields sqrt(v1*v1 + v2*v2 + v3*v3)
   As a matter of fact, this seemed simpler than adding yet another function to powerFUN program

   Compilation: 
   g++ -std=c++17 sticker.cpp -o sticker -lstdc++fs

   Usage: ./sticker.exe binary_file output_pathname NGRID NCOMPS mode
   where mode is a string variable equal to either: "average" or "norm"

   
   Example:
   ./sticker /home/topolski/DTFE512.a_velVort /home/topolski/DTFE512_component 512 3 average

DATE: 26.05.2021
@:    k.topolski2@student.uw.edu.pl
														       */


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include <math.h>

int main(int argc, char* argv[]){

std::cout << std::endl;

if(argc<6){
//throw std::invalid_argument("Not enough arguments supplied!");
std::cout << "Not enough arguments supplied!" << std::endl;
std::exit(1);
}


if(!std::filesystem::exists(argv[1])){
//throw std::invalid_argument("Invalid name of some input file! Check the name!");
std::cout << "Invalid name of some input file! Check the name!" << std::endl;
std::exit(1);
}

std::string mode = argv[5];
std::string average = "average";
std::string norm = "norm";


if(mode.compare(norm)!=0  & mode.compare(average)!=0 ){
//throw std::invalid_argument("Invalid mode of operation! It should be either average or norm!");
std::cout << "Invalid mode of operation! It should be either average or norm!" << std::endl;
std::exit(1);
}



   std::string binary_file = argv[1];
   std::string output= argv[2];
   int NGRID    = atoi(argv[3]);
   int NCOMPS = atoi(argv[4]);

std::cout << "The binary file name and location is: " << binary_file << std::endl;
std::cout << "The output name and location: " << output+"_"+mode << std::endl;
std::cout << "Grid size: " << NGRID << std::endl;
std::cout << "Number of field components to join: " << NCOMPS << std::endl;
std::cout << "The mode of operation is to calculate the: " << mode << std::endl;



std::ifstream fin(binary_file,std::ios::binary); 
std::vector<float> field(pow(NGRID,3),0);

std::cout << "Extracting individual components from the original data... " << std::endl;


float temp=0;	
for ( int i=0; i < int(pow(NGRID,3)) ; i++ ){
	for ( int j=0 ; j < NCOMPS ; j++ ){

		fin.read((char*)&temp, sizeof(float));

		if(mode=="average"){
			field[i]+= temp;
			}

		else if(mode=="norm"){
			field[i]+= temp*temp;
			}
		}
	if(mode=="average"){
		field[i]/=NCOMPS;
		}

	else if(mode=="norm"){
		field[i]=sqrt(field[i]);
		}
}

fin.close();



std::cout << "Saving the " << mode  << " of binary file " << binary_file << " to file: " << output+"_"+mode << std::endl;
std::ofstream fout(output+"_"+mode,std::ios::out | std::ios::binary);
fout.write((char*)&field[0], pow(NGRID,3)*sizeof(float));
fout.close();
std::cout << std::endl;


return 0;
}
