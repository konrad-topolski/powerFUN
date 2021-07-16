/* This simple program serves to load a binary file containing many fields and separate them into distinct binary files 
   As a matter of fact, this seemed simpler than adding yet another function to powerFUN program

   Compilation: g++ -std=c++17 separator.cpp -o separator -lstdc++fs

   Usage: ./separator binary_file output_pathname NGRID NCOMPS  <names of consecutive fields>
   Example:
   ./separator /home/topolski/DTFE512.a_velVort /home/topolski/DTFE512_component 512 3 v12 v13 v23

DATE: 06.04.2021
@:    k.topolski2@student.uw.edu.pl
														       */


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include <cmath>

int main(int argc, char* argv[]){
if(argc<6){
//throw std::invalid_argument("Not enough arguments supplied!");
std::cout << "Not enough arguments supplied!" << std::endl;
std::exit(1);
}

   std::string binary_file = argv[1];
   std::string output= argv[2];
   int NGRID    = atoi(argv[3]);
   int NCOMPS = atoi(argv[4]);

std::cout << "The binary file name and location is: " << binary_file << std::endl;
std::cout << "The output basename and location: " << output << " and the extensions will be given accordingly, as in provided arguments" << std::endl;
std::cout << "Grid size: " << NGRID << std::endl;
std::cout << "Number of fields to be extracted: " << NCOMPS << std::endl;

if(!std::filesystem::exists(binary_file)){
//throw std::invalid_argument("Invalid name of some input file! Check the name!");
std::cout << "Invalid name of some input file! Check the name!" << std::endl;
std::exit(1);
}


std::ifstream fin(binary_file,std::ios::binary); 
std::vector<float> fields(pow(NGRID,3)*NCOMPS);

std::cout << "Extracting individual components from the original data... " << std::endl;


float temp=0;	
for (int i=0; i<int(pow(NGRID,3)); i++){
	for(int j=0; j< NCOMPS; j++){
	fin.read((char*)&temp, sizeof(float));
	fields[j*pow(NGRID,3) + i]= temp;
	}
}
fin.close();



for(int i=5; i <argc; i++){
	std::cout << "Saving " << argv[i] << " to file: " << output+argv[i] << std::endl;
	std::ofstream fout(output+argv[i],std::ios::out | std::ios::binary);
	fout.write((char*)&fields[(i-5)*pow(NGRID,3)], pow(NGRID,3)*sizeof(float));
	fout.close();
}



return 0;
}
