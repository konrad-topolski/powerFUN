#include <iostream>
#include "constants.h"
#include "spectra.h"

using namespace boost;
namespace po = boost::program_options;


int main(int argc, char* argv[]){

/*=================================================================================================================================================
SPACE FOR OPTIONS
================================================================================================================================================*/

bool flag_DTFE = false;
po::options_description description("Usage:");

description.add_options()
("help,h", "Display this help message. The program will calculate the powerspectra of field1 and field2, it recognized three types of fields and scales them appropriately.")
("grid,gp",po::value<int>()->default_value(512),"Number of grid points used by both DTFE and powerFUN")
("filetype,ft",po::value<std::string>()->default_value("BinaryFile"),"Filetype (format) which is the output of DTFE and input of powerFUN. Available options:\nTextFile\nTextFile_xyz\nBinaryFile")
("box,length,l,l0",po::value<double>()->default_value(500),"The length of the box of the simulation (in Mpc/h) (check with the output of DTFE)")
("hubble,H0,h0",po::value<double>()->default_value(70.4),"The value of the Hubble constant at present time (in (km/s)/Mpc ) ")
("hubble1,H1,h1",po::value<double>()->default_value(70.4),"The value of the Hubble constant at present time (H_0) corresponding to field1 data suite (km/s)/Mpc ) ")
("hubble2,H2,h2",po::value<double>()->default_value(70.4),"The value of the Hubble constant at present time (H_0) corresponding to field2 data suite (in (km/s)/Mpc ) ")
("redshift,Z0,z0",po::value<double>()->default_value(0),"The value of redshift for considered data")
("redshift1,Z1,z1",po::value<double>()->default_value(0),"The value of redshift1 for considered data")
("redshift2,Z2,z2",po::value<double>()->default_value(0),"The value of redshift2 for considered data")
("FilterOneConvolve", po::value<std::vector<std::string>> () ->multitoken(), "Set the convolution filter for the field1, in the order of words: [FilterName Radius]. Filters: NGP CIC TSC TopHat Gaussian None.")
("FilterOneDeconvolve", po::value<std::vector<std::string>> () ->multitoken(), "Set the deconvolution filter for the field1, in the order of words: [FilterName Radius]. Filters: NGP CIC TSC TopHat Gaussian None.")
("FilterTwoConvolve", po::value<std::vector<std::string>> () ->multitoken(), "Set the convolution filter for the field2, in the order of words: [FilterName Radius]. Filters: NGP CIC TSC TopHat Gaussian None.")
("FilterTwoDeconvolve", po::value<std::vector<std::string>> () ->multitoken(), "Set the deconvolution filter for the field2, in the order of words: [FilterName Radius]. Filters: NGP CIC TSC TopHat Gaussian None.")
("f1",po::value<char>()->default_value('D'),"The type of field1 (D=density, T=theta, N=other) for considered data")
("f2",po::value<char>()->default_value('T'),"The type of field2 (D=density, T=theta, N=other) for considered data")
("omega_matter,OM",po::value<double>()->default_value(0.272),"Omega matter value for considered data")
("omega_lambda,OL",po::value<double>()->default_value(0.728),"Omega lambda value for considered data")
("calculate_all", po::value<std::vector<std::string>> () ->multitoken(), "A function taking two fields and returning their autocorrelations, crosscorrelation, Pearson correlation coefficient and ratio of recalled spectra")
("output_loc", po::value<std::string>()->default_value("/home/topolski"), "The location of the output for powerFUN")
("output_name", po::value<std::string>()->default_value("MyOutput"), "The name of the data files from DTFE (without extensions, i.e. without e.x. .a_den etc.)");




po::variables_map vm; //variables map contains the options values (2nd column above, in "description" object)
po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
po::notify(vm);

if (vm.count("help")){
std::cout << description;
}

std::cout << std::endl;
if (vm.count("grid")){
NGRID= vm["grid"].as<int>();
std::cout << "Number of grid points is set to: " << NGRID << std::endl;
}

NYQUIST=NGRID/2;
std::cout << "The Nyquist number is equal to: " << NYQUIST << std::endl;

if(vm.count("box")){
L=vm["box"].as<double>();
std::cout << "The box length is set to: " << L << " Mpc/h" << std::endl;
}

if(vm.count("hubble")){
H0=vm["hubble"].as<double>();
std::cout << "The Hubble constant is set to: " << H0 << " km/s / Mpc" << std::endl;
}

if(vm.count("hubble1")){
H1=vm["hubble1"].as<double>();
std::cout << "The hubble1 is set to: " << H1 << " km/s / Mpc" << std::endl;
}

if(vm.count("hubble2")){
H2=vm["hubble2"].as<double>();
std::cout << "The hubble2 is set to: " << H2 << " km/s / Mpc" << std::endl;
}



if(vm.count("redshift")){
redshift=vm["redshift"].as<double>();
std::cout << "The redshift value is set to: " << redshift << std::endl;
}
if(vm.count("redshift1")){
redshift1=vm["redshift1"].as<double>();
std::cout << "The redshift1 value is set to: " << redshift1 << std::endl;
}

if(vm.count("redshift2")){
redshift2=vm["redshift2"].as<double>();
std::cout << "The redshift2 value is set to: " << redshift2 << std::endl;
}


if(vm.count("omega_matter")){
omega_matter=vm["omega_matter"].as<double>();
std::cout << "The omega_matter value is set to: " << omega_matter << std::endl;
}

if(vm.count("omega_lambda")){
omega_lambda=vm["omega_lambda"].as<double>();
std::cout << "The omega_lambda value is set to: " << omega_lambda << std::endl;
}

std::string  FilterOneConvolveName = "None";
std::string  FilterOneDeconvolveName = "None";
std::string  FilterTwoConvolveName = "None";
std::string  FilterTwoDeconvolveName = "None";

double RadiusConvolveOne   = 0;
double RadiusDeconvolveOne = 0; 
double RadiusConvolveTwo   = 0; 
double RadiusDeconvolveTwo = 0;

std::vector<std::string> FilterNames(4);
std::vector<double> FilterRadii(4);


if(vm.count("FilterOneConvolve")){
std::vector<std::string> Data;
Data = vm["FilterOneConvolve"].as<std::vector<std::string>>();
FilterOneConvolveName=Data[0];
RadiusConvolveOne =std::stod(Data[1]);
std::cout << "The FilterOne for convolution is chosen to be: " << FilterOneConvolveName << " of radius " << RadiusConvolveOne << " Mpc/h" << std::endl;
}



if(vm.count("FilterOneDeconvolve")){
std::vector<std::string> Data;
Data = vm["FilterOneDeconvolve"].as<std::vector<std::string>>();
FilterOneDeconvolveName=Data[0];
RadiusDeconvolveOne=std::stod(Data[1]);
std::cout << "The FilterOne for deconvolution is chosen to be: " << FilterOneDeconvolveName << " of radius " << RadiusDeconvolveOne << " Mpc/h" << std::endl;
}



if(vm.count("FilterTwoConvolve")){
std::vector<std::string> Data;
Data = vm["FilterTwoConvolve"].as<std::vector<std::string>>();
FilterTwoConvolveName=Data[0];
RadiusConvolveTwo=std::stod(Data[1]);
std::cout << "The FilterTwo for convolution is chosen to be: " << FilterTwoConvolveName << " of radius " << RadiusConvolveTwo << " Mpc/h" << std::endl;
}


if(vm.count("FilterTwoDeconvolve")){
std::vector<std::string> Data;
Data = vm["FilterTwoDeconvolve"].as<std::vector<std::string>>();
FilterTwoDeconvolveName=Data[0];
RadiusDeconvolveTwo=std::stod(Data[1]);
std::cout << "The FilterTwo for deconvolution is chosen to be: " << FilterTwoDeconvolveName << " of radius " << RadiusDeconvolveTwo << " Mpc/h" << std::endl;
}

FilterNames[0] = FilterOneConvolveName;
FilterNames[1] = FilterOneDeconvolveName;
FilterNames[2] = FilterTwoConvolveName;
FilterNames[3] = FilterTwoDeconvolveName;


FilterRadii[0]= RadiusConvolveOne;
FilterRadii[1]= RadiusDeconvolveOne;
FilterRadii[2]= RadiusConvolveTwo;
FilterRadii[3]= RadiusDeconvolveTwo;



if(vm.count("filetype")){
filetype=vm["filetype"].as<std::string>();
std::cout << "The file format used is set to: " << filetype << std::endl;
}



std::string output_name=vm["output_name"].as<std::string>();
std::string output_loc=vm["output_loc"].as<std::string>();

std::cout << "The output location is: " << output_loc << std::endl;
std::cout << "The output filename is: " << output_name << std::endl;



/*================================================================================================================================================
END OF SPACE FOR PRELIMINARY OPTIONS
=========================================================================================================================================================*/

if (vm.count("calculate_all")){
std::vector<std::string> fieldnames = vm["calculate_all"].as<std::vector<std::string>>();
calculate_all(vm["f1"].as<char>(),vm["f2"].as<char>(),fieldnames[0],fieldnames[1], output_loc, output_name, FilterNames, FilterRadii);
std::cout << "All calculations finished." << std::endl;
}


return 0;
}




