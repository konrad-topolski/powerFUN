This script is designed to work on the slight modification of the DTFE program, where the following change has been made:
in /src/subpartition.h, in the place where indices of the subpartition are given
nx [    ]
ny [    ]
nz [    ]

copy and paste this:

//MY_OWN_CODE   
std::cout << "Indices will be saved to: " << "/home/topolski/DTFE/partNo"+std::to_string(userOptions.partNo)+"_"+userOptions.outputFilename.substr(((std::string)"/home/topolski/DTFE/partNoX").length()) << std::endl;
std::ofstream indexFile("/home/topolski/DTFE/partNo"+std::to_string(userOptions.partNo)+"_"+userOptions.outputFilename.substr(((std::string)"/home/topolski/DTFE/partNoX").length()), std::ios::binary | std::ios::out); 
for(int i=0; i<6; i++){
	int temp =subgrid[userOptions.partNo][i];
	if(i%2==1){temp-=1;}
	indexFile.write((char*)&temp, sizeof(int)); 
	}
indexFile.close();
//END_OF_MY_OWN_CODE

and remember to choose the correct directory to save the indices

also, at the top of subpartition.h, add 
#include <fstream>
