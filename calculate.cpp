#include "functions.h"
#include "filters.h"

void calculate_all( char field1, char field2, std::string filename1, std::string filename2,std::string output_loc, std::string output_name, std::vector<std::string> FilterNames, std::vector<double> FilterRadii){

std::string FilterOneConvolveName   = FilterNames[0];
std::string FilterOneDeconvolveName = FilterNames[1];
std::string FilterTwoConvolveName   = FilterNames[2];
std::string FilterTwoDeconvolveName = FilterNames[3];

double RadiusConvolveOne   = FilterRadii[0];
double RadiusDeconvolveOne = FilterRadii[1]; 
double RadiusConvolveTwo   = FilterRadii[2]; 
double RadiusDeconvolveTwo = FilterRadii[3];

double scale_factor1=1/(1+redshift1);
double scale_factor2=1/(1+redshift2);
double hubble1=H1*sqrt(omega_matter*pow(scale_factor1,-3)+omega_lambda);
double hubble2=H2*sqrt(omega_matter*pow(scale_factor2,-3)+omega_lambda);


if(!std::filesystem::exists(filename1) || !std::filesystem::exists(filename2)){
throw std::invalid_argument("Invalid name of some input file! Check the names");
std::exit(1);
}

std::cout << "The value of scale_factor1 is: " << scale_factor1 << std::endl;
std::cout << "The value of scale_factor2 is: " << scale_factor2 << std::endl;
std::cout << "The value of the Hubble parameter corresponding to scale_factor1 is: " <<  hubble1 << std::endl;
std::cout << "The value of the Hubble parameter corresponding to scale_factor2 is: " <<  hubble2 << std::endl;  

std::cout<< "Calculating all the relevant quantities from two given fields" << std::endl;
std::cout << "OpenMP multi-threading is " << (fftw_init_threads()? "working." : "not working.") << std::endl;


fftw_complex * input1; 
input1 = (fftw_complex *) fftw_malloc(NGRID*NGRID*NGRID * sizeof(fftw_complex));

fftw_complex * output1; //the output array just after executing the plan will be the DFT of the density array 
output1= (fftw_complex*) fftw_malloc(NGRID*NGRID*NGRID * sizeof(fftw_complex));

fftw_complex * input2; 
input2 = (fftw_complex *) fftw_malloc(NGRID*NGRID*NGRID * sizeof(fftw_complex));

fftw_complex * output2; //the output array just after executing the plan will be the DFT of the density array 
output2= (fftw_complex*) fftw_malloc(NGRID*NGRID*NGRID * sizeof(fftw_complex));


fftw_complex * output3; //for storing data for the crosscorrelation function 
output3= (fftw_complex*) fftw_malloc(NGRID*NGRID*NGRID * sizeof(fftw_complex));

fftw_complex * output4; //for storing data for the pearson correlation coefficient (real part of crosscorrelation)
output4= (fftw_complex*) fftw_malloc(NGRID*NGRID*NGRID* sizeof(fftw_complex));

//PLANNING THE FFTW TRANSFORMS

//enable planning with nthreads
fftw_plan_with_nthreads(nthreads); //if the number of threads is not defined by the user, it will use the maximum number 
std::cout << "The program is using: " << nthreads << " threads for calculating FFTs. " <<  std::endl;

std::cout << "Creating a DTF 3D plan..." << std::endl;
fftw_plan plan;
plan=fftw_plan_dft_3d(NGRID, NGRID, NGRID, input1, output1, FFTW_FORWARD, FFTW_ESTIMATE);
/*You must create the plan before initializing the input, because FFTW_MEASURE overwrites the in/out arrays.  (Technically,FFTW_ESTIMATE does not touch your arrays, but you should always create plans first just to be sure.)*/

std::cout<<"Reading data from file..." << std::endl;

/*========================================
Unfortunately, because the reading of the file via the fstream happens sequentially (one entry line by one)
it is not possible to parallelise the reading from file this simply - this operation is inherently nonparallelisable
==========================================*/



if(filetype=="TextFile_xyz"){
std::ifstream file(filename1); //taking the datafile into the stream
	for(int i=0; i<NGRID;i++){
			for(int j=0; j<NGRID;j++){
				for(int k=0; k<NGRID; k++){
					file>>input1[LOC(i,j,k)][0]; 	//accessing the x index 
					file>>input1[LOC(i,j,k)][0];	//accessing the y index
					file>>input1[LOC(i,j,k)][0];	//accessing the z index
					file>>input1[LOC(i,j,k)][0];	//accessing the density value at that points
			                input1[LOC(i,j,k)][1]=0;
						
				}
			}
		}
		

file.close();
}
else if(filetype=="TextFile"){
std::ifstream file(filename1); //taking the datafile into the stream
     for(int i=0; i<NGRID;i++){
     	for(int j=0; j<NGRID;j++){
	   for(int k=0; k<NGRID; k++){
		file>>input1[LOC(i,j,k)][0];

	   }
       }
    }
file.close();
}
else if(filetype=="BinaryFile"){

std::ifstream file(filename1,std::ios::binary); 
std::vector<float> buf(NGRID*NGRID*NGRID);
file.read(reinterpret_cast<char*>(buf.data()), buf.size()*sizeof(float));
	for(int i=0; i<NGRID*NGRID*NGRID;i++){
		input1[i][0]=buf[i];
		input1[i][1]=0;
	}

}

else{
std::cout << "File format not detected!" << std::endl;
}



if(filetype=="TextFile_xyz"){
std::ifstream file(filename2); //taking the datafile into the stream
	for(int i=0; i<NGRID;i++){
			for(int j=0; j<NGRID;j++){
				for(int k=0; k<NGRID; k++){
					file>>input2[LOC(i,j,k)][0]; 	//accessing the x index 
					file>>input2[LOC(i,j,k)][0];	//accessing the y index
					file>>input2[LOC(i,j,k)][0];	//accessing the z index
					file>>input2[LOC(i,j,k)][0];	//accessing the density value at that points
					input2[LOC(i,j,k)][1]=0;	
				}
			}
		}
		

file.close();
}
else if(filetype=="TextFile"){
std::ifstream file(filename2); //taking the datafile into the stream
     for(int i=0; i<NGRID;i++){
     	for(int j=0; j<NGRID;j++){
	   for(int k=0; k<NGRID; k++){
		file>>input2[LOC(i,j,k)][0];

	   }
       }
    }
file.close();
}

else if(filetype=="BinaryFile"){
std::ifstream file(filename2,std::ios::binary); 
std::vector<float> buf(NGRID*NGRID*NGRID);
file.read(reinterpret_cast<char*>(buf.data()), buf.size()*sizeof(float));
        for(int i=0; i<NGRID*NGRID*NGRID;i++){
                input2[i][0]=buf[i];
                input2[i][1]=0;
        }

}
else{
std::cout << "File format not detected!" << std::endl;
}






if(field1=='T'){
        std::cout <<"The chosen type of field1 is: " << field1 << " - theta" << std::endl;
        for( int i=0; i<NGRID*NGRID*NGRID; i++ ){
                input1[i][0]/=-scale_factor1*hubble1;
      }
}
else if (field1=='D'){
	std::cout <<"The chosen type of field1 is: " << field1 << " - density" << std::endl;
	std::cout << "Calculating the density average at redshift1..." << std::endl;
	double dens_avg1=density_average(input1); //computing the average of density values
	std::cout << "The density average at redshift1 is: " << dens_avg1 << std::endl;
	input1=dens_contrast(input1,dens_avg1);	  //now the input array is the density contrast array δ(x)
}
else{
	std::cout <<"The chosen type of field1 is: " << field1 << " - other" << std::endl;
	std::cout << "No preliminary scaling will be performed" << std::endl;
}

if(field2=='T'){
        std::cout <<"The chosen type of field2 is: " << field2 << " - theta" << std::endl;
        for( int i=0; i<NGRID*NGRID*NGRID; i++ ){
                input2[i][0]/=-scale_factor2*hubble2;
	}
}
else if (field2=='D'){
        std::cout <<"The chosen type of field2 is: " << field2 << " - density" << std::endl;
        std::cout << "Calculating the density average at redshift2..." << std::endl;
        double dens_avg2=density_average(input2); //computing the average of density values
        std::cout << "The density average at redshift2 is: " << dens_avg2 << std::endl;
        input2=dens_contrast(input2,dens_avg2);   //now the input array is the density contrast array δ(x)
}
else{
        std::cout <<"The chosen type of field2 is: " << field2 << " - other" << std::endl;
        std::cout << "No preliminary scaling will be performed" << std::endl;
}


std::cout << "Reading from files: " << filename1 << ", " << filename2 << " completed." << std::endl;


std::cout << "Executing the plan for the first input..." << std::endl;
fftw_execute(plan); //execute it with nthreads - now the output array is the density contrast in the fourier space - δ(k)
std::cout << "Executing the plan for the second input (reusing the plan)..." << std::endl;
fftw_execute_dft(plan,input2, output2);

fftw_destroy_plan(plan); //destroy the plan immediately afterwards, before cleaning up threads and planning new things
fftw_cleanup_threads();
fftw_plan_with_nthreads(1); 
std::cout<<"The plan executed successfully!" << std::endl;
 //come back to single-threaded operations (must be thread-safe!)



fftw_free(input1);
fftw_free(input2);


//AT THIS POINT WE HAVE FIELD1 AND FIELD2 IN FOURIER SPACE 
// Here we apply convolution and deconvolution 
std::cout << "Applying filters..." << std::endl;
FilterInterface(output1, FilterOneConvolveName, FilterOneDeconvolveName, RadiusConvolveOne, RadiusDeconvolveOne);
FilterInterface(output2, FilterTwoConvolveName, FilterTwoDeconvolveName, RadiusConvolveTwo, RadiusDeconvolveTwo);
// filters applied...

int ks;
double kx,ky,kz;
double fcorrect;

std::cout << "Getting f1f1*, f2f2*, f1f2*" << std::endl;


#pragma omp parallel for num_threads(nthreads)
for(int i=0; i<NGRID*NGRID*NGRID;i++){

output3[i][0] = pow(pow(output1[i][0]*output2[i][0],2) + pow(output1[i][0]*output2[i][1],2) + pow(output1[i][1]*output2[i][0],2) + pow(output1[i][1]*output2[i][1],2),0.5);
output3[i][1] = 0;
output4[i][0] = output1[i][0]*output2[i][0] + output1[i][1]*output2[i][1];//for an alternative way of calculating C(k)
output4[i][1] = 0;


output1[i][0]=output1[i][0]*output1[i][0] + output1[i][1]*output1[i][1];
output1[i][1]=0;

output2[i][0]=output2[i][0]*output2[i][0] + output2[i][1]*output2[i][1];
output2[i][1]=0;



}


std::vector<int> kmod(NGRID);
   for(int i=0; i<NGRID; i++)
     {
      if(i <= NYQUIST)
         kmod[i] = i;
      else
         kmod[i] = -(NGRID-i);
     }


std::vector<int> jpower(NYQUIST,0);		//counting the instances of a given wavenumber power contribution, initialize all entries to zero

std::vector<double> Pk1(NYQUIST,0);			//the normalised power spectrum
std::vector<double> Pk1_err(NYQUIST,0);			//the normalised power spectrum error
std::vector<double> dpower1(NYQUIST,0);		//unnormalised power spectrum value at a given frequency/wavenumber


std::vector<double> Pk2(NYQUIST,0);			//the normalised power spectrum
std::vector<double> Pk2_err(NYQUIST,0);			//the normalised power spectrum error
std::vector<double> dpower2(NYQUIST,0);		//unnormalised power spectrum value at a given frequency/wavenumber


std::vector<double> Crossk(NYQUIST,0);			//the normalised part of the cross-correlation function
std::vector<double> Crossk_err(NYQUIST,0);			//the error
std::vector<double> dcross(NYQUIST,0);		//unnormalised cross-correlation function at a given frequency/wavenumber


std::vector<double> Crossk_real(NYQUIST,0);			//the normalised real part of the cross-correlation function
std::vector<double> Crossk_real_err(NYQUIST,0);			//the error
std::vector<double> dcross_real(NYQUIST,0);		//unnormalised cross-correlation - real part - function at a given frequency/wavenumber


std::vector<double> Cofk(NYQUIST,0);			//the C(k) correlation coefficient
std::vector<double> Cofk_err(NYQUIST,0);			//the C(k) correlation coefficient errors


			
#pragma omp parallel for collapse(3) private(kx,ky,kz,ks,fcorrect) num_threads(nthreads) reduction(vec_int_plus:jpower) reduction(vec_double_plus:dcross,dcross_real,dpower1,dpower2)
for(int i=0; i<NGRID; i++){
      for(int j=0; j<NGRID; j++){
         for(int k=0; k<NGRID;k++)  {
           
            ks = (int) floor( sqrt( (double)pow(kmod[i],2) +
                                    (double)pow(kmod[j],2) +
                                    (double)pow(kmod[k],2) )); 
            fcorrect = 1.0;
            
 	    #ifdef TSC_CORRECT
            kx = M_PI*i/(double)NGRID;
            ky = M_PI*j/(double)NGRID;
            kz = M_PI*k/(double)NGRID;            
            if(kx > 0.) fcorrect *= win_func(kx);//sin(kx)/kx;
            if(ky > 0.) fcorrect *= win_func(ky);//sin(ky)/ky;
            if(kz > 0.) fcorrect *= win_func(kz);//sin(kz)/kz;            
            fcorrect = 1./pow(fcorrect,3);
	    #endif

            if(ks >= 1 && ks <= NYQUIST)
              {
               jpower[ks]  += 1;
               dpower1[ks] += fcorrect*output1[LOC(i,j,k)][0];
		dpower2[ks] +=fcorrect*output2[LOC(i,j,k)][0];
               dcross[ks]  += fcorrect*output3[LOC(i,j,k)][0];
               dcross_real[ks] += fcorrect*output4[LOC(i,j,k)][0];
           
              }

           }}}





double FFTnorm=NGRID*NGRID*NGRID;
double volume=pow(L,3);


#pragma omp parallel for num_threads(nthreads) 
for(int i =0; i<NYQUIST;i++){

if(jpower[i]!=0){

Pk1[i] =  (dpower1[i]/jpower[i])*(volume/pow(FFTnorm,2));
Pk1_err[i] = Pk1[i]/pow(jpower[i],0.5);


Pk2[i] =  (dpower2[i]/jpower[i]) *(volume/pow(FFTnorm,2));
Pk2_err[i] = Pk2[i]/pow(jpower[i],0.5);


Crossk[i] =  (dcross[i]/jpower[i]) * (volume/pow(FFTnorm,2));
Crossk_err[i]=Crossk[i]/pow(jpower[i],0.5);

Crossk_real[i] =  (dcross_real[i]/jpower[i]) * (volume/pow(FFTnorm,2));
Crossk_real_err[i]=Crossk_real[i]/pow(jpower[i],0.5);
}

else{
Pk1[i]=0;
Pk1_err[i]=0;
Pk2[i]=0;
Pk2_err[i]=0;
Crossk[i]=0;
Crossk_err[i]=0;
Crossk_real[i]=0;
Crossk_real_err[i]=0;}

if(Pk1[i]!=0 and Pk2[i]!=0){
Cofk[i] =Crossk_real[i]/sqrt(Pk1[i]*Pk2[i]);
Cofk_err[i]=0;
Cofk_err[i]+=pow(Pk2_err[i]/Pk2[i],2)+pow(Pk1_err[i]/Pk1[i],2);
Cofk_err[i]*=Crossk_real[i]*Crossk_real[i]/4;
Cofk_err[i]+=pow(Crossk_real_err[i],2);
Cofk_err[i]/=Pk1[i]*Pk2[i];
Cofk_err[i]=sqrt(Cofk_err[i]);
}
else{Cofk[i]=0;
Cofk_err[i]=0;}

}


std::cout << "Writing the data  to textfile..." << std::endl;
std::ofstream outfile;
outfile.open(output_loc+"powerspectra-"+output_name,std::ios::trunc | std::ios::out);

if(outfile.is_open()){
outfile << "#" << " k  P_f1(k)  deltaP_f1(k)  P_f2(k)  deltaP_f2(k)  P_f1f2(k)  deltaP_f1f2(k)  pearson(k)  deltapearson(k)  R(k) deltaR(k) \n";
for(int i=1;i<NYQUIST;i++){
outfile <<  (2*i+1)*M_PI/L 
<< " " 		//the errors are output in absolute value
<< Pk1[i] 
<< " " 
<< abs(Pk1_err[i])
<< " "  
<< Pk2[i]
<< " " 
<< abs(Pk2_err[i])
<< " " 
<< Crossk[i] 
<< " " 
<< abs(Crossk_err[i])
<< " "
<< Cofk[i] //Pearson correlation coefficient - its sign is important! 
<< " " 
<< abs(Cofk_err[i])
<< " " 
<< (field1=='D'? (field2=='T'? Pk2[i]/(Pk1[i]) : 0) : (field2=='D' ? Pk1[i]/(Pk2[i]): 0 )) 
<< " "
<< (field1=='D'? (field2=='T'? pow(pow(Pk2_err[i],2)+pow(Pk2[i]*Pk1_err[i]/Pk1[i],2),0.5)/(Pk1[i]) : 0) : (field2=='D' ? pow(pow(Pk1_err[i],2)+pow(Pk1[i]*Pk2_err[i]/Pk2[i],2),0.5)/(Pk2[i]) : 0 ))
<< "\n";
}}

outfile.close();


fftw_free(output1); 
fftw_free(output2); 
fftw_free(output3); 
fftw_free(output4); 

std::cout << "The calculations of (k,P_f1(k),deltaP_f1(k),P_f2(k),deltaP_f2(k), P_f1f2(k), deltaP_f1f2(k), pearson(k), deltapearson(k), R(k), deltaR(k)) are finished." << std::endl;

}





