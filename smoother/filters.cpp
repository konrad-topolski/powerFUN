#include "filters.h"
#include <iostream> 
#include <math.h> 

extern int NGRID;
extern double BOXSIZE;

#define LOC(ix,iy,iz)  (((ix)*NGRID*NGRID + (iy)*NGRID + (iz))    )


void FilterInterface(fftw_complex *field, std::string FilterConvolveName, std::string FilterDeconvolveName, std::string TypeConvolve, std::string TypeDeconvolve, double RadiusConvolve, double RadiusDeconvolve){
  
	std::cout << "The following options are used: " << std::endl;
	std::cout << "NGRID:" << NGRID << " BOXSIZE: " << BOXSIZE << " Mpc/h" << std::endl;
	std::cout << "Convolution: " << FilterConvolveName << " " << TypeConvolve << " " << RadiusConvolve << std::endl;
        std::cout << "Deconvolution: " << FilterDeconvolveName << " " << TypeDeconvolve << " " << RadiusDeconvolve << std::endl;
    
	if(FilterDeconvolveName.compare("TopHat")==0){
			std::cout << "Applying TopHat deconvolution of radius " << RadiusDeconvolve << " Mpc/h" << std::endl;
			TopHat(field,"Deconvolve", RadiusDeconvolve);	
		}
	else if(FilterDeconvolveName.compare("Gaussian")==0){
	                std::cout << "Applying Gaussian deconvolution of radius " << RadiusDeconvolve << " Mpc/h" << std::endl;
		        Gaussian(field,"Deconvolve", RadiusDeconvolve);   
		}
	else if(FilterDeconvolveName.compare("NGP")==0){
	                std::cout << "Applying NearestGridPoint deconvolution, type " << TypeConvolve  << ", of radius " << RadiusDeconvolve << " Mpc/h" << std::endl;
			if(TypeDeconvolve.compare("Radial")==0){ 
				 FilterFamilyRadial(field,1,"Deconvolve", RadiusDeconvolve);   
							}
	                if(TypeDeconvolve.compare("Cube")==0){ 
        	                 FilterFamily(field,1,"Deconvolve", RadiusDeconvolve);   
                                                }
		}
        else if(FilterDeconvolveName.compare("CIC")==0){
                	std::cout << "Applying CloudInCell deconvolution, type " << TypeConvolve  << ", of radius " << RadiusDeconvolve << " Mpc/h" << std::endl;
	                if(TypeDeconvolve.compare("Radial")==0){ 
	                         FilterFamilyRadial(field,2,"Deconvolve", RadiusDeconvolve);   
        	                                        }
	                if(TypeDeconvolve.compare("Cube")==0){ 
	                         FilterFamily(field,2,"Deconvolve", RadiusDeconvolve);   
        	                                        }
        	}
        else if(FilterDeconvolveName.compare("TSC")==0){
                	std::cout << "Applying TriangularShapedCloud deconvolution, type " << TypeConvolve  << ", of radius " << RadiusDeconvolve << " Mpc/h" << std::endl;
	                if(TypeDeconvolve.compare("Radial")==0){ 
	                         FilterFamilyRadial(field,3,"Deconvolve", RadiusDeconvolve);   
        	                                        }
	                if(TypeDeconvolve.compare("Cube")==0){ 
	                         FilterFamily(field,3,"Deconvolve", RadiusDeconvolve);   
	                                                }
        	}




	if(FilterConvolveName.compare("TopHat")==0){
                std::cout << "Applying TopHat convolution of radius " << RadiusConvolve << " Mpc/h" << std::endl;
	        TopHat(field,"Convolve", RadiusConvolve);   
        }
	else if(FilterConvolveName.compare("Gaussian")==0){
                std::cout << "Applying Gaussian convolution of radius " << RadiusConvolve << " Mpc/h" << std::endl;
	        Gaussian(field,"Convolve", RadiusConvolve);   
        }
        else if(FilterConvolveName.compare("NGP")==0){
                        std::cout << "Applying NearestGridPoint convolution, type " << TypeConvolve  << ", of radius " << RadiusConvolve << " Mpc/h" << std::endl;
                        if(TypeConvolve.compare("Radial")==0){ 
                                 FilterFamilyRadial(field,1,"Convolve", RadiusConvolve);   
                                                        }
                        if(TypeConvolve.compare("Cube")==0){ 
                                 FilterFamily(field,1,"Convolve", RadiusConvolve);   
                                                }
                }
        else if(FilterConvolveName.compare("CIC")==0){
                        std::cout << "Applying CloudInCell convolution, type " << TypeConvolve  << ", of radius " << RadiusConvolve << " Mpc/h" << std::endl;
                        if(TypeConvolve.compare("Radial")==0){ 
                                 FilterFamilyRadial(field,2,"Convolve", RadiusConvolve);   
                                                        }
                        if(TypeConvolve.compare("Cube")==0){ 
                                 FilterFamily(field,2,"Convolve", RadiusConvolve);   
                                                        }
                }
        else if(FilterConvolveName.compare("TSC")==0){
                        std::cout << "Applying TriangularShapedCloud convolution, type " << TypeConvolve  << ", of radius " << RadiusConvolve << " Mpc/h" << std::endl;
                        if(TypeConvolve.compare("Radial")==0){ 
                                 FilterFamilyRadial(field,3,"Convolve", RadiusConvolve);   
                                                        }
                        if(TypeConvolve.compare("Cube")==0){ 
                                 FilterFamily(field,3,"Convolve", RadiusConvolve);   
                                                        }
                }

}

void FilterFamilyRadial(fftw_complex* field, int P, std::string mode, double Radius){

   double fsize= Radius*double(NGRID)/double(BOXSIZE);

   std::function<double(const double&)> RadialFilter;

   	if(P==1){
			 RadialFilter = [](const double& x){ 
						return 3.0 * ( sin(x) - x * cos(x) ) /pow(x,3.0);
						};
			}
   else if(P==2){
			 RadialFilter = [](const double& x){ 
                                                return 12.0 * ( 2 - 2 * cos(x) - x * sin(x) )/pow(x,4.0);
                                                };
			}
   else if(P==3){
			 RadialFilter = [](const double& x){ 
                                                return (243.0/4.0) * ( 9 * sin(x/3) - 3 * sin(x) + x * cos(x) - x * cos(x/3) )/pow(x,5.0);
                                                };
			}

   else std::cout << "Wrong filter chosen!" << std::endl;

   std::function<double(const double& )> function;
   if(mode.compare("Convolve")==0){
                function = [](const double&  x){return x;};
        }
   else if(mode.compare("Deconvolve")==0){
                function = [](const double&  x){if(x!=0) return 1/x ;
                                                 else     return 1e6 ;} ;
        }

   int ni,nf,icx,icy,icz,isx,isy,isz;
   float irz,iry,irx, zk, rks, tmp, cfilter=0.0;

    field[LOC(0,0,0)][0] *= 1.0; field[LOC(0,0,0)][1] *= 1.0;
    rks = fsize*2.0*M_PI/(float)NGRID;
    ni = 2;
    nf = NGRID/2 +1;

// os x
   icz =1;
   icy =1;

#pragma omp parallel for private(isx,zk,tmp) shared(function, RadialFilter) num_threads(nthreads) 
   for(icx=ni; icx<=nf; icx++)
    {
     isx = NGRID - icx +2;
     zk = ((float)icx-1) * rks; // here we get k*R 
     tmp = RadialFilter(zk);
     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
    }

// os y
   icz =1;
   icx =1;

#pragma omp parallel for private(isy,zk,tmp) shared(function, RadialFilter) num_threads(nthreads)                                                                                             
   for(icy=ni; icy<=nf; icy++)
    {
     isy = NGRID - icy +2;
     zk = ((float)icy-1) * rks;
     tmp = RadialFilter(zk);

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
    }

// os z
   icx =1;
   icy =1;

#pragma omp parallel for  private(isz,zk,tmp) shared(function, RadialFilter) num_threads(nthreads)                                                                                             
   for(icz=ni; icz<=nf; icz++)
    {
     isz = NGRID - icz +2;
     zk = ((float)icz-1) * rks;
     tmp = RadialFilter(zk);

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
    }

// k_z =0 plane
  icz = 1;

#pragma omp parallel for collapse(2) private(isy,iry,isx,irx,zk,tmp) shared(function, RadialFilter) num_threads(nthreads)                                                                                             
  for(icy = ni; icy <=nf; icy++)
   {
    for(icx=ni; icx<=nf; icx++)
     {
      isy = NGRID- icy +2;
      iry = pow((float)icy-1.0,2.0);
      isx = NGRID - icx +2;
      irx = pow((float)icx-1.0,2.0);
      zk = sqrt(irx+iry) * rks;
      tmp = RadialFilter(zk);

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
     }
   }

// k_y=0 plane
  icy = 1;

#pragma omp parallel for collapse(2) private(isx, irx, isz, irz, zk, tmp)  shared(function, RadialFilter) num_threads(nthreads)                                                                                             
  for(icx = ni; icx <=nf; icx++)
   {
    for(icz=ni; icz<=nf; icz++)
     {

      isx = NGRID- icx +2;
      irx = pow((float)icx-1.0,2.0);
      isz = NGRID - icz +2;
      irz = pow((float)icz-1.0,2.0);
      zk = sqrt(irx+irz) * rks;
      tmp = RadialFilter(zk);

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
     }
   }

// k_x=0 plane
  icx = 1;

#pragma omp parallel for collapse(2) private(isy,iry,isz,irz,zk,tmp) shared(function, RadialFilter) num_threads(nthreads)                                                                                             
  for(icy = ni; icy <=nf; icy++)
   {
    for(icz=ni; icz<=nf; icz++)
     {


      isy = NGRID- icy +2;
      iry = pow((float)icy-1.0,2.0);
      isz = NGRID - icz +2;
      irz = pow((float)icz-1.0,2.0);
      zk = sqrt(iry+irz) * rks;
      tmp = RadialFilter(zk);

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
     }
   }

// the rest is volume-symetric

#pragma omp parallel for collapse(3) private(isz,irz,isy,iry,isx,irx,zk,tmp) shared(function, RadialFilter) num_threads(nthreads)                                                                                             
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
        tmp = RadialFilter(zk);
  
        field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,isz-1)][1] *= function(tmp);
       }
     }
    }

}

/*===================================
 FAMILY OF FILTERS, NGP, CIC, TSC  
====================================*/ 

void FilterFamily(fftw_complex* field, int P, std::string mode, double Radius){

   double fsize= Radius*double(NGRID)/double(BOXSIZE);
   double CUTOFF=1e-4;
   std::function<double(const double& )> function;
   if(mode.compare("Convolve")==0){
                function = [](const double&  x){return x;};
        }
   else if(mode.compare("Deconvolve")==0){
                function = [&r=Radius,&cutoff=CUTOFF](const double&  x){if(abs(x)>=cutoff) return 1/x;  
									else return 1/x; //turning off the cutoff option now... - this is dangerous?
									} ;
       }


   std::function<double(const double& )> powPsinc;
   powPsinc = [&p=P](const double&  x){
				return pow(sin(x)/x,double(p));
				} ;

   int ni,nf,icx,icy,icz,isx,isy,isz;
   float irz,iry,irx, zk, rks, tmp, cfilter=0.0;
   field[LOC(0,0,0)][0] *= 1.0; field[LOC(0,0,0)][1] *= 1.0;
   rks = fsize*2.0*M_PI/(double)NGRID;

   ni = 2;
   nf = NGRID/2 +1;
   tmp=1.0;


// os x
   icz =1;
   icy =1;

#pragma omp parallel for private(isx,irx,isy,isz,tmp) shared(function,powPsinc) num_threads(nthreads)
   for(icx=ni; icx<=nf; icx++)
    {
        isy = NGRID- icy +2;
        isz = NGRID - icz +2;
        isx = NGRID - icx +2;
        irx = (float)icx -1;

        irx*=rks;
	tmp=1.0; //dunno why, needed to initialize  tmp inside the loop???
        tmp *= function(powPsinc(irx));

     field[LOC(icx-1,icy-1,icz-1)][0] *= tmp; field[LOC(icx-1,icy-1,icz-1)][1] *= tmp;
     field[LOC(isx-1,icy-1,icz-1)][0] *= tmp; field[LOC(isx-1,icy-1,icz-1)][1] *= tmp;
    }

// os y
   icz =1;
   icx =1;

#pragma omp parallel for private(isy,iry,isz,isx,zk,tmp) shared(function) num_threads(nthreads)                                                                                             
   for(icy=ni; icy<=nf; icy++)
    {
        isy = NGRID- icy +2;
        iry = (float)icy-1.0;
        isz = NGRID - icz +2;
        isx = NGRID - icx +2;

        iry*=rks;
	tmp=1.0;
        tmp *= function(pow(sin(iry)/(iry),P));
        

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
    }

// os z
   icx =1;
   icy =1;

#pragma omp parallel for  private(isz,isy,isx,irz,zk,tmp) shared(function) num_threads(nthreads)                                                                                             
   for(icz=ni; icz<=nf; icz++)
    {
        isy = NGRID- icy +2;
        isz = NGRID - icz +2;
        irz = (float)icz-1.0;
        isx = NGRID - icx +2;

        irz*=rks;
        tmp=1.0;
        tmp *= function(pow(sin(irz)/(irz),P));
        

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
    }

// k_z =0 plane
  icz = 1;

#pragma omp parallel for collapse(2) private(isy,iry,isx,isz,irx,zk,tmp) shared(function) num_threads(nthreads)                                                                                             
  for(icy = ni; icy <=nf; icy++)
   {
    for(icx=ni; icx<=nf; icx++)
     {
        isy = NGRID- icy +2;
        iry = (float)icy-1.0;
        isz = NGRID - icz +2;
        isx = NGRID - icx +2;
        irx = (float)icx -1;

        irx*=rks;
        iry*=rks;
        tmp=1.0;

        tmp *= function(pow(sin(irx)/(irx),P));
        tmp *= function(pow(sin(iry)/(iry),P));
        

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
     }
   }

// k_y=0 plane
  icy = 1;

#pragma omp parallel for collapse(2) private(isx, irx, isz, irz, isy, zk, tmp) shared(function) num_threads(nthreads)                                                                                             
  for(icx = ni; icx <=nf; icx++)
   {
    for(icz=ni; icz<=nf; icz++)
     {

        isy = NGRID- icy +2;
        isz = NGRID - icz +2;
        irz = (float)icz-1.0;
        isx = NGRID - icx +2;
        irx = (float)icx -1;

        irx*=rks;
        irz*=rks;
        tmp=1.0;

        tmp *= function(pow(sin(irx)/(irx),P));
        tmp *= function(pow(sin(irz)/(irz),P));
        

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
     }
   }

// k_x=0 plane
  icx = 1;

#pragma omp parallel for collapse(2) private(isy,iry,isz,irz,isx,zk,tmp) shared(function) num_threads(nthreads)                                                                                             
  for(icy = ni; icy <=nf; icy++)
   {
    for(icz=ni; icz<=nf; icz++)
     {

        isy = NGRID- icy +2;
        iry = (float)icy-1.0;
        isz = NGRID - icz +2;
        irz = (float)icz-1.0;
        isx = NGRID - icx +2;

        iry*=rks;
        irz*=rks;
        tmp=1.0;

        tmp *= function(pow(sin(iry)/(iry),P));
        tmp *= function(pow(sin(irz)/(irz),P));
        

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
     }
   }

// the rest is volume-symetric

#pragma omp parallel for collapse(3) private(isz,irz,isy,iry,isx,irx,zk,tmp) shared(function) num_threads(nthreads)                                                                                             
  for(icz=ni; icz<=nf; icz++)
   {
    for(icy=ni; icy<=nf; icy++)
     {
      for(icx=ni; icx<=nf; icx++)
       {

	isy = NGRID- icy +2;
        iry = (float)icy-1.0;
	isz = NGRID - icz +2;
    	irz = (float)icz-1.0;
        isx = NGRID - icx +2;
        irx = (float)icx -1;

	irx*=rks;
	iry*=rks;
	irz*=rks;

        tmp=1.0;
        tmp *= function(pow(sin(irx)/(irx),P));
        tmp *= function(pow(sin(iry)/(iry),P));
	tmp *= function(pow(sin(irz)/(irz),P));
        


        field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,isz-1)][1] *= function(tmp);
       }
     }
    }

}



/*==============================
TOP HAT FILTER FOR K-SPACE
==============================*/


void TopHat(fftw_complex *field, std::string mode, double Radius)
{

   double fsize= Radius*double(NGRID)/double(BOXSIZE);
   std::function<double(const double& )> function;
   if(mode.compare("Convolve")==0){
                function = [&](const double&  x){return x;};
        }
   else if(mode.compare("Deconvolve")==0){
                function = [&](const double&  x){if(x!=0) return 1/x ;
                                                 else     return 1e6 ;} ;
        }

   int ni,nf,icx,icy,icz,isx,isy,isz;
   float irz,iry,irx, zk, rks, tmp, cfilter=0.0;

    field[LOC(0,0,0)][0] *= 1.0; field[LOC(0,0,0)][1] *= 1.0;
    rks = fsize*2.0*M_PI/(float)NGRID;
    ni = 2;
    nf = NGRID/2 +1;

// os x
   icz =1;
   icy =1;

//#pragma omp parallel for private(isx,zk,tmp) num_threads(nthreads)
   for(icx=ni; icx<=nf; icx++)
    {
     isx = NGRID - icx +2;
     zk = ((float)icx-1) * rks;
     tmp = 3.0*(sin(zk)-zk*cos(zk))/pow(zk,3.0);

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
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

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
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

     field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
     field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
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

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
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

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
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

      field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
      field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
      field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
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
  
        field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,isz-1)][1] *= function(tmp);
       }
     }
    }

}


/*==============================
GAUSSIAN FILTER FOR K-SPACE
==============================*/

void Gaussian(fftw_complex* field, std::string mode, double Radius)
{
   double fsize= Radius*double(NGRID)/double(BOXSIZE);
   std::function<double(const double& )> function;
   if(mode.compare("Convolve")==0){
		function = [&](const double&  x){return x;};
	}
   else if(mode.compare("Deconvolve")==0){
                function = [&](const double&  x){if(x!=0) return 1/x ;
						 else     return 1e6 ;} ;
	}

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

        field[LOC(icx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(icx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(icx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(icx-1,isy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,icz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,icz-1)][1] *= function(tmp);
        field[LOC(isx-1,icy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,icy-1,isz-1)][1] *= function(tmp);
        field[LOC(isx-1,isy-1,isz-1)][0] *= function(tmp); field[LOC(isx-1,isy-1,isz-1)][1] *= function(tmp);
       }
     }
   }
}



