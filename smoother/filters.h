#include <fftw3.h>
#include <functional>
#include <iostream>
#include <omp.h>

void FilterInterface(fftw_complex *field, std::string FilterConvolveName, std::string FilterDeconvolveName, std::string TypeConvolve, std::string TypeDeconvolve, double RadiusConvolve, double RadiusDeconvolve);

void FilterFamilyRadial(fftw_complex* field, int P, std::string mode, double Radius); //P determines the specific radial-generalized form of the filter, NGP, CIC, TSC;
void FilterFamily(fftw_complex* field, int P, std::string mode, double Radius); //P determines which generalized-by-multiplication filter to use,  NGP, TSC, CIC 
void Gaussian(fftw_complex* field, std::string mode, double Radius);
void TopHat(fftw_complex *field, std::string mode, double Radius);


inline int nthreads=omp_get_max_threads();
inline int NGRID =512;
inline double BOXSIZE=500;
// TopHat filter of Radius in real space is given, in spherical coordinates, by: W_R(r,theta,phi) = 1/(4/3 * Pi * R^3) if r<=R, and 0 otherwise
// performing a 3-dimensional Fourier transform, the resulting function is (FW_R)(k,theta_k,phi_k) = (FW_R)(k) = 3*j1(kR)/kR = 3* (sin(kR) - kR*cos(kR)) / (kR)^3

// Gaussian filter is G_R(r,theta,phi) = G_R(r) = (8*Pi^3*(R^3)^2)^(-1/2) * exp(-0.5*(r/R)^2) and in Fourier space, (FG_R)(k)= 0.5 exp(-(kR)^2 / 2),
// the integral in real space is normalized to 1 over cartesian product R^3, and normalized so that G(0)/G(R) = sqrt(e), and similarly G(0)/G(sqrt(2)*R)) = e 


//FilterFamilyRadial is a family of filters given by 3 functions, generalized from 1-dimensional cases, namely:
// NGP, p=1 is: W_R(r) = 1/((4/3)*M_PI*R^3) if r<=R and 0 otherwise
// in Fourier space: (FW_R)(k) =  3*j1(kR)/kR = 3* (sin(kR) - kR*cos(kR)) / (kR)^3
// CIC, p=2 is W_R(r) = 4/((4/3)*M_PI*R^3) * (1-r/R) if r<=R 
// in Fourier space (FW_R)(k) = 12 * ( 2 - 2*cos(kR) - kR*sin(kR) ) / (kR)^4
// TSC, p=3 is W_R(r) = { 1/3 - (r/R)^2 if r<=1/3 * R,  0.5*(1-(r/R))^2 if 1/3 R <= r <= R and 0 otherwise
// in Fourier space (FW_R)(k) = 243/4 * (9sin(kR/3) - 3sin(kR) + kRcos(kR) - kRcos(kR/3)) / (kR)^5

// all the above W_R integrate to 1 if integrated over a ball of radius R - hence the name filter of Radius=R

// FilterFamily is a generalization of 1-dimensional filters to a 3-dimensional situation by virtue of multiplying all the component filters, i.e. W_R(x,y,z) = Multiply_i W_R_xi(xi)
// as such, the support of such a filter is on a cube, not a sphere!
