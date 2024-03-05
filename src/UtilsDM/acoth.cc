#include <complex>


std::complex<double> acoth(std::complex<double> z){
  std::complex<double> zero(0.,0.);
  std::complex<double> imag(0.,1.);
  if (z==zero){
    return -imag*M_PI*(1./2);
  }
  else{
    return std::atanh(1./z);
  }
}
