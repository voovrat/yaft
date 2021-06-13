#ifndef YAFFT_H
#define YAFFT_H

#include <complex>

namespace yaft {
	
	void fft( const std::complex<double> * X 
		     , std::complex<double> * Y 
		     , int N
		     , int xstep = 1
		     , int ystep = 1);

   // dont forget to divide by N after conversion
	void ifft( const std::complex<double> * X 
		      , std::complex<double> * Y 
		      , int N
		      , int xstep = 1
		      , int ystep = 1);

  // i assume that the data lies row after row (c-style, not fortran or matlab-style)
	void fft2( const std::complex<double> * X
		      , std::complex<double> *Y  
		      , int nrow
		      , int ncol 
		      );

   // dont forget to divide by M*N after conversion
	void ifft2( const std::complex<double> *X
             , std::complex<double> * Y
             , int nrow
             , int ncol
		      );

}

#endif