#include "yafft.h"

#include <vector>

template<class Fn >
void xfft2( const std::complex<double> * X
               , std::complex<double> *Y  
		         , int nrow
		         , int ncol 
		         )
{
	std::vector< std::complex<double> > tmp(nrow*ncol);

	// do fft for all rows first
	const std::complex<double> *px = X;
	std::complex<double> *py = &tmp[0];
	for( int irow=0; irow<nrow; irow++, px += ncol, py += ncol )
	{
		Fn ( px, py, ncol, 1, 1 );
	}

	px = &tmp[0];
	py = Y;
	// and now for all cols 
	for( int icol=0; icol<ncol; icol++, px++, py++ )
	{
		Fn ( px, py, nrow, ncol, ncol );
	}

}


struct fftFn 
{
	 fftFn( const std::complex<double> * X
               , std::complex<double> *Y  
		         , int N
		         , int xstep
		         , int ystep 
		         )  
	 {  yaft::fft(X,Y,N,xstep,ystep); }

};


struct ifftFn 
{
	 ifftFn( const std::complex<double> * X
               , std::complex<double> *Y  
		         , int N
		         , int xstep
		         , int ystep 
		         )  
	 {  yaft::ifft(X,Y,N,xstep,ystep); }

};



void yaft::fft2( const std::complex<double> *X
                , std::complex<double> * Y
                , int nrow
                , int ncol
		          )
{
	xfft2<fftFn>( X,Y,nrow, ncol );

}


void yaft::ifft2( const std::complex<double> *X
                , std::complex<double> * Y
                , int nrow
                , int ncol
		          )
{
	xfft2<ifftFn>( X, Y, nrow, ncol );
}
