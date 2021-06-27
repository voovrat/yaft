#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
//#include <math.h>
#include <regex>
#include <functional>
#include <string>
#include <stdio.h>
#include <fstream>
#include <sstream>

#include <prime_numbers.h>

#ifndef M_PI
  #define M_PI            3.14159265358979323846
#endif


constexpr bool DEBUG = false;

typedef std::complex<double> cmplx;


std::string cmplx2str( const cmplx & c )
{
	char s[256];
	sprintf(s,"cmplx(%0.20lf,%0.20lf)",c.real(), c.imag());
	return s;
}


cmplx operator*( const cmplx & a, int b) 
{
	return a * (double)b;
}

cmplx operator*( int b, const cmplx & c)
{
	return (double)b * c;
}

cmplx operator/( const cmplx & a, int b) 
{
	return a / (double)b;
}

cmplx operator/( int b, const cmplx & c)
{
	return (double)b / c;
}


void dft( const cmplx * X, cmplx *Y, int N, int xstep, int ystep )
{
	std::vector<cmplx> w(N);
	for( int m=0; m<N; m++ ) 
	{
		w[m] = exp(-2*M_PI*cmplx(0,1)*(double)m/(double)N);
	}

	// y0 = sum(xi)
	const cmplx *px = X;
	*Y = 0;
	for( int n=0; n<N; n++, px += xstep )
	{
		*Y += *px;
	}

	// yk = sum_n exp( -2pi I/N kn) xn = exp( -2 pi I * (kn % N )) xn = w[ kn % N ] * xn
	Y += ystep;
	for( int k=1; k<N; k++, Y += ystep )
	{
	   *Y = *X;  //  w0 = 1 always
	   px = X + xstep;
		for( int n=1; n<N; n++, px += xstep )
		{
			*Y += w[ k*n % N ] * *px;
		}
	}

}


void idft( const cmplx * X, cmplx *Y, int N, int xstep, int ystep )
{
	std::vector<cmplx> w(N);
	for( int m=0; m<N; m++ ) 
	{
		w[m] = exp(2*M_PI*cmplx(0,1)*(double)m/(double)N);
	}

	// y0 = sum(xi)
	const cmplx *px = X;
	*Y = 0;
	for( int n=0; n<N; n++, px += xstep )
	{
		*Y += *px;
	}

	// yk = sum_n exp( -2pi I/N kn) xn = exp( -2 pi I * (kn % N )) xn = w[ kn % N ] * xn
	Y += ystep;
	for( int k=1; k<N; k++, Y += ystep )
	{
	   *Y = *X;  //  w0 = 1 always
	   px = X + xstep;
		for( int n=1; n<N; n++, px += xstep )
		{
			*Y += w[ k*n % N ] * *px;
		}
	}
}


// get generator of multiplicative group Fp
// if p is not prime and there is no generator - returns 0 
int get_generator(int p )
{
	for( int g=2; g<=p-1; g++ )
	{
	   int gn = g; // gn = g^n
  	   bool cycle = false;  
  	   for( int i=1; i<=p-3; i++ )
  	   {
    		gn = ( g * gn ) % p ;
    	   if( gn == 1 )  // cycle found 
    	   {
    	   	cycle = true;
    	   	break;
    	   }
  	   }
  		if( !cycle ) return g;
   }
   // no generator found 
	return 0;
}


// for N<=4 we define fft transformations explicitly 

void genfftN1( std::ostream & os )
{

	os << 
R"EOL(

inline void fftN1( const cmplx * X, cmplx *Y, int unused_xstep, int unused_ystep )
{
	*Y = *X;
}

inline void ifftN1(const cmplx * X, cmplx *Y, int unused_xstep, int unused_ystep )
{
	*Y = *X;
}

)EOL";

}


void genfftN2( std::ostream & os )
{
	os << 
R"EOL(

// yk = sum x0 exp(-2pi/2 * 0)
inline void fftN2( const cmplx *X, cmplx *Y, int xstep, int ystep )
{
	const cmplx x1 = X[xstep];
	Y[0]     = *X + x1;
	Y[ystep] = *X - x1; 
   return;
}

inline void ifftN2( const cmplx * X, cmplx *Y, int xstep, int ystep ) { return fftN2(X,Y,xstep,ystep); }

)EOL";


}


void genfftN3( std::ostream & os )
{
	os << 
R"EOL(

	// yk = sum n=0..2 xn exp( -2pi i/3 nk)  =  x0 + x1*exp(-2pi/3 k) + x2*exp(-4piI/3 k )

   //  U = exp(-2p/3)
	// exp( -4pi/3) = exp(-4pi/3+6pi/3) = exp( 2pi/3 ) = V
	// y0  =  x0 + x1 + x2
	// y1  =  x0 + x1 exp( -2pi/3) + x2*exp(-4pi/3)    = x0 + U*x1 + V*x2  
	// exp( -8pi/3 ) = exp( -2pi/3) = U
   // y2  =  x0 + x1 exp( -4pi/3) + x2* exp ( -8pi/3) = x0 + V*x1 + U*x2
inline void fftN3( const cmplx *X
					  , cmplx * Y
					  , int xstep
					  , int ystep )
{
		static const cmplx U = exp(-2*M_PI/3*cmplx(0,1));
		static const cmplx V = exp(+2*M_PI/3*cmplx(0,1));
		
		const cmplx x0 = *X;
		X += xstep;
		const cmplx x1 = *X;
		X += xstep;
		const cmplx x2 = *X;

		*Y              = x0  +   x1  +   x2;
		 Y += ystep;
		*Y              = x0  + U*x1  + V*x2;
		 Y += ystep;
		*Y              = x0 +  V*x1  + U*x2;
}

inline void ifftN3( const cmplx *X
					   , cmplx * Y
					   , int xstep
					   , int ystep )
{
		static const cmplx U = exp(-2*M_PI/3*cmplx(0,1));
		static const cmplx V = exp(+2*M_PI/3*cmplx(0,1));
		
		const cmplx x0 = *X;
		X += xstep;
		const cmplx x1 = *X;
		X += xstep;
		const cmplx x2 = *X;

		*Y              = x0  +   x1  +   x2;
		 Y += ystep;
		*Y              = x0  + V*x1  + U*x2;
		 Y += ystep;
		*Y              = x0 +  U*x1  + V*x2;
}
)EOL";
}


void genfftN4( std::ostream & os )
{
	os << 
R"EOL(

  // yk = sum n=0..3  xn exp( -2pi/4 nk ) = sum xn exp( -pi/2 nk)
   // y0 = x0 + x1 + x2 + x3
   // y1 = x0 + x1 exp(-pi/2) + x2 exp(-pi)  + x3 exp( -3pi/2) 
   //    = x0 - x1 * i    - x2  + x4 * i
   // y2 = sum n=0..3 xn exp( -2pi/4*2 n) = sum xn exp( -pi n) = x0 - x1 + x2 - x3
   // y3 = x0 + x1 exp(-3pi/2) + x2 exp(-6pi/4) + x3 exp( -9pi/2)
   //    = x0 + x1 * i - x2  - x3*i 
inline void fftN4( const cmplx *X
					  , cmplx * Y
					  , int xstep
					  , int ystep ) 
{
	const cmplx x0 = *X;
	X += xstep;
	const cmplx x1 = *X;
	X += xstep;
	const cmplx x2 = *X;
   X += xstep;
	const cmplx x3 = *X;

	const cmplx x0_plus_x2 = x0 + x2;
	const cmplx x0_minus_x2 = x0 - x2;
	const cmplx x1_plus_x3  = x1 + x3;
	const cmplx x1_minus_x3 = x1 - x3;
	//  cmplx x1_minus_x3_I = x1_minus_x3*cmplx(0,1);
   // (a + bi)*i = ai - b 
   const cmplx x1_minus_x3_I( -x1_minus_x3.imag(), x1_minus_x3.real());

	*Y      = x0_plus_x2  + x1_plus_x3;
	 Y += ystep;
	*Y      = x0_minus_x2 - x1_minus_x3_I ;
	 Y += ystep;
	*Y      = x0_plus_x2  - x1_plus_x3;
	 Y += ystep;
	*Y      = x0_minus_x2 + x1_minus_x3_I;  
}

inline void ifftN4( const cmplx *X
					   , cmplx * Y
					   , int xstep
					   , int ystep ) 
{ 
	const cmplx x0 = *X;
	X += xstep;
	const cmplx x1 = *X;
	X += xstep;
	const cmplx x2 = *X;
   X += xstep;
	const cmplx x3 = *X;

	const cmplx x0_plus_x2 = x0 + x2;
	const cmplx x0_minus_x2 = x0 - x2;
	const cmplx x1_plus_x3  = x1 + x3;
	const cmplx x1_minus_x3 = x1 - x3;
	//  cmplx x1_minus_x3_I = x1_minus_x3*cmplx(0,1);
   // (a + bi)*i = ai - b 
   const cmplx x1_minus_x3_I( -x1_minus_x3.imag(), x1_minus_x3.real());

	*Y      = x0_plus_x2  + x1_plus_x3;
	 Y += ystep;
	*Y      = x0_minus_x2 + x1_minus_x3_I ;
	 Y += ystep;
	*Y      = x0_plus_x2  - x1_plus_x3;
	 Y += ystep;
	*Y      = x0_minus_x2 - x1_minus_x3_I;  
}

)EOL";

}





void gensignature( std::ostream & os, int N, int sign, char tail = '{')
{
	const char * name = ( (sign == -1) ? "fft" : "ifft" );
	os << " void " << name << "N" << N << "( const cmplx * X, cmplx *Y, int xstep, int ystep ) " << tail << " \n";	
}


// discrete fourier transform with O(N^2) implementation "by definition"
// for small N (especially for prime numbers) it might still be faster then the rest 

void genDFTX( std::ostream & os, int N, int sign )
{
	gensignature( os, N, sign );

	std::vector<cmplx> w(N);
	for( int m=0; m<N; m++ ) 
	{
		w[m] = exp( sign * 2*M_PI*cmplx(0,1)*(double)m/(double)N);
	}


	for( int i=0; i<N; i++)
	{
		os << " const cmplx x" << i << " = *X;";
		os << ( (i==N-1) ? "\n" : "X+=xstep;\n");
	}

   if(DEBUG) os << "std::cout << \" Hello, I am direct DFT algorithm with N = " << N << "\\n\";\n";

	os << "// y0 = sum(xi)\n";
	
	os << "*Y = x0 ";
	for( int i=1; i<N; i++)
	{
		os << " +x" << i; 
	}
	os << ";\n";


	os << "// yk = sum_n exp( -+2pi I/N kn) xn = exp( -+2 pi I * (kn % N )) xn = w[ kn % N ] * xn\n";

	for( int k=1; k<N; k++ )
	{
	   os << "Y += ystep;\n";
	   os << "*Y = x0 ";
		for( int n=1; n<N; n++ )
		{
			//*Y += w[ k*n % N ] * *px;
			cmplx C = w[ k*n % N];
			
			os << "+" << cmplx2str(C) << " * x" << n << "\n" ;
		}
		os << ";\n";
	}

os << "}\n\n";

}


void genDFT( std::ostream & os, int N, bool geninline)
{
	if( geninline ) os << "inline ";
	genDFTX( os, N, -1 );
	if( geninline ) os << "inline ";
	genDFTX( os, N, 1  );
}



std::vector<cmplx> genWeights( std::ostream & os, int N, int count, int sign )
{

	std::vector<cmplx> w( count );

	std::string coma = "  ";
	os << " static const cmplx w[" << count << "] = {"; 
	for( int n=0; n< count ; n++ )
	{
		w[n] = exp( sign * 2*M_PI*cmplx(0,1)*n/N); 
		os << coma << cmplx2str( w[n] ) << "\n";
		coma = ", ";
	}
	os << "};\n";


	os << "\n";	

	return w;
}





inline std::string fftName( int sign, int N )
{
	std::string prefix = (sign == -1 ? "fftN" : "ifftN" ); 
	return prefix + std::to_string(N);
}

//  Cooley-Tukey algorithm for N = 2n
void genfft2nX( std::ostream & os, int N, int sign, bool inlinefor )
{
	gensignature( os, N, sign );
	std::vector<cmplx> w = genWeights( os, N, N/2, sign );
   
   os << "cmplx  E[" << N/2 << "];\n";  
   os << "cmplx  O[" << N/2 << "];\n";  

   os << "\n";

   int N2 = N /2;
   os << "int xstep2 = xstep<<1;\n";

   os << "\n";

   if( DEBUG ) os << "std::cout << \" Hello, I am FFT N2 with N = " << N << "\\n\";\n";

   os << fftName(sign, N2) << "(X       ,E, xstep2, 1 );\n";
   if( DEBUG) os << "std::cout << \" Hello, I am FFT N2 with N = " << N << "  : After first fft \\n\";\n";   

   os << fftName(sign, N2) << "(X+xstep ,O, xstep2, 1 );\n"; 

   if(DEBUG) os << "std::cout << \" Hello, I am FFT N2 with N = " << N << "  : After second fft \\n\";\n";   

   os << "\n";

   os << "cmplx * py = Y;\n";
   os << "cmplx * py1 = Y + " << N2 << "*ystep;\n";

   if( inlinefor )
   {
   	
	   for( int i=0; i<N2;  i++ )
  		{
  			if( i==0) os << "cmplx ";
   		os << "wO = ";
   		if( i>0 	) os <<  "w[" << i << "]* ";

      	os << "O[" << i << "]; \n"; 

   		os << "*py   = E[" << i << "] + wO;\n";
   		os << "*py1  = E[" << i << "] - wO;\n";
   		
   		if( i<N2-1)
   			os << "py+=ystep; py1+=ystep;\n";
   	}
   } else
   {
	   os << 
"for( int i=0; i<" << N2 << ";  py+=ystep, py1+=ystep, i++ )\n";
	   os << 
R"EOL(
{
   cmplx wO = w[i] * O[i];
   *py   = E[i] + wO;
   *py1  = E[i] - wO;
}
)EOL";
   }

 os << "} \n\n";
}


void genfft2n( std::ostream &os, int N, bool geninline )
{
	bool inlinefor = (N<=128);

	if( geninline ) os << "inline ";
	genfft2nX( os, N, -1, inlinefor );
	if( geninline ) os << "inline ";
	genfft2nX( os, N,  1,  inlinefor );
}

std::string make_index( int i, const std::string & ystep, const std::string & plus )
{
	if( i == 0 )
		return plus == "" ? "0" : "";
	else if( i==1)
		return plus + " " + ystep;

	int two_n = 2;
	for( int n=1, two_n=2; n<16; n++, two_n*=2 )
	{
		if( i==two_n)
		{
			if( plus != "" )
				return plus + " (" + ystep + " << " + std::to_string(n) + ")";
			else
				return  ystep + " << " + std::to_string(n);
		}
	}

	return plus + std::to_string(i) + " * " + ystep ;
}


// General Cooley-Tukey algorithm for N = d*N1
void genCooleyTukeyX( std::ostream & os, int d, int N1, int sign, bool inlinefor )  
{
	int N = d*N1;
	gensignature( os, d * N1, sign );
	std::vector<cmplx> w = genWeights( os, N, N, sign );

	os << "\n";

	for( int i=0; i<d; i++ )
	{
		os << "cmplx E" << i << "[" << N1 << "];\n";		
	}

	os << "\n";

   if(DEBUG) os << "std::cout << \" Hello, I am Cooley Tukey with N = " << d*N1 << " \\n\";\n"; 

	os << "int dxstep = " << d << "*xstep;\n";
	os << "const cmplx *px = X;\n";

	for( int i=0; i<d; i++ ) 
	{
		//fft( X + i*xstep,  E[i],  Y, w, N/d, xstep*d, ystep*d );
		os << fftName(sign, N1) << "( px,  E" << i << ", dxstep, 1 );\n";

     if(DEBUG) os << "std::cout << \" Hello, I am Cooley Tukey with N = " << d*N1 << " : After " << (i+1) << "th FFT \\n\";\n"; 

		if( i < d-1 )
			os << "px += xstep;\n";
	}

	os << "\n";

	if( inlinefor )
	{

	//	os << "cmplx * py_start = Y;\n";
		os << "int big_ystep = " << make_index( N1 , " ystep","") << ";\n";
		os << "int jump_back = ( " << make_index(d-1,"big_ystep","") << " ) - ystep;\n";
		os << "cmplx * py = Y;\n\n";
		for( int i=0; i<N1; i++ )
		{
	//		if( i>0)
	//			os << "py = py_start;\n";

			for( int r = 0; r<d; r++)
			{
				int ii = (i + r*N/d);
			

				//os << "Y[" << make_index(ii,"ystep","")  << "] = E0[" << i << "]"; 
				//Y[kk]  = E[0][k*d];
				os << "*py = E0[" << i << "]";

				for( int q=1; q<d; q++ )
				{
					//Y[kk] += w[ ((ii*q)%N)*ystep ] * E[q][k*d];
					int idx = ((ii*q)%N);

					if( ii == 0 || std::abs(w[idx]-1.0) < 1e-14)
						os << " + E" << q << "[" << i << "]/*w" << idx << "*/"; 
					else if ( std::abs( w[idx]+1.0) < 1e-14 )
						os << " - E" << q << "[" << i << "]/*w" << idx << "*/"; 
					else
						os << " + w[" << idx << "]*E" << q << "[" << i << "]";

					if ( q % 20 == 0 ) os << "\n";
				}
				os << ";\n";

				if( r < d-1 )
					os << "py += big_ystep;\n";

			}// for r

			if( i < N1-1)
				os << "py -= jump_back;\n\n";
		//		os << "py_start += ystep;\n\n";

		} // for i 

	}
	else
	{
		os << "cmplx * py_start = Y;\n";
		os << "int big_ystep = " << make_index( N1 , "ystep","") << ";\n";
		 
		os << "for( int i=0; i<" << N1 << "; i++, py_start += ystep  ) \n";
		os << "{\n";

		os << " int ii = i;\n";
		os << " cmplx * py = py_start;\n";

		for( int r = 0; r<d; r++)
		{
			os << " *py = E0[i] ";

			for( int q=1; q<d; q++ )
			{
			//	Y[kk] += w[ ((ii*q)%N)*ystep ] * E[q][k*d];

				if( q == 1 )
					os << " + w[ " << make_index(q, "ii","") << " ] * E" << q << "[i]";
				else 
					os << " + w[ (" << make_index(q, "ii","") << ") % " << N << " ] * E" << q << "[i]";

				if ( q % 20 == 0 ) os << "\n";
			} // for q 

			os << ";\n";
			if( r <  d-1 )
			{
				os << " py += big_ystep;     ii += " << N1 << ";\n";
			}

		}// for r 


		os << "} // for i \n";
	}

os << "}\n\n";
}


void genCooleyTukey( std::ostream & os, int d, int N1, bool geninline) 
{
	bool inlinefor = ( N1<=128 && d < 10 );
	if( geninline ) os << "inline ";
	genCooleyTukeyX( os, d ,N1, -1, inlinefor );
	if( geninline ) os << "inline ";
	genCooleyTukeyX( os, d ,N1,  1, inlinefor );
}

// Rader's algorithm for prime N=P
void genRaderX( std::ostream & os,  int P, int sign, bool inlinefor, bool inlineAssignment )
{
	gensignature( os, P, sign );

	std::vector<int> IY(P-1);
	std::vector<cmplx> C(P-1);

	int g = get_generator(P);

	int gn = g;
   for( int i=0; i<P-1; i++)
	{
	  	IY[gn-1] = i;
		C[i] = exp( sign * 2*M_PI*cmplx(0,1)/P*gn ) / (P-1);  
		gn = (g * gn) % P;
	}

  	std::vector<cmplx> FC(P-1);
  	dft( &C[0], &FC[0], P-1, 1, 1);

  	std::string coma = " ";

	os << " static const cmplx FC[" << P-1 << "] = {"; 
	for( int n=0; n< P-1 ; n++ )
	{
		//w[n] = exp( sign * 2*M_PI*cmplx(0,1)*n/N); 
		os << coma << cmplx2str( FC[n] ) << "\n";
		coma = ", ";
	}
	os << "};\n";


	os << "\n";	


	os << "//	y0 = sum(x);\n";

	
   if(DEBUG) os << "std::cout << \" Hello, I am Rader's algorithm with N = " << P << " \\n\";\n"; 


	if( ! inlineAssignment )
	{
		os << "*Y = *X;\n";
		os << "const cmplx *px = X + xstep;\n";
		os << "for ( int i=1; i<" << P << "; i++, px += xstep)  *Y += *px;\n";
	}
	else
	{
		os << " const cmplx * px = X;\n";
		os << "*Y = *px;  px += xstep; \n";
		for( int i=1; i<P; i++)
		{
			os << "*Y += *px;";
			if( i<P-1)
				os << "   px += xstep;";
			os << "\n";
		}

	}

   if(DEBUG) os << "std::cout << \" Hello, I am Rader's algorithm with N = " << P << " : After y0=sum(xn)\\n\";\n";

	os << "\n";

	os << "//	xx = x(2:P);\n";
	os << "const cmplx *xx = X + xstep;\n";

	os << "\n";

//	std::vector<int> IX(P-1);


	if( !inlinefor ) os << "int IY[" << P-1 << "];\n";  


	std::vector<cmplx> XX(P-1);

	os << "cmplx XX[" << P-1 << "];\n";

	if( inlinefor )
	{
		int gn = g;
	   for( int i=0; i<P-1; i++)
	   {
	  		os << "XX[" << P-i-2 << "] = xx[ " << make_index( gn-1, "xstep", "") << " ];\n";
	  		gn = (g * gn) % P;
	  	}
   }
   else
   {
   	os << "int gn = " << g <<";\n";

	   os << "for( int i=0; i<" << P-1 << "; i++)\n";
	   os << "{\n";

	   os << "  IY[gn-1] = i;\n";
	  	os << "  XX[" << P-2 << " - i ] = xx[ (gn-1) * xstep ];\n";
	  	os << "  gn = ( " << make_index(g , "gn", "" ) << " ) % " << P << ";\n";

	  	os << "}\n";

   }

   if(DEBUG) os << "std::cout << \" Hello, I am Rader's algorithm with N = " << P << " : After initialization of XX\\n\";\n";

   os << "\n";

 	//std::vector<cmplx> FXX(P-1);
 	os << "cmplx FXX[" << P-1 << "];\n";

  	//dft( &XX[0], &FXX[0], P-1, ystep, 1);
  	os << fftName( -1, P-1) << "( XX, FXX, 1, 1 );\n";

    if(DEBUG) os << "std::cout << \" Hello, I am Rader's algorithm with N = " << P << " : After first fourier\\n\";\n";

  	if( !inlineAssignment )
  	{
  		os <<  "for( int i=0; i<" << P-1 << "; i++ ) FXX[i] *= FC[i];\n";
   }
   else
   {
   	os << "\n";
   	for( int i=0; i<P-1; i++ )
   		os << "FXX[" << i << "] *= FC[ " << i << "];\n";
   	os << "\n";
   }

  	// cmplx * YY = &XX[0];

  	 //idft( &FXX[0], YY, P-1, 1, 1 );
  	os << fftName( 1, P-1 ) << "(FXX, XX, 1, 1);\n";
  	os << " // so now in XX we have convolution ifft( fft(XX) .* FC );\n";

   if(DEBUG) os << "std::cout << \" Hello, I am Rader's algorithm with N = " << P << " : After inverse fourier\\n\";\n";

  	os << "\n";

  	os << "Y += ystep;\n";

  	if( !inlinefor )
  	{
   	os << "for( int i=0; i<" << P-1 << "; i++, Y+=ystep )\n";
   	os << "  *Y = XX[ IY[i] ] + *X;\n";
   }
   else
   {
     	for( int i=0; i<P-1; i++ )
   	{
   		os << " *Y = XX[ " << IY[i] << " ] + *X;";
   		if( i < P-2 )
   			os << "    Y += ystep;";
   		os << "\n";
   	}
   }

   if(DEBUG) os << "std::cout << \" Hello, I am Rader's algorithm with N = " << P << " : After assignment to Y\\n\";\n";

   os << "}\n\n";
}


void genRader( std::ostream & os, int P, bool geninline)
{
	bool inlineAssignment = (P<20);
	bool inlinefor = (P<100);

	if( geninline ) os << "inline ";
   genRaderX( os,  P, -1, inlinefor, inlineAssignment );
	if( geninline ) os << "inline ";
   genRaderX( os,  P,  1, inlinefor, inlineAssignment );
}


void genFFT( std::ostream & os, int N, bool geninline  )
{
	if( N == 1 )
	{
		genfftN1( os ); 
		return;
	}

	if( N == 2 )
	{
		genfftN2( os  ); 
		return;
	}

	if( N == 3 )
	{
		genfftN3( os  ); 
		return;
	}

	if( N == 4 )
	{
		genfftN4( os  ); 
		return;
	}

	if( (N == 5) || (N == 7) )
	{
		os << "// N = " << N << " : dft by definition\n";
		genDFT( os, N, geninline  ); 
		return;
	}

	// if N = 2k we call Cooley–Tukey algorithm for N = 2n 
	if( (N & 1) == 0)
	{

		os << "// N = " << N << " fft2n \n";
		genfft2n( os, N, geninline );
		return;
	}

	// otherwise we find minimal prime divisor of N 
	int SQN = (int)sqrt(N);

	const std::vector<int> & primes = TPrimeNumbers::PrimeList();
	int d = 1;
	for( int p : primes )
	{
		if( p > SQN ) break;

		if( N % p == 0 ) 
		{
			d = p; break;
		}
	}

 
	if( d > 1 )  
	{
		os << " // N = " << N << " d = " << d << ", N/d = " << N/d << " :  Cooley-Tukey \n";
		// for composite N we generate general Cooley-Tukey fft 
		genCooleyTukey( os, d, N/d, geninline  );
	}
	else 
	{
		os << " // N = " << N << " Rader \n";
		// for prime N we generate Rader's fft 
		genRader( os, N, geninline  );
	}

}

void genAllFFT( std::ostream & os, int N)
{
	for(int n=1; n<=N; n++)
		genFFT( os, n, false); 

}





void genSimple( std::ostream & os)
{
	genfftN1(os);
	genfftN2(os);
	genfftN3(os);
	genfftN4(os);
}



void genTest( std::ostream & os, const std::function< void( std::ostream & os ) > & genfft, int N ) 
{

	os << 
R"EOL(
#include <vector>
#include <complex>
#include <iostream>
#include <stdio.h>
#include <string>

typedef std::complex<double> cmplx;

std::string cmplx2str( const cmplx & c )
{
	char s[256];
	sprintf(s,"(%10.5lf,%10.5lf)",c.real(), c.imag());
	return s;
}


)EOL";

	genfft( os);


	std::string text = 
R"EOL(

main()
{
	int N = @;
	std::vector<cmplx> X(N);
	std::vector<cmplx> Y(N);
	std::vector<cmplx> Z(N);

	int i=1;
	for( cmplx & z: X ) z = i++;

	fftN@(&X[0], &Y[0], 1, 1);
	ifftN@( &Y[0], &Z[0], 1, 1);

	for( int i=0; i<N; i++ )
	{
		std::cout << cmplx2str( X[i] ) << "  " 
		          << cmplx2str( Y[i] ) << "  " 
		          << cmplx2str( Z[i] ) << "  " 
		          << cmplx2str( X[i] - Z[i]/(double)N )
		          << "\n";
	}

}

)EOL";

	os << std::regex_replace(text, std::regex("@"), std::to_string(N));

}


void genTestList( std::ostream & os, int N )
{

	os << 
R"EOL(
#include <iostream>
#include <vector>

#include "yaft.h"


)EOL";



	std::string text = 
R"EOL(

double do_test@()
{
	int N = @;
	std::vector<cmplx> X(N);
	std::vector<cmplx> Y(N);
	std::vector<cmplx> Z(N);

	int i=1;
	for( cmplx & z: X ) z = i++;

	fftN@(&X[0], &Y[0], 1, 1);
	ifftN@( &Y[0], &Z[0], 1, 1);

	double S = 0;
	for( int i=0; i<N; i++ )
	{
		S += std::norm( X[i] - Z[i]/(double)N );
	}

	return sqrt(S)/(double)N;
}

)EOL";


for( int n=1; n<=N; n++)
{
	os << std::regex_replace(text, std::regex("@"), std::to_string(n));
}

os << "main() { \n";

for( int n=1; n<=N; n++)
{
	os << "std::cout << " << n << " << \"  \" << do_test" << n << "() << \"\\n\"; \n"; 
}

os << "}\n";

}



void test()
{
   genTest(std::cout
   		 , [](std::ostream & os) { 
   		 	//	genSimple(os);
   		 	//	genfft2n(os, 8); 
   		 	//	genfft2n(os,16);
   		 	//	genCooleyTukey( os, 97, 97 );
   		 		genAllFFT(os, 2048 );
   		 	//	genRaderX( os, 7, -1, true, true );
   		 	//	genRaderX( os, 7, 1, true, false );
   		 }
   		 , 121
   		 );
//	genTest( std::cout, [](std::ostream &os) { genfftN4(os);},4 );
//	genDFT( std::cout, 5);

}


void createHeader( int nmax, int ninline )
{
	std::ofstream fs("yaft.h");

	fs << 
R"EOL(#ifndef YAFT_H
#define YAFT_H

#include <complex>


typedef std::complex<double> cmplx;

)EOL";

	if( DEBUG ) fs << "#include <iostream> \n\n";

	for( int N = 1; N <= nmax; N++ )
	{
		if( N >= ninline )
		{
			gensignature( fs, N, -1, ';' );
			gensignature( fs, N,  1, ';' );
		}
		else
		{
			genFFT( fs, N, true);
		}
		fs << "\n";		

	}

	fs << "#endif\n";
}




void createCPP( int N )
{
	std::ofstream fs( "yaft_fftN" + std::to_string(N)  + ".cpp");

	fs << "#include \"yaft.h\"\n\n";
	if(DEBUG) fs << "#include <iostream>\n\n";

	genFFT(fs, N, false);
}



void createInterface( int Ninline, int NMAX )
{
	std::ofstream fs("yafft.cpp");

	fs <<
R"EOL(
#include "yafft.h"
#include "yaft.h"


typedef void (*FnFFTN)( const cmplx *X, cmplx *Y, int xstep, int ystep);

)EOL";


for( int N = 1; N<Ninline; N++)
{
	fs << "void runfftN" << N << "( const cmplx *X, cmplx *Y, int xstep, int ystep) { fftN" << N << "(X,Y,xstep,ystep); }\n";
	fs << "void runifftN" << N << "( const cmplx *X, cmplx *Y, int xstep, int ystep) { ifftN" << N << "(X,Y,xstep,ystep); }\n";
	fs << "\n";
}

fs << "\n";

for( int sign = -1; sign < 2; sign+=2)
{
	std::string fftName = (sign == -1 ? "fft" : "ifft");

	fs << "void yaft::" << fftName << R"EOL( ( const std::complex<double> * X 
		     , std::complex<double> * Y 
		     , int N
		     , int xstep
		     , int ystep )
{
)EOL";

	fs << "const static FnFFTN  fn[" << NMAX << "] = {\n";
	std::string coma = "  ";

	for( int N = 1; N<Ninline; N++)
	{
		fs << coma << "&run" << fftName << "N" << N << "\n";
		coma = ", ";
	}

	for( int N = Ninline; N<=NMAX; N++)
	{
		fs << ", &" << fftName << "N" << N << "\n";
   }

	fs << "  };\n";


	fs << "  if( N > " << NMAX << ") throw \"fft not implemented for N>" << NMAX << ". If you need it you may recompile yaft with larger NMAX.\"; ";
	fs << "\n";
	fs << "  fn[N-1](X,Y,xstep,ystep);\n"; 

	fs << "}\n\n\n";

}

}

main()
{
	constexpr int Ninline = 5;
	constexpr int  NMAX = 2048;

	std::ofstream fs( "make_yaft.sh");

	std::stringstream ss;


	fs << "rm *.o\n";

	createHeader( NMAX, Ninline );
	for( int N = Ninline; N<=NMAX; N++)
	{
		std::cout << "yaft_fftN" << N << "\n";
		createCPP(N);
		fs << "echo yaft_fftN" << N << "\n";
		fs << "g++ -O3 -c yaft_fftN" << N << ".cpp \n";
//		ss << " yaft_fftN" << N << ".o";
	}

	fs << "g++ -O3 -c yafft.cpp\n";
	fs << "g++ -O3 -c yafft2.cpp\n";

	fs << "rm libyaft.a\n";
	fs << "ar -rcs libyaft.a *.o\n";
	fs << "g++ -o test_yaft test_yaft.cpp -g libyaft.a --std=c++14\n";


	createInterface(Ninline,  NMAX);


//	fs << ss.str() << "\n";

   std::ofstream ftest("test_yaft.cpp") ;
	genTestList( ftest, NMAX );
}
