#include "yafft.h"

#include <iostream>

main()
{
	std::complex<double> X[6] = { 1, 2, 3, 4, 5, 6};
	std::complex<double> Y[6];
	std::complex<double> Z[6];

	yaft::fft2( X, Y, 3, 2);
	yaft::ifft2(Y,Z, 3, 2) ;

	std::cout << " Y \n";
	for( auto & y : Y ) std::cout << y << "  ";
	std::cout << "\n";

	std::cout << " Z \n";
	for( auto & z : Z ) std::cout << z << "  ";
	std::cout << "\n";



}