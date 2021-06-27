YAFT stays for Yet Another Fourier Transform

This is a fourier transform library, which uses the Cooley-Turkey and Rader algorithms.
The library generates the hard-coded algorithms for all possible transofration sizes up to some constant (by default 2048)
It uses the prime_numbers.h from voovrat_utils (see https://github.com/voovrat/voovrat-utils ).

The library is compiled in linux with g++ standard c++14 
In windows you may use cygwin.
If you need something else - feel free to update the make_genfft.sh and genfft.cpp to meet your needs. 

How to compile:
  + download voovrat-utils https://github.com/voovrat/voovrat-utils
  + run make_genfft.sh  ( edit if necessary. Edit genfft.cpp if you need for example more than 2048 functions etc )
  + run ./genfft  ( the make_yaft.sh and cpp files will be created)
  + chmod 777 make_yaft.sh
  + run ./make_yaft.sh ( it make take long time to compile)
  + run ./test_yaft. It performs forward and inverse fourier transformations and checks the difference between the initial and forward/inverse converted vectors

How to use: 
  include yafft.h and call fft( input, output, N, xstep, ystep)
  compile your code with libyaft.a library

