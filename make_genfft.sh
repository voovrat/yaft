#!/bin/bash

DEBUG_FLAGS="-O0 -g"
RELEASE_FLAGS="-O3"

#FLAGS=$DEBUG_FLAGS 
FLAGS=$RELEASE_FLAGS

STD="--std=c++17"

UTILS=../utils

SRC=$(  cat << EOL
$UTILS/prime_numbers.cpp
genfft.cpp
EOL
)

INC=$( cat << EOL
-I$UTILS
EOL
)

LIB=$( cat << EOL
EOL
)

g++ -o genfft $SRC $FLAGS $INC $STD $LIB


