#!/bin/sh
ardir="${drbdir}/Program/CPP_Library"
srcdir="${drbdir}/Program/CPP_Library/SML"
ariccname="libSMLi++.a"
argccname="libSMLg++.a"
echo icpc -std=c++17 -Wall -Wextra -O3 -xHost -qopenmp -parallel -mkl -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c $srcdir/\*.cpp
icpc -std=c++17 -Wall -Wextra -O3 -xHost -qopenmp -parallel -mkl -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c $srcdir/*.cpp
echo rm $ardir/$ariccname
rm $ardir/$ariccname
ar r $ardir/$ariccname *.o
rm ./*.o
echo g++ -std=c++17 -Wall -Wextra -O3 -Xpreprocessor -fopenmp -c $srcdir/\*.cpp
g++ -std=c++17 -Wall -Wextra -O3 -Xpreprocessor -fopenmp -c $srcdir/*.cpp
echo rm $ardir/$argccname
rm $ardir/$argccname
ar r $ardir/$argccname *.o
rm ./*.o
