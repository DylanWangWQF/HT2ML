#!/bin/sh

#--------------------------------------------------------------------------------
# Download GMP, NTL and HElib

# GMP gmp-6.2.1/ 
# Note: gmp-6.2.1 is installed in apkman
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz;
lzip -d gmp-6.2.1.tar.lz;
tar -xvf gmp-6.2.1.tar;

#gmp-6.1.2 is automatically installed on my Ubuntu18.04 computer, if I use gmp-6.2.1, will fail to install NTL on host due to mismatch (6.2.1/6.1.2)
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.lz;
lzip -d gmp-6.1.2.tar.lz;
tar -xvf gmp-6.1.2.tar;

#NTL ntl/src/
git clone https://github.com/libntl/ntl.git;

#HElib HElib/
git clone https://github.com/homenc/HElib.git;
# helib requires NTL to be built with SHARED=on, or modify Line 231 in CMakeLists.txt: add_library(ntl_external STATIC IMPORTED)
(
  cd HElib || exit 1;
  sed -i 's/add_library(ntl_external SHARED IMPORTED)/add_library(ntl_external STATIC IMPORTED)/g' CMakeLists.txt;
)

SCRIPT_DIR=$(dirname $(readlink -f "$0"))
echo $SCRIPT_DIR

#--------------------------------------------------------------------------------

echo "Host: Build GMP start"

mkdir host_gmp_install;

(
  cd gmp-6.1.2 || exit 1;
  make clean;
  ./configure CC=clang CFLAGS="-g -Og" --prefix=$SCRIPT_DIR/host_gmp_install;
  make -j16;
  make install;
)

cp $SCRIPT_DIR/host_gmp_install/lib/libgmp.a $SCRIPT_DIR/../host/;
cp $SCRIPT_DIR/host_gmp_install/include/gmp.h $SCRIPT_DIR/../host/;

echo "Host: Build GMP end"

# #--------------------------------------------------------------------------------

echo "Enclave: Build GMP start"

mkdir enclave_gmp_install;

(
  cd gmp-6.2.1 || exit 1;
  apkman exec make clean;
  apkman exec ./configure CC=clang CFLAGS="-g -Og" --prefix=$SCRIPT_DIR/enclave_gmp_install;
  apkman exec make -j16;
  apkman exec make install;
)

cp $SCRIPT_DIR/enclave_gmp_install/lib/libgmp.a $SCRIPT_DIR/../enclave/;
cp $SCRIPT_DIR/enclave_gmp_install/include/gmp.h $SCRIPT_DIR/../enclave/;

echo "Enclave: Build GMP end"

#--------------------------------------------------------------------------------

echo "Host: Build NTL start"

mkdir host_ntl_install;

(
  cd ntl/src || exit 1;
  make clean
  ./configure CXX=clang++ CXXFLAGS="-g -Og" NTL_GMP_LIP=on NTL_THREADS=on NTL_THREAD_BOOST=on DEF_PREFIX=$SCRIPT_DIR/host_ntl_install GMP_PREFIX=$SCRIPT_DIR/host_gmp_install;
  make -j16;
  make install;
)

# echo "copy .a file"

cp $SCRIPT_DIR/host_ntl_install/lib/libntl.a $SCRIPT_DIR/../host/;
cp -r $SCRIPT_DIR/host_ntl_install/include/NTL/ $SCRIPT_DIR/../host;


echo "Host: Build NTL end"

#--------------------------------------------------------------------------------

echo "Enclave: Build NTL start"

mkdir enclave_ntl_install;

(
  cd ntl/src || exit 1
  apkman exec make clean;
  apkman exec ./configure CXX=clang++ CXXFLAGS="-g -Og" NTL_GMP_LIP=on SHARED=off NTL_THREADS=off NTL_THREAD_BOOST=off DEF_PREFIX=$SCRIPT_DIR/enclave_ntl_install GMP_PREFIX=$SCRIPT_DIR/enclave_gmp_install;
  apkman exec make -j16
  apkman exec make install
)

# echo "copy .a file"

cp $SCRIPT_DIR/enclave_ntl_install/lib/libntl.a $SCRIPT_DIR/../enclave/;
cp -r $SCRIPT_DIR/enclave_ntl_install/include/NTL/ $SCRIPT_DIR/../enclave/;

echo "Enclave: Build NTL end"

#--------------------------------------------------------------------------------

echo "Host: Build helib start"

mkdir host_helib_build host_helib_install;

(
  cd host_helib_build || exit 1;
  cmake ../HElib/ -DCMAKE_CXX_COMPILER=clang++ -DENABLE_THREADS=ON -DGMP_DIR="$SCRIPT_DIR/host_gmp_install" -DNTL_DIR="$SCRIPT_DIR/host_ntl_install" -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/host_helib_install" -DCMAKE_BUILD_TYPE=Debug;
  make -j16;
  make install;
)

# echo "copy .a file"

# lib/ or lib64/ , it depends on machine, so duplicate the cp command
cp $SCRIPT_DIR/host_helib_install/lib/libhelib.a $SCRIPT_DIR/../host;
cp $SCRIPT_DIR/host_helib_install/lib64/libhelib.a $SCRIPT_DIR/../host;
cp -r $SCRIPT_DIR/host_helib_install/include/helib/ $SCRIPT_DIR/../host;

echo "Host: Build helib end"

#--------------------------------------------------------------------------------

echo "Enclave: Build helib start"

mkdir enclave_helib_build enclave_helib_install;

(
  cd enclave_helib_build || exit 1;
  apkman exec make clean;
  apkman exec cmake ../HElib/ -DCMAKE_CXX_COMPILER=clang++ -DENABLE_THREADS=OFF -DGMP_DIR="$SCRIPT_DIR/enclave_gmp_install" -DNTL_DIR="$SCRIPT_DIR/enclave_ntl_install" -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/enclave_helib_install" -DCMAKE_BUILD_TYPE=Debug;
  apkman exec make -j16;
  apkman exec make install;
)

# echo "copy .a file"

# lib/ or lib64/ , it depends on machine, so duplicate the cp command
cp $SCRIPT_DIR/enclave_helib_install/lib/libhelib.a $SCRIPT_DIR/../enclave;
cp $SCRIPT_DIR/enclave_helib_install/lib64/libhelib.a $SCRIPT_DIR/../enclave;
cp -r $SCRIPT_DIR/enclave_helib_install/include/helib/ $SCRIPT_DIR/../enclave;

echo "Enclave: Build helib end"

#--------------------------------------------------------------------------------

rm gmp-6.2.1.tar gmp-6.1.2.tar;
rm -rf gmp-6.2.1 gmp-6.1.2 ntl HElib;
rm -rf host_gmp_install enclave_gmp_install;
rm -rf host_ntl_install enclave_ntl_install;
rm -rf host_helib_build  enclave_helib_build;
rm -rf host_helib_install enclave_helib_install;

#--------------------------------------------------------------------------------