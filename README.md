# HETEE

## Table of Contents
- [HETEE](#hetee)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
    - [Dependencies](#dependencies)
    - [PPML Examples](#ppml-examples)
  - [Installation](#installation)
    - [Install Microsoft Open Enclave (OE) SDK](#install-microsoft-open-enclave-oe-sdk)
    - [Install HElib](#install-helib)
  - [Run HETEE](#run-hetee)

## Introduction

**An Efficient Hybrid Framework for Privacy-preserving Machine Learning Using HE and TEE**

HETEE is a C++-based framework for privacy-preserving machine learning (PPML) based on Homomorphic Encryption (HE) and Intel SGX.
To accelerate the HE-based computations, HETEE selectively outsources HE-unfriendly computations to the SGX enclave while preseving the integrity and privacy of the computation.

### Dependencies

HETEE runs with three dependencies: 
- [HElib](https://github.com/homenc/HElib): an open-source software library that implements two HE schemes, focusing mostly on effective use of the [Smart-Vercauteren ciphertext packing](https://eprint.iacr.org/2011/133) techniques and the [Gentry-Halevi-Smart optimizations](https://eprint.iacr.org/2012/099): 
  - The implementations of the [Brakerski-Gentry-Vaikuntanathan (BGV)](https://eprint.iacr.org/2011/277) scheme with bootstrapping.
  - The Approximate Number scheme of [Cheon-Kim-Kim-Song (CKKS)](https://eprint.iacr.org/2016/421).
- [Open Enclave SDK](https://github.com/openenclave/openenclave): a hardware-agnostic open source library for developing applications that utilize Hardware-based Trusted Execution Environments (a.k.a, Enclaves). Open Enclave (OE) is an SDK for building enclave applications in C and C++. An enclave application partitions itself into two components: untrusted component (*host*) and trusted component (*enclave*).
- [HEMat](https://github.com/K-miran/HEMat): a software package for performing a secure outsourced matrix computation using HE. HEMat is implemented based on the HE library [HEAAN](https://eprint.iacr.org/2016/421.pdf). It is described in more detail in the [CCS2018](https://dl.acm.org/doi/10.1145/3243734.3243837) paper. **Note**: to support both integer-based and rational number-based HE matrix computations, HETEE re-implements HEMat with HElib based on BGV and CKKS schemes.

### PPML Examples

In this implementation, HETEE provides two examples: linear regression and convolutional neural network (CNN) inference.

**Linear Regression**
- *HETEE* vs HE only baseline

To explore the performance gains of processing integer-based task that HETEE achieves, we compare HETEE with HE only baseline. HE only baseline is an implementation of linear regression with the protection of only HE. We re-implement Wu et al's work [SI-HE](https://github.com/dwu4/fhe-si) with the latest version of HElib and set it as HE only baseline.

**CNN Inference**
- *HETEE* vs HE only baseline (E2DM)
- *HETEE* vs Oblivious baseline
- *HETEE* vs HCNN

In terms of CNN inference, we perform rational number-based computations using CKKS and compare HETEE with three baselines. HE only baseline is adopted from the [CCS2018](https://dl.acm.org/doi/10.1145/3243734.3243837) (called E2DM). We implement E2DM with HElib according to the algorithms designed in the paper. Oblivious baseline is implemented with the usage of [Oblivious Primitives](https://github.com/mc2-project/secure-xgboost) inside the enclave. [HCNN](https://ieeexplore.ieee.org/document/9546527) is the most recent and similar work that utilizing HE and SGX to protect the evaluation of CNN. 

## Installation

The following instructions will create an environment for HETEE. Note that HETEE has only been tested on **Ubuntu 18.04**, so we recommend that you install everything on Ubuntu 18.04. 

### Install Microsoft Open Enclave (OE) SDK

OE provides several options such as Ubuntu 18.04 or 20.04 with SGX hardware or simulation mode. You can check [which SGX level your machine support](https://github.com/openenclave/openenclave/blob/master/docs/GettingStartedDocs/Contributors/building_oe_sdk.md#1-determine-the-sgx-support-level-on-your-developmenttarget-system)  before installing OE. 
After confirming SGX support levels, install OE according to the corresponding [instructions](https://github.com/openenclave/openenclave/blob/master/docs/GettingStartedDocs/install_oe_sdk-Ubuntu_18.04.md). 

Note that HETEE is tested on Ubuntu 18.04 with SGX1+FLC mode. In SGX1+FLC mode, the Open Enclave SDK takes advantage of the Flexible Launch Control mode for better managing architectural enclaves.

### Install HElib

HETEE performs HE-friendly computations (e.g., matrix/vector multiplications) in the host while performing HE-unfriednly computations (e.g., calculate non-linear functions or refresh the HE ciphertexts) inside the enclave. Therefore, we need to build the required libraries [NTL](https://github.com/libntl/ntl) and [GMP](https://gmplib.org/) against GLIBC (host) and MUSL (enclave) C library, respectively. 

**Build the NTL and GMP against GLIBC (i.e., in the host)**

Take the GMP as an example:
```
# Download GMP
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.lz;
lzip -d gmp-6.1.2.tar.lz;
tar -xvf gmp-6.1.2.tar;

# Build GMP
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
```

**Build the NTL and GMP against MUSL (i.e., in the enclave)**

To build GMP, NTL and HElib against MUSL, we resort to build them in MUSL-based Alpine Linux using Microsoft Package Manager [apkman](https://github.com/anakrish/apkman). `apkman` is a self-contained bash script which helps us install the packages in Alpine Linux.

Install `apkman`:
```
wget https://raw.githubusercontent.com/openenclave/openenclave/feature/apkman/tools/apkman/apkman
chmod +x apkman
```
Build GMP using `apkman`:
```
mkdir enclave_gmp_install;

(
  cd gmp-6.1.2 || exit 1;
  apkman exec make clean;
  apkman exec ./configure CC=clang CFLAGS="-g -Og" --prefix=$SCRIPT_DIR/enclave_gmp_install;
  apkman exec make -j16;
  apkman exec make install;
)

cp $SCRIPT_DIR/enclave_gmp_install/lib/libgmp.a $SCRIPT_DIR/../enclave/;
cp $SCRIPT_DIR/enclave_gmp_install/include/gmp.h $SCRIPT_DIR/../enclave/;
```

We provide the script, `build-lib.sh`, for installing all the libraries in the host and enclave. 

**Configure `apkman`**

To configure the `apkman`, we add the following commands in the `CMakeLists.txt` of each directory:
```
# Create an imported executable for use by tests.
add_executable(apkman IMPORTED)
set_target_properties(
  apkman PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/../../../../alpineapkman/apkman)

execute_process(COMMAND ${PROJECT_SOURCE_DIR}/../../../../alpineapkman/apkman help
                        COMMAND_ERROR_IS_FATAL ANY)

# Fetch apkman root folder for accessing includes and libraries.
execute_process(
  COMMAND ${PROJECT_SOURCE_DIR}/../../../../alpineapkman/apkman root COMMAND_ERROR_IS_FATAL
          ANY OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE APKMAN_ROOT)

# Ensure that apkman is initialized.
add_custom_target(
  apkman-init # Setup so that packages can be built from source.
  COMMAND apkman add build-base clang gdb cmake git
  COMMAND apkman --wrap-ld)
```

Note that you may have to change the path of `apkman` in the configuration according to your own local installation path of `apkman`.

## Run HETEE

Take the linear regression `HETEE_LR` as an exmpale: 
```
# Source the openenclaverc file
. /opt/openenclave/share/openenclave/openenclaverc
cd HETEE_LR
mkdir build && cd build
cmake ..
make
make run
```