# HETEE

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
- [HEMat](https://github.com/K-miran/HEMat): a software package for performing a secure outsourced matrix computation using HE. HEMat is implemented based on the HE library [HEAAN](https://eprint.iacr.org/2016/421.pdf). It is described in more detail in the [CCS2018](https://dl.acm.org/doi/10.1145/3243734.3243837) paper. **Note**: to support both integer-based and rational number-based HE matrix computations, HETEE re-implement HEMat with HElib based on BGV and CKKS schemes.

### PPML Examples

In this implementation, HETEE provides two examples: linear regression and convolutional neural network (CNN) inference.

**Linear Regression**
- ***HETEE*** vs HE only baseline

To explore the performance gains that HETEE achieves, we compare HETEE with HE only baseline. HE only baseline is an implementation of linear regression with the protection of only HE. We re-implement Wu et al's work [SI-HE](https://github.com/dwu4/fhe-si) with the latest version of HE and set it as HE only baseline.

**CNN Inference**
- ***HETEE*** vs HE only baseline (E2DM)
- ***HETEE*** vs Oblivious baseline
- ***HETEE*** vs HCNN

In terms of CNN inference, we compare HETEE with three baselines. HE only baseline is adopted from the [CCS2018](https://dl.acm.org/doi/10.1145/3243734.3243837) (called E2DM). We implement E2DM with HElib according to the algorithms designed in the paper. Oblivious baseline is implemented with the usage of [Oblivious Primitives](https://github.com/mc2-project/secure-xgboost) inside the enclave. [HCNN](https://ieeexplore.ieee.org/document/9546527) is the most recent and similar work that utilizing HE and SGX to protect the evaluation of CNN. 

## Installation