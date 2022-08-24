cmake_minimum_required(VERSION 3.11)

# set (CMAKE_C_COMPILER "/usr/local/bin/gcc")
# set (CMAKE_CXX_COMPILER "/usr/local/bin/g++")

# If the CC environment variable has been specified or if the CMAKE_C_COMPILER
# cmake variable has been passed to cmake, use the C compiler that has been
# specified. Otherwise, prefer clang. Same for C++ compiler.
# This must be done before the `project` command.
if (UNIX)
  if (NOT DEFINED ENV{CC} AND NOT DEFINED CMAKE_C_COMPILER)
    find_program(CMAKE_C_COMPILER clang-10 clang)
    # find_program(CMAKE_C_COMPILER gcc-10 gcc)
  endif ()
  if (NOT DEFINED ENV{CXX} AND NOT DEFINED CMAKE_CXX_COMPILER)
    find_program(CMAKE_CXX_COMPILER clang++-10 clang++)
    # find_program(CMAKE_CXX_COMPILER g++-10 g++)
  endif ()
endif ()

project("HElib Sample" LANGUAGES C CXX)

find_package(OpenEnclave CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(OE_CRYPTO_LIB
        mbedtls
        CACHE STRING "Crypto library used by enclaves.")

option(USE_AVX2 "Use AVX2 instructions to speed up oblivious primitives" ON)

# Create an imported executable for use by tests.
add_executable(apkman IMPORTED)
set_target_properties(
  apkman PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/../../../../../alpineapkman/apkman)

execute_process(COMMAND ${PROJECT_SOURCE_DIR}/../../../../../alpineapkman/apkman help
                        COMMAND_ERROR_IS_FATAL ANY)

# Fetch apkman root folder for accessing includes and libraries.
execute_process(
  COMMAND ${PROJECT_SOURCE_DIR}/../../../../../alpineapkman/apkman root COMMAND_ERROR_IS_FATAL
          ANY OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE APKMAN_ROOT)

# Ensure that apkman is initialized.
add_custom_target(
  apkman-init # Setup so that packages can be built from source.
  COMMAND apkman add build-base clang gdb cmake git
  COMMAND apkman --wrap-ld)

add_subdirectory(enclave)
add_subdirectory(host)

# Generate key
add_custom_command(
  OUTPUT private.pem public.pem
  COMMAND openssl genrsa -out private.pem -3 3072
  COMMAND openssl rsa -in private.pem -pubout -out public.pem)

# Sign enclave
add_custom_command(
  OUTPUT enclave/enclave.signed
  DEPENDS enclave enclave/fosdsc.conf private.pem
  COMMAND openenclave::oesign sign -e $<TARGET_FILE:enclave> -c
          ${CMAKE_SOURCE_DIR}/enclave/fosdsc.conf -k private.pem)

add_custom_target(sign ALL DEPENDS enclave/enclave.signed)

if ((NOT DEFINED ENV{OE_SIMULATION}) OR (NOT $ENV{OE_SIMULATION}))
  add_custom_target(
    run
    DEPENDS fosdsc_host sign
    COMMAND fosdsc_host ${CMAKE_BINARY_DIR}/enclave/enclave.signed)
endif ()

add_custom_target(
  simulate
  DEPENDS fosdsc_host sign
  COMMAND fosdsc_host ${CMAKE_BINARY_DIR}/enclave/enclave.signed --simulate)