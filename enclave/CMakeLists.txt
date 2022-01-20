add_custom_command(
  OUTPUT fosdsc_t.h fosdsc_t.c fosdsc_args.h
  DEPENDS ${CMAKE_SOURCE_DIR}/fosdsc.edl
  COMMAND
    openenclave::oeedger8r --trusted ${CMAKE_SOURCE_DIR}/fosdsc.edl
    --search-path ${OE_INCLUDEDIR} --search-path
    ${OE_INCLUDEDIR}/openenclave/edl/sgx)

# add_custom_command(
#   OUTPUT dispatcher.o
#   COMMAND apkman exec clang++ -std=c++17 -fPIE -g -c -I${CMAKE_CURRENT_SOURCE_DIR}/include -I/opt/openenclave/include ${CMAKE_CURRENT_SOURCE_DIR}/dispatcher.cpp
# )

add_custom_command(
  OUTPUT dispatcher.o
  COMMAND apkman exec clang++ -std=c++17 -fPIE -g -c -I${CMAKE_CURRENT_SOURCE_DIR}/include -I/opt/openenclave/include ${CMAKE_CURRENT_SOURCE_DIR}/dispatcher.cpp
)

add_executable(enclave ecalls.cpp dispatcher.o stubs.c ${CMAKE_CURRENT_BINARY_DIR}/fosdsc_t.c)


if (WIN32)
  maybe_build_using_clangw(enclave)
endif ()

target_compile_definitions(enclave PUBLIC OE_API_VERSION=2)

# APKMAN_ROOT /home/dylan/.apkman/alpine-fs
message("APKMAN_ROOT ${APKMAN_ROOT}")

target_include_directories(
  enclave
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
          #lib include path
          ${CMAKE_CURRENT_SOURCE_DIR}
          # C++
          ${APKMAN_ROOT}/usr/include/c++/10.2.1
          ${APKMAN_ROOT}/usr/include/c++/10.2.1/x86_64-alpine-linux-musl
          # C
          ${APKMAN_ROOT}/usr/include)

target_compile_options(
  enclave
  PRIVATE
  -g
  -Wno-shorten-64-to-32
  -Wno-bitwise-op-parentheses
  -Wno-shift-op-parentheses
  -Wno-sign-conversion
  -Wno-unused-parameter
  -Wno-implicit-int-conversion
  -ftls-model=local-exec)

target_link_libraries(
    enclave 
    -nostdlib
    # Add library path and libraries.
    -L${CMAKE_CURRENT_SOURCE_DIR}/lib
    libhelib.a
    libntl.a
    libgmp.a
    openenclave::oecore
    -L${APKMAN_ROOT}/usr/lib
    libstdc++.a
    -L${APKMAN_ROOT}/usr/lib/gcc/x86_64-alpine-linux-musl/10.2.1
    libgcc.a
    openenclave::oeenclave
    openenclave::oecrypto${OE_CRYPTO_LIB}
    openenclave::oelibcxx
    openenclave::oelibc
    # Use apkman's ld to avoid bug with Ubuntu 18.04's binutils.
    -fuse-ld=${CMAKE_CURRENT_SOURCE_DIR}/ld
    # Emit warnings for unresolved symbols.
    -Wl,--warn-unresolved-symbols
    -Wl,--wrap,setlocale
    -Wl,--wrap,oe_SYS_clock_gettime_impl)