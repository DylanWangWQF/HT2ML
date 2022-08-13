add_custom_command(
  OUTPUT fosdsc_u.h fosdsc_u.c fosdsc_args.h
  DEPENDS ${CMAKE_SOURCE_DIR}/fosdsc.edl
  COMMAND
    openenclave::oeedger8r --untrusted ${CMAKE_SOURCE_DIR}/fosdsc.edl
    --search-path ${OE_INCLUDEDIR} --search-path
    ${OE_INCLUDEDIR}/openenclave/edl/sgx)

set(CMAKE_CXX_FLAGS "-g -O3 -pthread")

add_executable(fosdsc_host
               host.cpp LoadData.cpp hematrix.cpp ${CMAKE_CURRENT_BINARY_DIR}/fosdsc_u.c)

if (WIN32)
  copy_oedebugrt_target(fosdsc_host_oedebugrt)
  add_dependencies(fosdsc_host fosdsc_host_oedebugrt)
endif ()

# Needed for the generated file switchless_sample_u.h 
target_include_directories(
    fosdsc_host PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(
  fosdsc_host 
  openenclave::oehost
  -L${CMAKE_CURRENT_SOURCE_DIR}/lib
  libhelib.a
  libntl.a
  libgmp.a)