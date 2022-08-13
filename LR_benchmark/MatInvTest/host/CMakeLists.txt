add_custom_command(
  OUTPUT hteva_u.h hteva_u.c hteva_args.h
  DEPENDS ${CMAKE_SOURCE_DIR}/hteva.edl
  COMMAND
    openenclave::oeedger8r --untrusted ${CMAKE_SOURCE_DIR}/hteva.edl
    --search-path ${OE_INCLUDEDIR} --search-path
    ${OE_INCLUDEDIR}/openenclave/edl/sgx)

set(CMAKE_CXX_FLAGS "-g -O3 -pthread")

add_executable(hteva_host
               host.cpp ${CMAKE_CURRENT_BINARY_DIR}/hteva_u.c)

if (WIN32)
  copy_oedebugrt_target(hteva_host_oedebugrt)
  add_dependencies(hteva_host hteva_host_oedebugrt)
endif ()

# Needed for the generated file switchless_sample_u.h 
target_include_directories(
    hteva_host PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(
  hteva_host 
  openenclave::oehost)