add_custom_command(
  OUTPUT teehe_u.h teehe_u.c teehe_args.h
  DEPENDS ${CMAKE_SOURCE_DIR}/teehe.edl
  COMMAND
    openenclave::oeedger8r --untrusted ${CMAKE_SOURCE_DIR}/teehe.edl
    --search-path ${OE_INCLUDEDIR} --search-path
    ${OE_INCLUDEDIR}/openenclave/edl/sgx)

set(CMAKE_CXX_FLAGS "-g -O3 -pthread")

# message("CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(teehe_host
               host.cpp ${CMAKE_CURRENT_SOURCE_DIR}/data/DatasetOperations.cpp ${CMAKE_CURRENT_SOURCE_DIR}/HEMat/CKKSmatrix.cpp ${CMAKE_CURRENT_SOURCE_DIR}/utils/LoadData.cpp ${CMAKE_CURRENT_SOURCE_DIR}/utils/MatOpe.cpp ${CMAKE_CURRENT_BINARY_DIR}/teehe_u.c)

if (WIN32)
  copy_oedebugrt_target(teehe_host_oedebugrt)
  add_dependencies(teehe_host teehe_host_oedebugrt)
endif ()

# Needed for the generated file switchless_sample_u.h 
target_include_directories(teehe_host PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_SOURCE_DIR}/data/mnist/include/mnist ${CMAKE_CURRENT_SOURCE_DIR}/HEMat ${CMAKE_CURRENT_SOURCE_DIR}/utils  ${CMAKE_CURRENT_BINARY_DIR})
# /data ${CMAKE_CURRENT_SOURCE_DIR}/data/mnist/include/mnist ${CMAKE_CURRENT_SOURCE_DIR}/HEMat ${CMAKE_CURRENT_SOURCE_DIR}/utils

target_link_libraries(
  teehe_host 
  openenclave::oehost
  -L${CMAKE_CURRENT_SOURCE_DIR}/lib
  libhelib.a
  libntl.a
  libgmp.a)