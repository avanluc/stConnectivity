set(PROJECT "STCONN")
project(${PROJECT} CXX)

cmake_minimum_required(VERSION 3.0)
FIND_PACKAGE(CUDA REQUIRED)

message("\nTest on cmake v" ${CMAKE_VERSION} "  System: " ${CMAKE_SYSTEM} "\n")

#cmake -D ARCH=35 -D CMAKE_BUILD_TYPE=Release ..
# ------------------------------------------------------------------
# CHECK C++ 11 SUPPORT AND CUDA >= 6.0

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

if(CUDA_VERSION_MAJOR LESS 6)
	message(FATAL_ERROR "\nRequired Cuda Toolkit >= v6.0   Founded: Cuda Toolkit v${CUDA_VERSION_MAJOR}")
endif()

# ------------------------------------------------------------------

include_directories(${PROJECT_SOURCE_DIR}/cub/
					${PROJECT_SOURCE_DIR}/include/cuda/
					${PROJECT_SOURCE_DIR}/Util/)

set(HOST_SRC 	${PROJECT_SOURCE_DIR}/src/cuda/graph.cpp
				${PROJECT_SOURCE_DIR}/src/cuda/stConn.cpp)

set(DEVICE_SRC 	${PROJECT_SOURCE_DIR}/src/cuda/ptx.cu)
			
set(UTIL_SRC 	${PROJECT_SOURCE_DIR}/Util/fUtil.cpp
				${PROJECT_SOURCE_DIR}/Util/readGraphHeader.cpp
				${PROJECT_SOURCE_DIR}/Util/readGraphSTD.cpp
				${PROJECT_SOURCE_DIR}/Util/readUtil.cpp
				${PROJECT_SOURCE_DIR}/Util/mersenne-twister.cpp)

# ------------------------------------------------------------------
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS "")

set(ARCH "" CACHE STRING "architectures")

if (ARCH STREQUAL "20")
	message("\nGPU architecture: 2.0\n")
	list(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_20,code=sm_20,sm_21" "-DARCH=200")
elseif (ARCH STREQUAL "30")
	message("\nGPU architecture: 3.0\n")
	list(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" "-DARCH=300")
elseif (ARCH STREQUAL "35")
	message("\nGPU architecture: 3.5\n")
	list(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_35,code=sm_35" "-DARCH=350" "-rdc=true")
elseif (ARCH STREQUAL "50")
	message("\nGPU architecture: 5.0\n")
	list(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_50,code=sm_50" "-DARCH=500" "-rdc=true")
elseif (ARCH STREQUAL "52")
	message("\nGPU architecture: 5.2\n")
	list(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_52,code=sm_52" "-DARCH=520" "-rdc=true")
else()
	message(FATAL_ERROR "\nplease select GPU architecture: \"cmake -D ARCH=[ 20 | 30 | 35 | 50 | 52 ]\"")
endif()

# ------------------------------------------------------------------
if (CMAKE_BUILD_TYPE STREQUAL "")
	set (CMAKE_BUILD_TYPE Release)
endif()

message("Selected CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} "\n")
if (    (NOT CMAKE_BUILD_TYPE STREQUAL "Release")
	AND (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
	AND (NOT CMAKE_BUILD_TYPE STREQUAL "DebugLine")
	AND (NOT CMAKE_BUILD_TYPE STREQUAL "Info"))
	message(FATAL_ERROR "\nwrong BUILD_TYPE: \"-D CMAKE_BUILD_TYPE=[ Release | Debug | Info | (void) ]\"")
endif()

#debug options
set(CUDA_NVCC_FLAGS_DEBUG "-g" "-G")	# "--ptxas-options=-v"
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")

#debugline
list(APPEND CUDA_NVCC_FLAGS_DEBUGLINE "-g" "-lineinfo")
set(CMAKE_CXX_FLAGS_DEBUGLINE "-g -Wall")

#info
set(CUDA_NVCC_FLAGS_INFO "--ptxas-options=-v")	# "-maxrregcount=32"
set(CMAKE_CXX_FLAGS_INFO "-g -Wall")

#release options
set(CUDA_NVCC_FLAGS_RELEASE "-maxrregcount=32" "-use_fast_math")		# "-w"
set(CMAKE_CXX_FLAGS_RELEASE "-funroll-loops -O3 -march=native -m64")
# ------------------------------------------------------------------

find_library(
	DEVRT_LIBRARY cudadevrt
	ENV LD_LIBRARY_PATH
	PATHS ${CUDA_TOOLKIT_ROOT_DIR}/lib64 ${CUDA_TOOLKIT_ROOT_DIR}/lib 
)

cuda_add_executable(stConnectivity ${PROJECT_SOURCE_DIR}/src/cuda/stConnectivity.cu ${HOST_SRC} ${DEVICE_SRC} ${UTIL_SRC})

target_link_libraries(stConnectivity ${DEVRT_LIBRARY})

set(EXECUTABLE_OUTPUT_PATH ./)
