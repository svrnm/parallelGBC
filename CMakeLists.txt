# parallelGBC — CMake build (library + test-f4 + install / CMake package)
cmake_minimum_required(VERSION 3.16)
project(parallelGBC VERSION 0.9.1 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(ENABLE_MPI "Build with MPI (Boost.MPI + MPI C++)" OFF)

set(PGBC_COEFF_BITS "16" CACHE STRING "Coefficient storage bits (typically 8, 16, or 32)")
set(PGBC_USE_SSE "1" CACHE STRING "Use SIMDe/SIMD field ops (0 or 1)")
set(PGBC_SORTING "1" CACHE STRING "Reduction sorting mode (0 or 1)")
set(PGBC_POST_REDUCE "0" CACHE STRING "Post-reduce with simplify (0 or 1)")
set(PGBC_PARALLEL_SETUP "1" CACHE STRING "Parallel matrix setup (0 or 1)")

if(ENABLE_MPI)
	set(_PGBC_MPI 1)
else()
	set(_PGBC_MPI 0)
endif()

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
# configure_file @ONLY: substitute cache strings and MPI flag (0/1)
set(PGBC_WITH_MPI ${_PGBC_MPI})
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/Definitions.H.in"
	"${CMAKE_CURRENT_BINARY_DIR}/generated/Definitions.H"
	@ONLY)

set(PGBC_SOURCES
	src/CoeffField.C
	src/F4Algorithm.C
	src/F4DefaultReducer.C
	src/F4Simplify.C
	src/F4SimplifyDB.C
	src/F4Utils.C
	src/BuchbergerVerify.C
	src/Polynomial.C
	src/TMonoid.C
	src/TOrdering.C
	src/Term.C)

add_library(f4 STATIC ${PGBC_SOURCES})
add_library(parallelGBC::f4 ALIAS f4)

target_compile_features(f4 PUBLIC cxx_std_20)

target_include_directories(f4
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/parallelGBC>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/simde)

target_compile_options(f4 PRIVATE
	-Wall -Wextra -Wno-unused-parameter -O2)

target_compile_definitions(f4 PUBLIC
	PGBC_COEFF_BITS=${PGBC_COEFF_BITS}
	PGBC_USE_SSE=${PGBC_USE_SSE}
	PGBC_SORTING=${PGBC_SORTING}
	PGBC_POST_REDUCE=${PGBC_POST_REDUCE}
	PGBC_WITH_MPI=${_PGBC_MPI}
	PGBC_PARALLEL_SETUP=${PGBC_PARALLEL_SETUP})

# Apple Clang does not ship OpenMP; Homebrew's libomp lives outside the default search path.
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	execute_process(
		COMMAND brew --prefix libomp
		OUTPUT_VARIABLE _PGBC_LIBOMP_PREFIX
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET
		RESULT_VARIABLE _PGBC_BREW_LIBOMP_RC)
	if(_PGBC_BREW_LIBOMP_RC EQUAL 0 AND EXISTS "${_PGBC_LIBOMP_PREFIX}/include/omp.h")
		set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${_PGBC_LIBOMP_PREFIX}/include")
		set(OpenMP_CXX_LIB_NAMES "omp")
		set(OpenMP_omp_LIBRARY "${_PGBC_LIBOMP_PREFIX}/lib/libomp.dylib")
	endif()
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(Boost REQUIRED COMPONENTS regex)
find_library(PGBC_TBB_LIBRARY NAMES tbb REQUIRED)
find_path(PGBC_TBB_INCLUDE_DIR NAMES tbb/tbb.h REQUIRED)

add_library(parallelGBC::TBB INTERFACE IMPORTED)
set_target_properties(parallelGBC::TBB PROPERTIES
	INTERFACE_LINK_LIBRARIES "${PGBC_TBB_LIBRARY}"
	INTERFACE_INCLUDE_DIRECTORIES "${PGBC_TBB_INCLUDE_DIR}")

target_link_libraries(f4 PUBLIC
	OpenMP::OpenMP_CXX
	Boost::regex
	parallelGBC::TBB)

if(ENABLE_MPI)
	find_package(MPI COMPONENTS CXX REQUIRED)
	find_package(Boost REQUIRED COMPONENTS mpi serialization)
	target_link_libraries(f4 PUBLIC Boost::mpi Boost::serialization MPI::MPI_CXX)
endif()

add_executable(test-f4 test/test-f4.C)
target_include_directories(test-f4 PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_BINARY_DIR}/generated
	${CMAKE_CURRENT_SOURCE_DIR}/simde)
target_link_libraries(test-f4 PRIVATE f4)

enable_testing()
add_test(
	NAME gb_regression
	COMMAND ${CMAKE_COMMAND} -E env
		"TEST_F4_BIN=$<TARGET_FILE:test-f4>"
		"VERIFY_GB=1"
		"VERIFY_MAX_GB=100"
		"VERIFY_TIMEOUT=900"
		"CORE_LIST=1 4"
		bash "${CMAKE_SOURCE_DIR}/test/RunTests.sh"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

file(GLOB PGBC_PUBLIC_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.H")
list(FILTER PGBC_PUBLIC_HEADERS EXCLUDE REGEX ".*/Definitions\\.H$")

install(TARGETS f4 test-f4
	EXPORT parallelGBCTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/generated/Definitions.H"
	${PGBC_PUBLIC_HEADERS}
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/parallelGBC)

install(EXPORT parallelGBCTargets
	FILE parallelGBCTargets.cmake
	NAMESPACE parallelGBC::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parallelGBC)

if(ENABLE_MPI)
	set(PARALLELOGBC_MPI_CONFIG_BLOCK [[find_dependency(MPI COMPONENTS CXX)
find_dependency(Boost COMPONENTS mpi serialization)]])
else()
	set(PARALLELOGBC_MPI_CONFIG_BLOCK "")
endif()

configure_package_config_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/parallelGBCConfig.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/parallelGBCConfig.cmake"
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parallelGBC)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/parallelGBCConfigVersion.cmake"
	COMPATIBILITY SameMajorVersion)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/parallelGBCConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/parallelGBCConfigVersion.cmake"
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parallelGBC)
