cmake_minimum_required(VERSION 3.0)
project(HijingFortranCXX)

# Compiler options
set(CMAKE_CXX_STANDARD 17)  # Set C++ standard to 14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(HEPMC_INSTALL_DIR "dire/to/hepmc3/")

# Find Fortran compiler
enable_language(Fortran)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Add executable
add_executable(Hij hijing_wrapper.cpp hipyset1.36.f hijing1.420.f)

# Set Fortran standard
set_target_properties(Hij PROPERTIES LINKER_LANGUAGE Fortran)

# Specify linker flags
target_link_libraries(Hij PRIVATE gfortran ${ROOT_LIBRARIES} ${HEPMC_INSTALL_DIR}/lib64/libHepMC3.so)

# Set include directories
target_include_directories(Hij PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE ${ROOT_INCLUDE_DIRS} ${HEPMC_INSTALL_DIR}/include)

# Workaround for ROOT with C++14
target_compile_definitions(Hij PRIVATE ROOT_FORCE_STL=1)
