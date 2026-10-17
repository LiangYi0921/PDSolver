cmake_minimum_required(VERSION 3.20)
project(PDSolverOpenSource LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find Boost (header-only Spirit Qi for parsing MST logs)
find_package(Boost REQUIRED)

add_executable(PDSolverDemo
    main.cpp
    include/Vector3.h
    include/Matrix3x3.h
    include/Particle.h
    include/Element.h
    include/TrussElement.h
    include/PDSolver.h
    include/MstModel.h
    include/MstParser.h
)

target_include_directories(PDSolverDemo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(PDSolverDemo PRIVATE
    Eigen3::Eigen
)

if(MSVC)
    target_compile_options(PDSolverDemo PRIVATE /utf-8 /EHsc /DNOMINMAX /O2)
else()
    target_compile_options(PDSolverDemo PRIVATE -O3 -Wall -Wextra)
endif()
