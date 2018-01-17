# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(ScatteredInterpolationDemo)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})



add_executable(ScatteredInterpolation 
		main.cpp 
		SimulateNoisyImage.h
		TestBSplineScatteredInterpolation.h
		common.h
		sort.h)

target_link_libraries(ScatteredInterpolation ${ITK_LIBRARIES})
