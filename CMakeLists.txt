cmake_minimum_required(VERSION 3.10)
project(objectdetection)
find_package(OpenCV REQUIRED PATHS "C:/Users/andre/Downloads/opencv/build")

include_directories("C:/Users/andre/Downloads/opencv/build/")
message(STATUS "OPENCV DIREC: ${OpenCV_INCLUDE_DIRS}")
add_executable(${PROJECT_NAME} src/hw2prob4.cpp)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})


