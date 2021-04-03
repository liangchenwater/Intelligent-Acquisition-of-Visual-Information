cmake_minimum_required(VERSION 3.17)
project(iavi)
add_definitions(-std=c++11)

find_package(OpenCV REQUIRED)

include_directories(/opt/X11/include $ENV{OPENCV_INCLUDE})

link_directories(/opt/X11/lib)

add_executable(iavi main.cpp camera.cpp camera.h group.cpp group.h hit.h image.cpp image.h matrix.cpp matrix.h object3d.h ray.h scene_parser.cpp scene_parser.h sphere.cpp sphere.h vectors.h light.h light.cpp plane.h plane.cpp triangle.h triangle.cpp transform.h transform.cpp glcanvas.h glcanvas.cpp material.h material.cpp raytracer.h raytracer.cpp rayTree.h rayTree.cpp)

TARGET_LINK_LIBRARIES(iavi ${OpenCV_LIBS})
target_link_libraries(iavi libGL.1.dylib)
target_link_libraries(iavi libglapi.0.dylib)
target_link_libraries(iavi libGLESv1_CM.1.dylib)
target_link_libraries(iavi libGLESv2.2.dylib)
target_link_libraries(iavi libGLU.1.dylib)
target_link_libraries(iavi libglut.3.dylib)