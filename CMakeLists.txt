cmake_minimum_required(VERSION 3.0)

project("FluidRigidCoupling2D")

file(GLOB LIB_HEADERS "pcgsolver/*.h"
                      "array1.h"
                      "array2.h"
                      "array2_utils.h"
                      "box2dgeometry.h"
                      "coordframe2d.h"
                      "fluidsim.h"
                      "levelset.h"
                      "mat.h"
                      "rigidbody.h"
                      "rigidgeometry.h"
                      "util.h"
                      "vec.h"
                      "vector_math.h")
                      
file(GLOB LIB_SOURCES "pcgsolver/*.cpp"
                      "box2dgeometry.cpp"
                      "fluidsim.cpp"
                      "levelset.cpp"
                      "rigidbody.cpp")

file(GLOB EXEC_HEADERS "gluvi.h" "openglutils.h")
file(GLOB EXEC_SOURCES "gluvi.cpp" "openglutils.cpp" "main.cpp")

add_library(fluidrigidcoupling2d_lib STATIC ${LIB_HEADERS} ${LIB_SOURCES})
target_include_directories(fluidrigidcoupling2d_lib PUBLIC ./)

option(ENABLE_TEST_APP "Build the opengl test app" OFF)

if (ENABLE_TEST_APP)
  find_package(GLUT REQUIRED)
  find_package(OpenGL REQUIRED)

  add_executable(fluidrigidcoupling2d ${EXEC_HEADERS} ${EXEC_SOURCES})
  target_link_libraries(fluidrigidcoupling2d fluidrigidcoupling2d_lib GLUT::GLUT OpenGL::GL OpenGL::GLU)
endif ()
