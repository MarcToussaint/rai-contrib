cmake_minimum_required(VERSION 3.4)

project(LISRobotOperation)

set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")

option(USE_PYBIND "build python bindings" ON)
option(USE_BULLET "link to bullet" OFF)
option(USE_PHYSX "link to physx" ON)
option(USE_OPENCV "link to opencv" OFF)
option(USE_OPENVR "link to openvr" ON)
option(USE_OPENCV_LOCAL "link to opencv" ON)
option(USE_REALSENSE "link to realsense2 driver library" OFF)
option(USE_LIBFRANKA "link to libfranka" ON)
option(USE_PCL "link to pcl library" OFF)
option(CMAKE_INSTALL_PREFIX "install directory" "/tmp")
set(PYBIND11_PYTHON_VERSION 3.10 CACHE STRING "python version")

add_compile_options(
  -O3 -g -Wall -Wpedantic -fPIC -I/usr/include/openvr -L/usr/lib/openvr -lopenvr_api
  $<$<COMPILE_LANGUAGE:CXX>:-std=c++14>
  $<$<COMPILE_LANGUAGE:CXX>:-Wno-terminate>
  )

################################################################################

add_definitions(
  -DRAI_ROOT_PATH=\"${CMAKE_SOURCE_DIR}/rai\"
  -DRAI_X11 -DRAI_LAPACK -DRAI_EIGEN -DRAI_PTHREAD
  -DRAI_GL
  -DRAI_GLFW -DRAI_ANN -DRAI_QHULL -DRAI_ASSIMP -DNO_CYCLE_DETECTION -DRAI_PNG -DRAI_PLY
  -DRAI_GJK -DRAI_CCD -DRAI_FCL -DRAI_Lewiner #-DRAI_SWIFT
  -DRAI_PORTAUDIO
  -BOOST_ERROR_CODE_HEADER_ONLY
  -DRAI_ROBOTIQ
  -DRAI_VIVECONTROLLER
  )

################################################################################

include_directories(
  rai/rai
  src/
  /usr/include/eigen3
  $ENV{HOME}/opt/include
  )

link_directories(
  $ENV{HOME}/opt/lib
  $ENV{HOME}/opt/bin/linux.clang/release
  )

file(GLOB SRC_Core rai/rai/Core/*.cpp)
file(GLOB SRC_Optim rai/rai/Optim/*.cpp)
file(GLOB SRC_Geo rai/rai/Geo/*.cpp)
file(GLOB SRC_GeoGJK rai/rai/Geo/GJK/*.c)
file(GLOB SRC_GeoPly rai/rai/Geo/ply/*.c)
file(GLOB SRC_GeoCcd rai/rai/Geo/ccd_rai/*.c)
file(GLOB SRC_GeoLewiner rai/rai/Geo/Lewiner/*.cpp)
file(GLOB SRC_Algo rai/rai/Algo/*.cpp)
file(GLOB SRC_Search rai/rai/Search/*.cpp)
file(GLOB SRC_PathAlgos rai/rai/PathAlgos/*.cpp)
file(GLOB SRC_Gui rai/rai/Gui/*.cpp)
file(GLOB SRC_Kin rai/rai/Kin/*.cpp)
#file(GLOB SRC_KinSwift rai/rai/Kin/SWIFT/*.cpp)
file(GLOB SRC_KOMO rai/rai/KOMO/*.cpp)
file(GLOB SRC_Logic rai/rai/Logic/*.cpp)
file(GLOB SRC_LGP rai/rai/LGP/*.cpp)
file(GLOB SRC_Perception rai/rai/Perception/*.cpp)
file(GLOB SRC_Control rai/rai/Control/*.cpp)
file(GLOB SRC_ry rai/rai/ry/*.cpp)

file(GLOB SRC_BotOp src/BotOp/*.cpp)
file(GLOB SRC_Franka src/Franka/*.cpp)
file(GLOB SRC_OptiTrack src/OptiTrack/*.cpp)
file(GLOB SRC_ViveController src/ViveController/*.cpp)
file(GLOB SRC_RealSense src/RealSense/*.cpp)
file(GLOB SRC_Robotiq src/Robotiq/*.cpp)
file(GLOB SRC_Audio src/Audio/*.cpp)
file(GLOB SRC_MarkerVision src/MarkerVision/*.cpp)

add_library(rai SHARED
  ${SRC_Core}
  ${SRC_Optim}
  ${SRC_Geo}
  ${SRC_GeoGJK}
  ${SRC_GeoPly}
  ${SRC_GeoCcd}
  ${SRC_GeoLewiner}
  ${SRC_Algo}
  ${SRC_Search}
  ${SRC_PathAlgos}
  ${SRC_Gui}
  ${SRC_Kin}
  ${SRC_KOMO}
  ${SRC_Logic}
  #${SRC_LGP}
  ${SRC_Perception}
  ${SRC_Control}
  ${SRC_BotOp}
  ${SRC_Franka}
  ${SRC_OptiTrack}
  ${SRC_ViveController}
  ${SRC_RealSense}
  ${SRC_Robotiq}
  ${SRC_Audio}
  ${SRC_MarkerVision}
  )

target_link_libraries(rai
  jsoncpp lapack blas pthread rt X11
  assimp ccd fcl qhull ann glfw GLEW glut GLU GL png
  portaudio
  boost_system 
  )

################################################################################

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#target_link_libraries (rai Eigen3::Eigen)

if(USE_BULLET)
  find_package(Bullet REQUIRED)
  add_definitions(-DRAI_BULLET)
  include_directories(${BULLET_INCLUDE_DIRS})
  target_link_libraries(rai ${BULLET_LIBRARIES})
  message(STATUS "[rai] using bullet libs:" ${BULLET_LIBRARIES})
endif()

if(USE_PHYSX)
  add_definitions(-DRAI_PHYSX -D_DEBUG -DPX_DISABLE_FLUIDS -DCORELIB -DPX32 -DLINUX)
  link_directories($ENV{HOME}/opt/bin/linux.clang/release)
  target_link_libraries(rai PhysXExtensions_static_64 PhysX_static_64 PhysXCooking_static_64 PhysXCommon_static_64 PhysXFoundation_static_64)
  message(STATUS "[rai] using physx libs from ~/opt:")
endif()

if(USE_OPENCV)
  find_package(OpenCV REQUIRED)
  add_definitions(-DRAI_OPENCV)
  target_link_libraries(rai ${OpenCV_LIBS})
  message(STATUS "[rai] using OpenCV libs:" ${OpenCV_LIBS})
endif()

if(USE_OPENCV_LOCAL)
  add_definitions(-DRAI_OPENCV)
  include_directories($ENV{HOME}/opt/include/opencv4)
  target_link_libraries(rai opencv_core opencv_highgui opencv_imgproc opencv_videoio)
  message(STATUS "[rai] using OpenCV libs from ~/opt:")
endif()

#if(USE_OPENVR)
#  find_package(OpenVR REQUIRED)
#  add_definitions(-DRAI_VIVECONTROLLER)
#  target_link_libraries(rai ${OpenVR_LIBS})
#  message(STATUS "[rai] using OpenVR libs:" ${OpenVR_LIBS})
#endif()

target_link_libraries(rai openvr_api)
link_directories(
  $ENV{HOME}/opt/lib/vrqtlibs
)

if(USE_REALSENSE)
  add_definitions(-DRAI_REALSENSE)
  target_link_libraries(rai realsense2)
  message(STATUS "[rai] using realsense libs")
endif()

if(USE_LIBFRANKA)
  add_definitions(-DRAI_FRANKA)
  target_link_libraries(rai franka)
  message(STATUS "[rai] using libfranka")
endif()

################################################################################

if(USE_PYBIND)
  set(PYTHON_EXECUTABLE "python${PYBIND11_PYTHON_VERSION}")
  message(STATUS "[rai] using python executable " ${PYTHON_EXECUTABLE} )
  execute_process(
    COMMAND bash "-c" "${PYTHON_EXECUTABLE} -m pybind11 --cmakedir" 
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE pybind11_DIR
  )
  message(STATUS "[rai] using pybind11_DIR " ${pybind11_DIR} )
  
  find_package(pybind11)
  add_definitions(-DRAI_PYBIND)
  include_directories(${pybind11_INCLUDE_DIRS})
  
  pybind11_add_module(libry SHARED
    ${SRC_ry}
    src/BotOp/pyBot.cxx
  )

  target_link_libraries(libry PRIVATE
    rai
  )
endif()

################################################################################

install(DIRECTORY rai/rai DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hxx" PATTERN "*.ipp"
)

install(TARGETS rai libry DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ COMPONENT pkgconfig)

################################################################################

add_executable(kinEdit rai/bin/src_kinEdit/main.cpp)
target_link_libraries(kinEdit rai)

add_executable(bot bin/bot/main.cpp)
target_link_libraries(bot rai)


