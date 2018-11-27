cmake_minimum_required(VERSION 3.12)
project(ex6)

set(CMAKE_CXX_STANDARD 14)


if (APPLE)
    set (CMAKE_CXX_FLAGS "-lX11 -L/usr/X11/lib -I/usr/X11/include")
    include_directories(ImageStitching)

    add_executable(ex6
            ImageStitching/CImg.h
            ImageStitching/main.cpp
            ImageStitching/Projection.cpp
            ImageStitching/Projection.h
            ImageStitching/vl/sift.h
            ImageStitching/vl/sift.c
            ImageStitching/vl/generic.h
            ImageStitching/vl/generic.c
            ImageStitching/vl/host.h
            ImageStitching/vl/host.c
            ImageStitching/vl/random.h
            ImageStitching/vl/random.c
            ImageStitching/vl/imopv.h
            ImageStitching/vl/imopv.c
            ImageStitching/ImageProcess.cpp
            ImageStitching/ImageProcess.h
            )
endif (APPLE)

if(UNIX AND NOT APPLE)
    include_directories(ImageStitching)

    add_executable(ex6
            ImageStitching/CImg.h
            ImageStitching/main.cpp
            ImageStitching/Projection.cpp
            ImageStitching/Projection.h
            ImageStitching/vl/sift.h
            ImageStitching/vl/sift.c
            ImageStitching/vl/generic.h
            ImageStitching/vl/generic.c
            ImageStitching/vl/host.h
            ImageStitching/vl/host.c
            ImageStitching/vl/random.h
            ImageStitching/vl/random.c
            ImageStitching/vl/imopv.h
            ImageStitching/vl/imopv.c
            ImageStitching/ImageProcess.cpp
            ImageStitching/ImageProcess.h
            )
    set (CMAKE_CXX_FLAGS "-lX11 -pthread")
    message(STATUS "Looking for X11...")
    find_package(X11 REQUIRED)
    include_directories(${X11_INCLUDE_DIR})
    target_link_libraries(a4detect-linux ${X11_LIBRARIES})
endif(UNIX AND NOT APPLE)

if (WIN32)
    set (CMAKE_CXX_FLAGS "-lgdi32 -static-libstdc++ -static-libgcc")
    set (CMAKE_SHARED_LINKER_FLAGS "-lgdi32 -static-libstdc++ -static-libgcc")

    include_directories(ImageStitching)

    add_executable(ex6
            ImageStitching/CImg.h
            ImageStitching/main.cpp
            ImageStitching/Projection.cpp
            ImageStitching/Projection.h
            ImageStitching/vl/sift.h
            ImageStitching/vl/sift.c
            ImageStitching/vl/generic.h
            ImageStitching/vl/generic.c
            ImageStitching/vl/host.h
            ImageStitching/vl/host.c
            ImageStitching/vl/random.h
            ImageStitching/vl/random.c
            ImageStitching/vl/imopv.h
            ImageStitching/vl/imopv.c
            ImageStitching/ImageProcess.cpp
            ImageStitching/ImageProcess.h
            )
endif(WIN32)
#
#ImageStitching/MyBlending.cpp
#ImageStitching/MyBlending.h
#ImageStitching/MyMatching.cpp
#ImageStitching/MyMatching.h
#ImageStitching/MySift.cpp
#ImageStitching/MySift.h
include_directories(.)