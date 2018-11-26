cmake_minimum_required(VERSION 3.12)
project(ex6)

set(CMAKE_CXX_STANDARD 14)


if (APPLE)
    set (CMAKE_CXX_FLAGS "-lX11 -L/usr/X11/lib -I/usr/X11/include")
    include_directories(ImageStitching)

    add_executable(ex6
            ImageStitching/CImg.h
            ImageStitching/main.cpp
            ImageStitching/MyBlending.cpp
            ImageStitching/MyBlending.h
            ImageStitching/MyMatching.cpp
            ImageStitching/MyMatching.h
            ImageStitching/MySift.cpp
            ImageStitching/MySift.h
            ImageStitching/stdafx.cpp
            ImageStitching/stdafx.h
            ImageStitching/targetver.h)
endif (APPLE)

if(UNIX AND NOT APPLE)
    include_directories(ImageStitching)

    add_executable(ex6
            ImageStitching/CImg.h
            ImageStitching/main.cpp
            ImageStitching/MyBlending.cpp
            ImageStitching/MyBlending.h
            ImageStitching/MyMatching.cpp
            ImageStitching/MyMatching.h
            ImageStitching/MySift.cpp
            ImageStitching/MySift.h
            ImageStitching/stdafx.cpp
            ImageStitching/stdafx.h
            ImageStitching/targetver.h)
    set (CMAKE_CXX_FLAGS "-lX11 -pthread")
    message(STATUS "Looking for X11...")
    find_package(X11 REQUIRED)
    include_directories(${X11_INCLUDE_DIR})
    target_link_libraries(a4detect-linux ${X11_LIBRARIES})
endif(UNIX AND NOT APPLE)

if (WIN32)
    set (CMAKE_CXX_FLAGS "-lgdi32")
    set (CMAKE_SHARED_LINKER_FLAGS "-lgdi32")
    include_directories(ImageStitching)

    add_executable(ex6
            ImageStitching/CImg.h
            ImageStitching/main.cpp
            ImageStitching/MyBlending.cpp
            ImageStitching/MyBlending.h
            ImageStitching/MyMatching.cpp
            ImageStitching/MyMatching.h
            ImageStitching/MySift.cpp
            ImageStitching/MySift.h
            ImageStitching/stdafx.cpp
            ImageStitching/stdafx.h
            ImageStitching/targetver.h)
endif(WIN32)

include_directories(.)