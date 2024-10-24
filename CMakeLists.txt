cmake_minimum_required(VERSION 3.22)

project(ipc)

add_library(${PROJECT_NAME} SHARED
    src/BarrierEnergy.cpp
    src/GravityEnergy.cpp
    src/InertiaEnergy.cpp
    src/MassSpringEnergy.cpp
    src/main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC ./include)
