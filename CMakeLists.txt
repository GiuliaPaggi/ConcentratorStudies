cmake_minimum_required(VERSION 3.16)
project(statistics VERSION 0.1.0)

# Activate test support,including BUILD_TESTING option
include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# use C++17, without non-standard extension 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Activate warnings
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

# Activate address sanitizer and undefined-behaviour sanitizer in debug mode
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer")

find_package(ROOT 6.20 COMPONENTS CONFIG REQUIRED)
#find_package(Boost REQUIRED)

#add_executable(runMynalisis main.cpp Digi.cpp TriggerPrimitive.cpp Segment.cpp Cluster.cpp)
#target_link_libraries(runMyAnalisis PUBLIC ROOT::Core ROOT::Hist ROOT::Physics ROOT::TreePlayer)

# If testing is enabled...
#   to disable testing, run cmake configuration adding-DBUILD_TESTING=OFF
if (BUILD_TESTING)

  # add runTestAnalisis executable
  add_executable(analisis bin/analysis.cpp src/Analyser.cpp src/Cluster.cpp src/TriggerPrimitive.cpp src/Segment.cpp src/Digi.cpp )
  target_include_directories(analisis PRIVATE ./)
  target_link_libraries(analisis PRIVATE ROOT::Core ROOT::Hist ROOT::Physics ROOT::TreePlayer)
  # add runTestAnalisis to test list
  add_test(NAME analisis COMMAND analisis)

endif()
