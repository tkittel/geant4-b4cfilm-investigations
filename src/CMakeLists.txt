################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

cmake_minimum_required(VERSION 3.25...3.31)

project( G4Film VERSION 0.0.1 LANGUAGES CXX )

find_package( Geant4 11.3.2 REQUIRED )

if( NOT DEFINED NCrystalGeant4_DIR )
  execute_process(
    COMMAND ncrystal-geant4 --cmakedir
    OUTPUT_VARIABLE NCrystalGeant4_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()
find_package( NCrystalGeant4 0.1.0 REQUIRED )

add_executable( g4filmlib "main.cc" )
target_compile_features( g4filmlib PRIVATE cxx_std_17 )

target_link_libraries(
  g4filmlib
  NCrystalGeant4::NCrystalGeant4
  ${Geant4_LIBRARIES}
)

install(TARGETS "g4filmlib" DESTINATION bin )
