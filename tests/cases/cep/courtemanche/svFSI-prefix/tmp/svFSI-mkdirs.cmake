# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/Code"
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-build"
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix"
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix/tmp"
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix/src/svFSI-stamp"
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix/src"
  "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix/src/svFSI-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix/src/svFSI-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/mnt/c/Users/Luisa María Gavier/svFSIplus/courtemanche/svFSI-prefix/src/svFSI-stamp${cfgdir}") # cfgdir has leading slash
endif()
