@echo off

set yaml-cpp_ROOT=C:\Users\9blo\miniconda3\envs\rkt\Library\share\cmake\yaml-cpp
 
cmake -G "Visual Studio 16 2019" -A x64 -S . -B build