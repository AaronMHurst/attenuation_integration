#!/bin/bash

rm -rf build
mkdir -p build
cd build
cmake ../
make -j$nproc

echo "\n Go into the build directory to run the application and follow the prompts:"
echo "\n\tcd build"
echo "\t./attenuation\n"
