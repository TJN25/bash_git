#!/bin/bash
brew install gcc@6
export CC=gcc-6
export CXX=g++-6

git clone https://github.com/soedinglab/WIsH.git
cd WIsH
cmake .
make