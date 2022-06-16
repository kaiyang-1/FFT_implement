#!/bin/sh
SRC=fft_conv_1d.cc
TARGET=fft_conv_1d.exe
CXX=g++
CXX_FLAGS="-O2 -Wall -lm -std=c++11 "
CXX_FLAGS="$CXX_FLAGS -g"

rm -rf $TARGET
$CXX $SRC $CXX_FLAGS -o $TARGET
