# Makefile
#
# Copyright (c) 2025, 2026, Constantine Sideris (sideris@stanford.edu) and Jui-Hung Sun
# (juihungs@usc.edu)
# 
# This program is free software: you can redistribute it and/or modify it under the terms 
# of the GNU Affero General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
# PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License along with 
# this program. If not, see <https://www.gnu.org/licenses/>. 
# 

UNAME := $(shell uname)
MKLROOT = /opt/intel/oneapi/2025.0/mkl

# Replace with appropriate paths to
# - Eigen (include only)
# - If on Linux, your Intel OneAPI installation (include, lib), or
# - If on MacOS, the Accelerate framework
ifeq ($(UNAME), Linux)
CC := g++
CFLAGS := -Ofast -fopenmp -std=c++17 -I/usr/include/eigen3/Eigen -I/usr/include/eigen3 -c
LDFLAGS := -fopenmp -lopenblas -lpthread -lm -ldl
endif
ifeq ($(UNAME), Darwin)
CC=g++-12
CFLAGS=-Ofast -std=c++17 -fopenmp -framework Accelerate -flax-vector-conversions -DACCELERATE_NEW_LAPACK -I ../Eigen/ -c
LDFLAGS= -fopenmp -framework Accelerate
endif

SOURCES := main.cpp
OBJECTS := $(SOURCES:.cpp=.o)
EXECUTABLE := pngf-opt

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
