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

# Replace with appropriate paths to
# - Open MPI (include, lib)
# - Eigen (include only)
# - If on Linux, your Intel OneAPI installation (include, lib)
# - MUMPS (include, lib)
ifeq ($(UNAME), Linux)
CC := g++
CFLAGS := -Ofast -fopenmp -std=c++14 -I/usr/mpi/gcc/openmpi-4.1.5a1/include -I/home/USERNAME/.local/Eigen -I/opt/intel/oneapi/2025.0/include -I/home/USERNAME/MUMPS_5.7.3/include
LDFLAGS := -fopenmp -L/usr/mpi/gcc/openmpi-4.1.5a1/lib -L/usr/local/lib -L/opt/intel/oneapi/2025.0/lib -L/home/USERNAME/MUMPS_5.7.3/lib -lmpi -lmpi_cxx -lcmumps -lpord -liomp5 -lmumps_common -lmpi_mpifh
endif
ifeq ($(UNAME), Darwin)
CC=g++-12
CFLAGS=-Ofast -fopenmp -std=c++14 -I ../Eigen/ -I /opt/homebrew/Cellar/open-mpi/5.0.7/include -I /Users/USERNAME/MUMPS_5.7.3/include
LDFLAGS= -fopenmp -L /opt/homebrew/Cellar/open-mpi/5.0.7/lib -L /Users/USERNAME/MUMPS_5.7.3/lib -lmpi -lcmumps 
endif

SOURCES := main.cpp
OBJECTS := $(SOURCES:.cpp=.o)
EXECUTABLE := acmefdfd

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
