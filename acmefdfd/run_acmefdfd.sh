#!/bin/bash

# run_acmefdfd.sh
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
#
#
# This script generates the 5 G matrices for the 5 different frequencies
# of optimization {25GHz, 27.5GHz, 30GHz, 32.5GHz, 35GHz} used for optimization
# with the substrate antenna design example
#
# USAGE
#       ./run_acmefdfd.sh
# No arguments are required
#
# NOTES
#       - An output directory OUTPUT_DIR will be created, in which the 
#         G matrices will be stored

OUTPUT_DIR="../pngf-opt"
if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir ${OUTPUT_DIR}
fi

# run acmefdfd for each frequency
./acmefdfd ${OUTPUT_DIR}/Gmat_sub_01.bin 25.0E9
./acmefdfd ${OUTPUT_DIR}/Gmat_sub_02.bin 27.5E9
./acmefdfd ${OUTPUT_DIR}/Gmat_sub_03.bin 30E9
./acmefdfd ${OUTPUT_DIR}/Gmat_sub_04.bin 32.5E9
./acmefdfd ${OUTPUT_DIR}/Gmat_sub_05.bin 35E9
