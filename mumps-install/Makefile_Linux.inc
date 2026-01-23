# Makefile_Linux.inc
#
# Constantine Sideris (sideris@stanford.edu), Jui-Hung Sun (juihungs@usc.edu)
# 16 March 2025
# 
# This file is derived from the MUMPS package, distributed under the 
# CeCILL-C license. limitation of warranty and liability of the licensors are 
# set forth in Articles 8 and 9 of the CeCILL-C license.
#

#
#  This file is part of MUMPS 5.7.1, released
#  on Thu May  2 10:15:09 UTC 2024
#
################################################################################
#
#   Makefile.inc.generic
#
#   This defines some parameters dependent on your platform; you should
#   look for the approriate file in the directory ./Make.inc/ and copy it
#   into a file called Makefile.inc. For example, from the MUMPS root
#   directory, use
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc"
#   (see the main README file for details)
#
#   If you do not find any suitable Makefile in Makefile.inc, use this file:
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc" and modify it according
#   to the comments given below. If you manage to build MUMPS on a new platform,
#   and think that this could be useful to others, you may want to send us
#   the corresponding Makefile.inc file.
#
################################################################################


########################################################################
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. It is recommended to
# install other orderings. For that, you need to obtain the corresponding package
# and modify the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 6.0 and later) orderings are recommended.
#

#SCOTCHDIR  = ${HOME}/scotch_6.0
#ISCOTCH    = -I$(SCOTCHDIR)/include
#
# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord$(PLAT)

# It is strongly recommended to install METIS for 3D calculations. The 3D calculations can be
# extremely slow and time consuming w`ithout METIS. You should specify the library and include
# directory for METIS in the following.

LMETISDIR = /home/username/.local/metis/lib
IMETIS    = -I/home/username/.local/metis/include

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

LMETIS    = -L$(LMETISDIR) -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
# If you want to use Metis 4.X or an older version, you should use -Dmetis4 instead of -Dmetis
# or in addition with -Dparmetis (if you are using parmetis 3.X or older).
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis

ORDERINGSF  = -Dpord -Dmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################

########################################################################

# DEFINE HERE SOME COMMON COMMANDS, THE COMPILER NAMES, ETC...

# PLAT : use it to add a default suffix to the generated libraries
PLAT    =
# Library extension, + C and Fortran "-o" option
# may be different under Windows
LIBEXT  = .a
LIBEXT_SHARED  = .so
FPIC_OPT = -fPIC
SONAME = -soname
# Adapt/uncomment RPATH_OPT to avoid modifying
# LD_LIBRARY_PATH in case of shared libraries
RPATH_OPT = -Wl,-rpath,/home/username/MUMPS_5.7.3/lib/
OUTC    = -o
OUTF    = -o
# RM : remove files
RM = /bin/rm -f
# CC : C compiler
CC = gcc
# FC : Fortran 90 compiler
FC = gfortran
# FL : Fortran linker
FL = gfortran
# AR : Archive object in a library
#      keep a space at the end if options have to be separated from lib name
AR = ar vr
#RANLIB = ranlib
RANLIB  = echo

# Make this variable point to the path where the Intel MKL library is
# installed. It is set to the default install directory for Intel MKL.
MKLROOT=/opt/intel/oneapi/mkl

# DEFINE HERE YOUR LAPACK LIBRARY
LAPACK = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

# SCALAP should define the SCALAPACK and BLACS libraries.
SCALAP = -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

LIBPAR = $(SCALAP) $(LAPACK)

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq$(PLAT)

# DEFINE HERE YOUR BLAS LIBRARY
# BLAS library is under the folder of $(MKLROOT)/lib/intel64
LIBBLAS = -L/usr/mpi/gcc/openmpi-4.1.5a1/lib -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmpi
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -I /usr/mpi/gcc/openmpi-4.1.5a1/include -O -fPIC -fopenmp -fallow-argument-mismatch # -DBLR_MT -fopenmp -DGEMMT_AVAILABLE
OPTL    = -shared -O -fPIC  -fopenmp
OPTC    = -I /usr/mpi/gcc/openmpi-4.1.5a1/include -O -I. -fPIC -fopenmp
#End Optimized options

INCS = $(INCPAR)
LIBS = $(LIBPAR)

LIBSEQNEEDED =
