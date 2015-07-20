
# Makefile for Torus_py

#------------------------------------------------------
HOST=$(shell hostname)
ifneq (,$(findstring .ast.cam.ac.uk,$(HOST)))
	CCOMPILER = /opt/ioa/software/gcc/4.7.2/bin/g++
	DEBUG = # -ggdb3 #-ffast-math -Ofast
	OPTS = -std=c++0x -Wall -O2 -fopenmp -fPIC $(DEBUG)
	GSLPATH = /opt/ioa/software/gsl/1.16/
	BOOSTINCPATH =/opt/ioa/software/boost/1.55.0/include
	BOOSTLIBPATH =-L/opt/ioa/software/boost/1.55.0/lib
	PYTHONINCPATH=/opt/ioa/software/python/2.7.8/include/python2.7/
	PYTHONLIBPATH=-L/opt/ioa/software/python/2.7.8/lib/
else
	CCOMPILER = g++-4.9
	OPTS = -std=c++0x -Wall -O3 -Ofast -ffast-math -fPIC -fopenmp
	BOOSTINCPATH =.
	BOOSTLIBPATH =
	PYTHONINCPATH=/Users/jls/anaconda/include/python2.7/ -I/Users/jls/anaconda/lib/python2.7/site-packages/numpy/core/include/
	PYTHONLIBPATH=-L/Users/jls/anaconda/lib/
	EXTRACLIB = -lc++
endif


#------------------------------------------------------
# General compiler & options
SRC = src/
OBJ = obj/
#------------------------------------------------------
# Set paths to other code

TORUSPATH = ../Torus_github/Torus/
TORUSINCPATH = -I$(TORUSPATH)src/ -I$(TORUSPATH)src/utils/
TORUSLIBPATH = -L$(TORUSPATH)obj/ -L$(TORUSPATH)WDlib/obj/

#------------------------------------------------------
# Set include and library paths

INC = -I../../ -Iinc/ $(TORUSINCPATH) -I$(PYTHONINCPATH) -I$(BOOSTINCPATH)

LIB = $(TORUSLIBPATH) -L$(GSLPATH)lib/ -lTorus -lOther -lPot -lWD -lebf_cpp -lgsl -lgslcblas -lm $(EXTRACLIB)
BOOSTLIB = -lboost_python -lpython2.7
LIBPY = $(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)

#------------------------------------------------------

default: all

all:  Torus_py.so

forceall:	clean all

#------------------------------------------------------

$(OBJ)%.o : $(SRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

Torus_py.so: $(OBJ)Torus_py.o
	$(CCOMPILER) -fPIC -shared $(OPTS) $(OBJ)Torus_py.o $(INC) $(LIB) $(LIBPY) -o Torus_py.so
