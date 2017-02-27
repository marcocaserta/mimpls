#####################################################
# Author : Marco Caserta
# E-mail : marco dot caserta at uni-hamburg dot de
# Date   : 10.03.11
#
# This makefile is used to compile the 
# sources (c++) that implements the
# algorithm for the multi-item multi-period 
# capacitated lot sizing problem with setups.
#
#####################################################
SRCDIR        = src
BINDIR        = bin
OBJDIR        = obj

# ---------------------------------------------------------------------
# ILOG CPLEX stuff
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

# ---------------------------------------------------------------------                                                    
# CPLEXDIR      = /home/ubuntu/opt/ibm/ILOG/cplex1261/cplex
# CONCERTDIR    = /home/ubuntu/opt/ibm/ILOG/cplex1261/concert
CPLEXDIR      = /home/marco/opt/ibm/cplex1261/cplex
CONCERTDIR    = /home/marco/opt/ibm/cplex1261/concert
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
# ---------------------------------------------------------------------                   
# END ILOG CPLEX stuff                                        

CCOPT     = -m64 -O -fPIC -fexceptions  -DIL_STD  -DOPTIMAL
#-DNDEBUG                                                     
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CCFLAGS   = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
#CTIME     = -pg -fprofile-arcs                               

# set default name for the executable
EXEC ?= benders
EXECCPLEX ?= ls

# set compiler
CC = g++

# set debug options
DEBUG = -ggdb

#set optimization level
OPTLEVEL = -O -DEBIAN_BUILDARCH=pentium

#set flags
FLAGS =  -fomit-frame-pointer -pipe -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings
#-Wconversion

default: $(OBJDIR)/clsp.o $(OBJDIR)/options.o $(OBJDIR)/timer.o $(OBJDIR)/mimpls.o
	$(CC) $(CCFLAGS) $(OBJDIR)/options.o $(OBJDIR)/timer.o $(OBJDIR)/mimpls.o $(OBJDIR)/clsp.o -o $(BINDIR)/$(EXEC) $(CCLNFLAGS) 
$(OBJDIR)/clsp.o: $(SRCDIR)/clsp.cpp $(SRCDIR)/mimpls.cpp  
	$(CC) -c $(CCFLAGS) $(SRCDIR)/clsp.cpp -o $(OBJDIR)/clsp.o
$(OBJDIR)/mimpls.o: $(SRCDIR)/mimpls.cpp  
	$(CC) -c $(CCFLAGS) $(SRCDIR)/mimpls.cpp -o $(OBJDIR)/mimpls.o
$(OBJDIR)/options.o: $(SRCDIR)/options.cpp
	$(CC) -c $(FLAGS) $(SRCDIR)/options.cpp -o $(OBJDIR)/options.o
$(OBJDIR)/timer.o: $(SRCDIR)/timer.cpp
	$(CC) -c $(FLAGS) $(SRCDIR)/timer.cpp -o $(OBJDIR)/timer.o

# clean backup files
clean:
	rm *.*~ 
	rm *.o

# create doxygen documentation using "doxygen.conf" file
# the documentation is put into the directory Doc	
doc: $(SRCDIR)/mimplsCplex.cpp doxygen.conf
	doxygen doxygen.conf


#	
