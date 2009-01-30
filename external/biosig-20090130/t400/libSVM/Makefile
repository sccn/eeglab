# This Makefile is used under Linux

MATLABDIR ?= /usr/local/matlab
CXX ?= g++
CC ?= gcc
#CXX = g++-3.3
#CC = gcc-3.3
CFLAGS = -Wall -O3 -fPIC -I$(MATLABDIR)/extern/include

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CC) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)"
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)

all: svmpredict.$(MEX_EXT) svmtrain.$(MEX_EXT)

svmpredict.$(MEX_EXT):     svmpredict.c svm.h svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmpredict.c svm.o svm_model_matlab.o

svmtrain.$(MEX_EXT):       svmtrain.c svm.h svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmtrain.c svm.o svm_model_matlab.o

svm_model_matlab.o:     svm_model_matlab.c svm.h
	$(CC) $(CFLAGS) -c svm_model_matlab.c

svm.o:  svm.cpp svm.h
	$(CXX) $(CFLAGS) -c svm.cpp

clean:
	rm -f *~ *.o; rm -f *.mex*; rm -f *obj

