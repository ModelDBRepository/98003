# Make file for generating rallpack.exe

# This is a straightforward make file that can
# be customized to meet system requirements.
# Linux is used as an example system here.

# A Program Development Environment (PDE) will greatly
# simplify using BNSF and will in many cases create its
# own version of a makefile instead of using this file.

# Preliminary defines for compile options.
CPP = g++
CPPFLAGS = -O -I../Include -I../Include/Rallpack
LD = g++
LDFLAGS = 

# Define a list of all the object modules generated
OBJS= \
	bnsf_base.o \
	bnsf_liaf.o \
	bnsf_math.o \
	bnsf_math_3rd_party.o \
	bnsf_nmod.o \
	bnsf_sim.o \
	ionchan_k_rallpack.o \
	ionchan_na_rallpack.o \
	neuron_rallpack.o \
	test_rallpack.o \
	test_rallpack1.o \
	test_rallpack2.o \
	test_rallpack3.o \
	test_rallpack_ab.o

# Housekeeping dependencies
all : rallpack.exe

clean :
	rm -f $(OBJS) rallpack.exe


# Link the main module
rallpack.exe : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o rallpack.exe 


# Compile source into object modules
bnsf_base.o : \
	../Src/bnsf_base.cpp \
	../Include/bnsf_base.h
	$(CPP) $(CPPFLAGS) -c ../Src/bnsf_base.cpp

bnsf_liaf.o : \
	../Src/bnsf_liaf.cpp \
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/bnsf_liaf.cpp

bnsf_math.o : \
	../Src/bnsf_math.cpp\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h
	$(CPP) $(CPPFLAGS) -c ../Src/bnsf_math.cpp

bnsf_math_3rd_party.o : \
	../Src/bnsf_math_3rd_party.cpp\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h
	$(CPP) $(CPPFLAGS) -c ../Src/bnsf_math_3rd_party.cpp

bnsf_nmod.o : \
	../Src/bnsf_nmod.cpp\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/bnsf_nmod.cpp

bnsf_sim.o : \
	../Src/bnsf_sim.cpp\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/bnsf_sim.cpp

ionchan_k_rallpack.o : \
	../Src/Rallpack/ionchan_k_rallpack.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_k_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Src/Rallpack/ionchan_k_rallpack.cpp

ionchan_na_rallpack.o : \
	../Src/Rallpack/ionchan_na_rallpack.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_na_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Src/Rallpack/ionchan_na_rallpack.cpp

neuron_rallpack.o : \
	../Src/Rallpack/neuron_rallpack.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_k_rallpack.h\
	../Include/Rallpack/ionchan_na_rallpack.h\
	../Include/Rallpack/neuron_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Src/Rallpack/neuron_rallpack.cpp

test_rallpack1.o : \
	../Testcases/test_rallpack1.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_k_rallpack.h\
	../Include/Rallpack/ionchan_na_rallpack.h\
	../Include/Rallpack/neuron_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_rallpack1.cpp

test_rallpack2.o : \
	../Testcases/test_rallpack2.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_k_rallpack.h\
	../Include/Rallpack/ionchan_na_rallpack.h\
	../Include/Rallpack/neuron_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_rallpack2.cpp

test_rallpack3.o : \
	../Testcases/test_rallpack3.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_k_rallpack.h\
	../Include/Rallpack/ionchan_na_rallpack.h\
	../Include/Rallpack/neuron_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_rallpack3.cpp	

test_rallpack_ab.o : \
	../Testcases/test_rallpack_ab.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h\
	../Include/Rallpack/ionchan_k_rallpack.h\
	../Include/Rallpack/ionchan_na_rallpack.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_rallpack_ab.cpp

test_rallpack.o : ../Testcases/test_rallpack.cpp
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_rallpack.cpp

# end makefile for rallpack.exe
