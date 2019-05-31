# Make file for generating utility.exe.

# This is a straightforward make file that can
# be customized to meet system requirements.
# Linux is used as an example system here.

# Preliminary defines for compile options.
CPP = g++
CPPFLAGS = -I../Include
LD = g++
LDFLAGS = 

# Define a list of all the object modules generated
OBJS= \
	swc_reader.o \
	test_utility.o

# Housekeeping dependencies
all : utility.exe

clean :
	rm -f $(OBJS) utility.exe


# Link the main module
utility.exe : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o utility.exe

# Compile source into object modules
test_utility.o : ../Testcases/test_utility.cpp
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_utility.cpp

swc_reader.o : ../Src/Utility/swc_reader.cpp
	$(CPP) $(CPPFLAGS) -c ../Src/Utility/swc_reader.cpp

# end makefile for utility.exe
