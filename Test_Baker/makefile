# Make file for generating test_baker.exe

# This is a straightforward make file that can
# be customized to meet system requirements.
# Linux is used as an example system here.

# A Program Development Environment (PDE) will greatly
# simplify using BNSF and will in many cases create its
# own version of a makefile instead of using this file.

# Preliminary defines for compile options.
CPP = g++
CPPFLAGS = -O -I../Include -I../Include/Baker
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
	cell_l56a_25_micron.o \
	cell_l56a_50_micron.o \
	cell_l56a_5_micron.o \
	hippoform_baker_2003.o \
	ionchan_ca_baker_2003.o \
	ionchan_ih_baker_2003.o \
	ionchan_k_a_baker_2003.o \
	ionchan_k_ahp_baker_2003.o \
	ionchan_k_c_baker_2003.o \
	ionchan_k_dr_baker_2003.o \
	ionchan_k_m_baker_2003.o \
	ionchan_na_baker_2003.o \
	maze_baker_2003.o \
	mouse_baker_2003.o \
	neuron_baker_2003.o \
	placecell_baker_2003.o \
	plasticity_gaba_baker_2003.o \
	plasticity_glu_baker_2003.o \
	subject_baker_2003.o \
	synapse_gaba_baker_2003.o \
	synapse_glu_baker_2003.o \
	test_baker.o \
	test_baker_010.o \
	test_baker_020.o \
	test_baker_100.o \
	test_baker_110.o \
	test_baker_210.o \
	test_baker_300.o \
	test_baker_310.o \
	test_baker_320.o \
	test_baker_350.o \
	test_baker_600.o \
	test_baker_610.o \
	test_baker_670.o

# Housekeeping dependencies
all : test_baker.exe

clean :
	rm -f $(OBJS) test_baker.exe


# Link the main module
test_baker.exe : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o test_baker.exe 


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

cell_l56a_25_micron.o : \
	../Src/Baker/cell_l56a_25_micron.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/cell_l56a_25_micron.cpp

cell_l56a_50_micron.o : \
	../Src/Baker/cell_l56a_50_micron.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/cell_l56a_50_micron.cpp

cell_l56a_5_micron.o : \
	../Src/Baker/cell_l56a_5_micron.cpp\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/cell_l56a_5_micron.cpp

hippoform_baker_2003.o : \
	../Src/Baker/hippoform_baker_2003.cpp\
	../Include/Baker/hippoform_baker_2003.h\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/mouse_baker_2003.h\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/placecell_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/hippoform_baker_2003.cpp

ionchan_ca_baker_2003.o : \
	../Src/Baker/ionchan_ca_baker_2003.cpp\
	../Include/Baker/ionchan_ca_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_ca_baker_2003.cpp

ionchan_ih_baker_2003.o : \
	../Src/Baker/ionchan_ih_baker_2003.cpp\
	../Include/Baker/ionchan_ih_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_ih_baker_2003.cpp

ionchan_k_a_baker_2003.o : \
	../Src/Baker/ionchan_k_a_baker_2003.cpp\
	../Include/Baker/ionchan_k_a_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_k_a_baker_2003.cpp

ionchan_k_ahp_baker_2003.o : \
	../Src/Baker/ionchan_k_ahp_baker_2003.cpp\
	../Include/Baker/ionchan_k_ahp_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_k_ahp_baker_2003.cpp

ionchan_k_c_baker_2003.o : \
	../Src/Baker/ionchan_k_c_baker_2003.cpp\
	../Include/Baker/ionchan_k_c_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_k_c_baker_2003.cpp

ionchan_k_dr_baker_2003.o : \
	../Src/Baker/ionchan_k_dr_baker_2003.cpp\
	../Include/Baker/ionchan_k_dr_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_k_dr_baker_2003.cpp

ionchan_k_m_baker_2003.o : \
	../Src/Baker/ionchan_k_m_baker_2003.cpp\
	../Include/Baker/ionchan_k_m_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_k_m_baker_2003.cpp

ionchan_na_baker_2003.o : \
	../Src/Baker/ionchan_na_baker_2003.cpp\
	../Include/Baker/ionchan_na_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/ionchan_na_baker_2003.cpp

maze_baker_2003.o : \
	../Src/Baker/maze_baker_2003.cpp\
	../Include/Baker/maze_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/maze_baker_2003.cpp

mouse_baker_2003.o : \
	../Src/Baker/mouse_baker_2003.cpp\
	../Include/Baker/hippoform_baker_2003.h\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/mouse_baker_2003.h\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/placecell_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/mouse_baker_2003.cpp

neuron_baker_2003.o : \
	../Src/Baker/neuron_baker_2003.cpp\
	../Include/Baker/ionchan_ca_baker_2003.h\
	../Include/Baker/ionchan_ih_baker_2003.h\
	../Include/Baker/ionchan_k_a_baker_2003.h\
	../Include/Baker/ionchan_k_ahp_baker_2003.h\
	../Include/Baker/ionchan_k_c_baker_2003.h\
	../Include/Baker/ionchan_k_dr_baker_2003.h\
	../Include/Baker/ionchan_k_m_baker_2003.h\
	../Include/Baker/ionchan_na_baker_2003.h\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_gaba_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/neuron_baker_2003.cpp

placecell_baker_2003.o : \
	../Src/Baker/placecell_baker_2003.cpp\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/placecell_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/placecell_baker_2003.cpp

plasticity_gaba_baker_2003.o : \
	../Src/Baker/plasticity_gaba_baker_2003.cpp\
	../Include/Baker/plasticity_gaba_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/plasticity_gaba_baker_2003.cpp

plasticity_glu_baker_2003.o : \
	../Src/Baker/plasticity_glu_baker_2003.cpp\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/plasticity_glu_baker_2003.cpp

subject_baker_2003.o : \
	../Src/Baker/subject_baker_2003.cpp\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/subject_baker_2003.cpp

synapse_gaba_baker_2003.o : \
	../Src/Baker/synapse_gaba_baker_2003.cpp\
	../Include/Baker/plasticity_gaba_baker_2003.h\
	../Include/Baker/synapse_gaba_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/synapse_gaba_baker_2003.cpp

synapse_glu_baker_2003.o : \
	../Src/Baker/synapse_glu_baker_2003.cpp\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Src/Baker/synapse_glu_baker_2003.cpp

test_baker_010.o : \
	../Testcases/test_baker_010.cpp\
	../Include/Baker/ionchan_ih_baker_2003.h\
	../Include/Baker/ionchan_k_a_baker_2003.h\
	../Include/Baker/ionchan_k_dr_baker_2003.h\
	../Include/Baker/ionchan_k_m_baker_2003.h\
	../Include/Baker/ionchan_na_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_010.cpp

test_baker_020.o : \
	../Testcases/test_baker_020.cpp\
	../Include/Baker/ionchan_ca_baker_2003.h\
	../Include/Baker/ionchan_k_ahp_baker_2003.h\
	../Include/Baker/ionchan_k_c_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_020.cpp

test_baker_100.o : \
	../Testcases/test_baker_100.cpp\
	../Include/Baker/neuron_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_100.cpp

test_baker_110.o : \
	../Testcases/test_baker_110.cpp\
	../Include/Baker/neuron_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_110.cpp

test_baker_210.o : \
	../Testcases/test_baker_210.cpp\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_210.cpp
test_baker_300.o : \
	../Testcases/test_baker_300.cpp\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_gaba_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/Baker/test_baker_300.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_300.cpp

test_baker_310.o : \
	../Testcases/test_baker_310.cpp\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_gaba_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/Baker/test_baker_300.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_310.cpp

test_baker_320.o : \
	../Testcases/test_baker_320.cpp\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_gaba_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/Baker/test_baker_300.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_320.cpp

test_baker_350.o : \
	../Testcases/test_baker_350.cpp\
	../Include/Baker/plasticity_glu_baker_2003.h\
	../Include/Baker/synapse_gaba_baker_2003.h\
	../Include/Baker/synapse_glu_baker_2003.h\
	../Include/Baker/test_baker_300.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_350.cpp

test_baker_600.o : \
	../Testcases/test_baker_600.cpp\
	../Include/Baker/hippoform_baker_2003.h\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/mouse_baker_2003.h\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/placecell_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_600.cpp

test_baker_610.o : \
	../Testcases/test_baker_610.cpp\
	../Include/Baker/hippoform_baker_2003.h\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/mouse_baker_2003.h\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/placecell_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_610.cpp

test_baker_670.o : \
	../Testcases/test_baker_670.cpp\
	../Include/Baker/hippoform_baker_2003.h\
	../Include/Baker/maze_baker_2003.h\
	../Include/Baker/mouse_baker_2003.h\
	../Include/Baker/neuron_baker_2003.h\
	../Include/Baker/placecell_baker_2003.h\
	../Include/Baker/subject_baker_2003.h\
	../Include/bnsf.h\
	../Include/bnsf_base.h\
	../Include/bnsf_liaf.h\
	../Include/bnsf_math.h\
	../Include/bnsf_nmod.h\
	../Include/bnsf_sim.h
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker_670.cpp

test_baker.o : ../Testcases/test_baker.cpp
	$(CPP) $(CPPFLAGS) -c ../Testcases/test_baker.cpp

# end makefile for test_baker.exe
