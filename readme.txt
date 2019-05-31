This readme.txt file is for model code associated with the paper:

Baker JL, Olds JL (2007) Theta Phase Precession Emerges from a 
Hybrid Computational Model of a CA3 Place Cell. 
Cognitive Neurodynamics 3: 237-248. DOI 10.1007/s11571-007-9018-9.

ABSTRACT

The origins and functional significance of theta phase precession in the 
hippocampus remain obscure, in part, because of the difficulty of 
reproducing hippocampal place cell firing in experimental settings 
where the biophysical underpinnings can be examined in detail. The 
present study concerns a neurobiologically based computational model 
of the emergence of theta phase precession in which the responses of a 
single model CA3 pyramidal cell are examined in the context of 
stimulation by realistic afferent spike trains including those of place 
cells in entorhinal cortex, dentate gyrus, and other CA3 pyramidal 
cells. Spike-timing dependent plasticity in the model CA3 pyramidal 
cell leads to a spatially correlated associational synaptic drive that 
subsequently creates a spatially asymmetric expansion of the model 
cell's place field. Following an initial training period, theta phase 
precession can be seen in the firing patterns of the model CA3 
pyramidal cell. Through selective manipulations of the model it is 
possible to decompose theta phase precession in CA3 into the separate 
contributing factors of inheritance from upstream afferents in the 
dentate gyrus and entorhinal cortex, the interaction of synaptically 
controlled increasing afferent drive with phasic inhibition, and the 
theta phase difference between dentate gyrus granule cell and CA3 
pyramidal cell activity. In the context of a single CA3 pyramidal cell, 
the model shows that each of these factors plays a role in theta phase 
precession within CA3 and suggests that no one single factor offers a 
complete explanation of the phenomenon. The model also shows parallels
between theta phase encoding and pattern completion within the CA3 
autoassociative network.

COPYRIGHT AND SOFTWARE LICENSE

Except as noted below, usage rights for this software and associated 
materials are provided under the terms of the open source MIT License
contained below. A copy of this license is also found at
http://www.opensource.org/licenses/mit-license.php.

Copyright 2007 John L Baker. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

THIRD-PARTY SOFTWARE

File bsnf_math_3rd_party.cpp contains software provided by third-parties
with open source license(s) and disclaimer(s) as described within the file.  

MT19937 portions are copyright 1997 - 2002, Makoto Matsumoto and 
Takuji Nishimura, all rights reserved.

INSTALLATION

Unpack this software from the .zip distribution file into a suitable
location. The following directory structure should then be found:

BNSF_1_0 - top-level directory for programs and (optionally) data results.
|
|-- readme.txt  - this file
|
|-- /Include    - contains C++ header files used. Files beginning with
|                 bnsf pertain to the basic neural simulation framework
|                 (BNSF), which is a general framework for simulating
|                 neuronal models in C++.
|
|---- /Baker    - C++ header files specific to the current model
|
|---- /Rallpack - C++ header files specific to Rallpack test cases
|
|-- /Src        - contains C++ source files for BNSF and the current
|                 model. Files beginning with bnsf are the source
|                 files for the framework.
|
|---- /Baker    - C++ source files for the current model
|
|---- /Rallpack - C++ source files for Rallpack test cases
|
|---- /Utility  - C++ source files for utility programs used in 
|                 converting neuron morphologies from SWC and other
|                 formats to the tables used in the model.
|
|-- /Testcases  - source code for various test case programs. These
|                 potentially could incorporate multiple models but
|                 at present do not do so. File test_baker.cpp is a
|                 driver module for executing test cases associated
|                 with the current model. Only test cases related to
|                 the current publication are included here.
|
|-- /Test_Baker - for MS VC++ users, this is the project directory.
|                 For UNIX users, this is the directory where make
|                 outputs would be directed. Files of the form *.dsw,
|                 *.dsp, and *.opt were created using MS VC++ 6.0.
|                 makefile is a basic make file for generating
|                 the project executable under UNIX. 
|
|-- /Rallpack   - either MS VC++ project directory or Unix make
|                 output directory as above for the Rallpack tests.
|
|-- /Utility    - either MS VC++ project directory or Unix make
|                 output directory as above for the utility programs.
|
|-- /Matlab     - MATLAB source files used in plots of results.
|                 Results data pathnames will need to be adjusted
|                 to conform with the installation directory names.
|                 Files in /Src/Baker are referenced to obtain
|                 cell morphology information used in plots.
|
|-- /R          - R functions used in circular regression statistics,
|                 a plotting function for rallpack test outputs,
|                 and a subset of the types of plots that can
|                 be generated with the MATLAB code noted above. 
|                 See http://www.r-project.org for more about R.
|
|-- /Data       - Default data directory. File l56a.swc is a copy of
|                 morphology data downloaded from the Duke/Southampton
|                 archive of neuronal morphology (originally from URL
|                 http://neuron.duke.edu/cells/).
|
|---- /Default  - Default working directory during execution when
                  using MS VC++ projects. Other locations can be used,
                  but Matlab and R files will need to be adjusted,
                  which will probably be the case anyway. The results
                  of executing test_rallpack3 and test_rallpack_ab
                  are included as example output. Output from the
                  theta phase precession test cases is much larger,
                  on the order of 300 Mb per test execution, and 
                  thus not suitable for inclusion here.

The framework and model code were originally developed using 
MS VC++ 6.0. Project workspaces can be loaded from the corresponding
.dsw files. Minor changes to the original MS VC++ 6.0 code were needed
to facilitate usage with UNIX and GNU compilers.

Use of some form of Program Development Environment (PDE) will 
considerably simply use and inspection of model code. UNIX-style
makefiles supplied here are elementary and PDE tools may be able 
to create a more complete version of the make files to meet their 
specific needs.

Execution of the model code involves building file test_baker.exe. 
Source code for the corresponding C++ file test_baker.cpp will need to 
be adjusted to execute the test case desired. Each test case module 
includes a number of parameters that may similarly be adjusted to meet 
the specific needs of the simulation being performed.

Rallpack test cases are included here for installation verification.
The Rallpack implementation can serve as example of a model that
may be more familiar to many users and can serve as a basic illustration
of the power of using object-oriented models. Rallpack test cases do
not play a role in the model of theta phase precession.

BNSF framework classes are documented using a variation of the Class 
Responsibility Collaborator (CRC) method as described in: Beck and 
Cunningham (1989) A laboratory for teaching object-oriented thinking, 
Proceeding of OOPSLA 1989: 1-6. Class comments are basically in the 
form of a CRC card without explicitly stating collaborators. See file
bnsf_base.h in the Include directory for further information and 
examples. 

Users of BNSF are assumed to be programmers or readers with access to 
source code. As noted above, use of a PDE will make BNSF much more 
accessible. As frameworks go, BNSF should be considered a still immature
member of the species with considerable room for further development. 

For questions or comments regarding this software, please contact 
John L Baker at jbakerb-at-gmu.edu (replace -at- with @).
