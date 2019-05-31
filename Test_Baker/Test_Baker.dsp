# Microsoft Developer Studio Project File - Name="Test_Baker" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Test_Baker - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Test_Baker.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Test_Baker.mak" CFG="Test_Baker - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Test_Baker - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Test_Baker - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Test_Baker - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /G6 /W3 /GX /Zi /O2 /I "..\Include" /I "..\Include\Baker" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /map /debug /machine:I386
# SUBTRACT LINK32 /profile

!ELSEIF  "$(CFG)" == "Test_Baker - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /G6 /W3 /GX /ZI /Od /I "..\Include" /I "..\Include\Baker" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# SUBTRACT CPP /X
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /map /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Test_Baker - Win32 Release"
# Name "Test_Baker - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\Src\bnsf_base.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\bnsf_liaf.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\bnsf_math.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\bnsf_math_3rd_party.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\bnsf_nmod.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\bnsf_sim.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\cell_l56a_25_micron.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\cell_l56a_50_micron.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\cell_l56a_5_micron.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\hippoform_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_ca_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_ih_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_k_a_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_k_ahp_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_k_c_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_k_dr_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_k_m_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\ionchan_na_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\maze_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\mouse_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\neuron_baker_2003.cpp

!IF  "$(CFG)" == "Test_Baker - Win32 Release"

!ELSEIF  "$(CFG)" == "Test_Baker - Win32 Debug"

# ADD CPP /W3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\Src\Baker\placecell_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\plasticity_gaba_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\plasticity_glu_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\subject_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\synapse_gaba_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Src\Baker\synapse_glu_baker_2003.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_010.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_020.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_100.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_110.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_210.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_300.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_310.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_320.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_350.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_600.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_610.cpp
# End Source File
# Begin Source File

SOURCE=..\Testcases\test_baker_670.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\Include\bnsf.h
# End Source File
# Begin Source File

SOURCE=..\Include\bnsf_base.h
# End Source File
# Begin Source File

SOURCE=..\Include\bnsf_liaf.h
# End Source File
# Begin Source File

SOURCE=..\Include\bnsf_math.h
# End Source File
# Begin Source File

SOURCE=..\Include\bnsf_nmod.h
# End Source File
# Begin Source File

SOURCE=..\Include\bnsf_sim.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\hippoform_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_ca_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_ih_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_k_a_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_k_ahp_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_k_c_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_k_dr_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_k_m_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\ionchan_na_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\maze_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\mouse_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\neuron_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\placecell_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\plasticity_gaba_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\plasticity_glu_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\subject_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\synapse_gaba_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\synapse_glu_baker_2003.h
# End Source File
# Begin Source File

SOURCE=..\Include\Baker\test_baker_300.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
