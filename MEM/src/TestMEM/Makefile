#------------------------------------------------------------------------------
#Madgraph definitions

#TTH
TTH_MADLIBDIR=../Madgraph/PROC_SA_CPP_sm_4/lib
TTH_MADINCDIR=../Madgraph/PROC_SA_CPP_sm_4/src
TTH_MADPROCDIR=../Madgraph/PROC_SA_CPP_sm_4/SubProcesses/P0_Sigma_sm_gg_ttxh
TTH_MADMODELLIB=model_sm

#T->bW->bjj
THAD_MADLIBDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_tbwjj/lib
THAD_MADINCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_tbwjj/src
THAD_MADPROCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_tbwjj/SubProcesses/P0_Sigma_sm_t_budx
THAD_MADMODELLIB=model_sm_tbwjj

#Anti T->bW->bjj
ATHAD_MADLIBDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_antitbwjj/lib
ATHAD_MADINCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_antitbwjj/src
ATHAD_MADPROCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_antitbwjj/SubProcesses/P0_Sigma_sm_tx_bxdux
ATHAD_MADMODELLIB=model_sm_antitbwjj

#T->bW->blnu
TLEP_MADLIBDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_tbwlnu/lib
TLEP_MADINCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_tbwlnu/src
TLEP_MADPROCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_tbwlnu/SubProcesses/P0_Sigma_sm_t_bepve
TLEP_MADMODELLIB=model_sm_tbwlnu

#H->WW->2l2nu
H2L2NU_MADLIBDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_hw2l2nu/lib
H2L2NU_MADINCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_hw2l2nu/src
H2L2NU_MADPROCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_hw2l2nu/SubProcesses/P0_Sigma_sm_h_epveemvex
H2L2NU_MADMODELLIB=model_sm_hw2l2nu

#TTGamma*
TTLL_MADLIBDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_ggttll/lib
TTLL_MADINCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_ggttll/src
TTLL_MADPROCDIR=../Madgraph/PROC_SA_CPP_sm_DECAY_ggttll/SubProcesses/P0_Sigma_sm_gg_ttxepem
TTLL_MADMODELLIB=model_sm_ggttll


#------------------------------------------------------------------------------
# GoSam definitions

GOSAM_DIR=

#------------------------------------------------------------------------------
# ROOT definitions
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lMathMore

#------------------------------------------------------------------------------
# LHAPDF definitions
#LHAPDFINCDIR=-ILHAPDF-6.1.5-install/include -I/opt/local/include/
#LHAPDFLIB=-LLHAPDF-6.1.5-install/lib/ -lLHAPDF
#LHAPDFINCDIR=-I/grid_mnt/home/nchanon/LHAPDF-6.1.5-install/include
#LHAPDFLIB=-L/grid_mnt/home/nchanon/LHAPDF-6.1.5-install/lib/ -lLHAPDF 
LHAPDFINCDIR=-I/afs/cern.ch/work/c/chanon/LHAPDF-6.1.5-install/include
LHAPDFLIB=-L/afs/cern.ch/work/c/chanon/LHAPDF-6.1.5-install/lib -lLHAPDF 

#------------------------------------------------------------------------------
# g++ options

CXX = g++
CXXFLAGS = -g -O -fPIC -Wall -Wno-deprecated  
LD = g++

SOFLAGS = #-shared #-dynamiclib #-single_module #-undefined dynamic_lookup

CXXFLAGS     += $(ROOTCFLAGS) -I${TTH_MADPROCDIR} -I$(TTH_MADINCDIR) -I${THAD_MADPROCDIR} -I$(THAD_MADINCDIR) -I${ATHAD_MADPROCDIR} -I$(ATHAD_MADINCDIR) -I${TLEP_MADPROCDIR} -I$(TLEP_MADINCDIR) -I${H2L2NU_MADPROCDIR} -I$(H2L2NU_MADINCDIR) -I${TTLL_MADPROCDIR} -I$(TTLL_MADINCDIR) ${LHAPDFINCDIR}

LIBS         = ${ROOTGLIBS} ${ROOTLIBS} -L$(TTH_MADLIBDIR) -l$(TTH_MADMODELLIB) -L$(THAD_MADLIBDIR) -l$(THAD_MADMODELLIB)  -L$(ATHAD_MADLIBDIR) -l$(ATHAD_MADMODELLIB) -L$(TLEP_MADLIBDIR) -l$(TLEP_MADMODELLIB) -L$(H2L2NU_MADLIBDIR) -l$(H2L2NU_MADMODELLIB) -L$(TTLL_MADLIBDIR) -l$(TTLL_MADMODELLIB) ${LHAPDFLIB}
#------------------------------------------------------------------------------
# Compilation

all : ../test/test

clean:
	rm *.o
	rm test*

objects=${TTH_MADPROCDIR}/CPPProcess.o ${THAD_MADPROCDIR}/CPPProcess_tbwjj.o ${ATHAD_MADPROCDIR}/CPPProcess_antitbwjj.o ${TLEP_MADPROCDIR}/CPPProcess_tbwlnu.o ${H2L2NU_MADPROCDIR}/CPPProcess_hw2l2nu.o ${TTLL_MADPROCDIR}/CPPProcess_ggttll.o TestMEM.o

../test/test: $(objects) 
	$(CXX) -o ../test/test $(objects) ${LIBS}  $(CXXFLAGS) ${SOFLAGS}
#	$(CXX) -o test $(objects) ${ROOTLIBS} ${ROOTGLIBS} $(CXXFLAGS) ${SOFLAGS} 


TestMEM.o: TestMEM.cpp HypIntegrator.h MEPhaseSpace.h MultiLepton.h $(TTH_MADLIBDIR)/lib$(TTH_MADMODELLIB).a $(THAD_MADLIBDIR)/lib$(THAD_MADMODELLIB).a $(ATHAD_MADLIBDIR)/lib$(ATHAD_MADMODELLIB).a $(TLEP_MADLIBDIR)/lib$(TLEP_MADMODELLIB).a $(H2L2NU_MADLIBDIR)/lib$(H2L2NU_MADMODELLIB).a $(TTLL_MADLIBDIR)/lib$(TTLL_MADMODELLIB).a 
	$(CXX) ${LIBS} ${CXXFLAGS} -c TestMEM.cpp 


${TTH_MADPROCDIR}/CPPProcess.o: ${TTH_MADPROCDIR}/CPPProcess.cc $(TTH_MADLIBDIR)/lib$(TTH_MADMODELLIB).a
	$(CXX) -l$(LIBS) ${CXXFLAGS} -c ${TTH_MADPROCDIR}/CPPProcess.cc 

${THAD_MADPROCDIR}/CPPProcess_tbwjj.o: ${THAD_MADPROCDIR}/CPPProcess_tbwjj.cc $(THAD_MADLIBDIR)/lib$(THAD_MADMODELLIB).a
	$(CXX) -l$(LIBS) ${CXXFLAGS} -c ${THAD_MADPROCDIR}/CPPProcess_tbwjj.cc

${ATHAD_MADPROCDIR}/CPPProcess_antitbwjj.o: ${ATHAD_MADPROCDIR}/CPPProcess_antitbwjj.cc $(ATHAD_MADLIBDIR)/lib$(ATHAD_MADMODELLIB).a
	$(CXX) -l$(LIBS) ${CXXFLAGS} -c ${ATHAD_MADPROCDIR}/CPPProcess_antitbwjj.cc

${TLEP_MADPROCDIR}/CPPProcess_tbwlnu.o: ${TLEP_MADPROCDIR}/CPPProcess_tbwlnu.cc $(TLEP_MADLIBDIR)/lib$(TLEP_MADMODELLIB).a
	$(CXX) -l$(LIBS) ${CXXFLAGS} -c ${TLEP_MADPROCDIR}/CPPProcess_tbwlnu.cc

${H2L2NU_MADPROCDIR}/CPPProcess_hw2l2nu.o: ${H2L2NU_MADPROCDIR}/CPPProcess_hw2l2nu.cc $(H2L2NU_MADLIBDIR)/lib$(H2L2NU_MADMODELLIB).a
	$(CXX) -l$(LIBS) ${CXXFLAGS} -c ${H2L2NU_MADPROCDIR}/CPPProcess_hw2l2nu.cc

${TTLL_MADPROCDIR}/CPPProcess_ggttll.o: ${TTLL_MADPROCDIR}/CPPProcess_ggttll.cc $(TTLL_MADLIBDIR)/lib$(TTLL_MADMODELLIB).a
	$(CXX) -l$(LIBS) ${CXXFLAGS} -c ${TTLL_MADPROCDIR}/CPPProcess_ggttll.cc


$(TTH_MADLIBDIR)/lib$(TTH_MADMODELLIB).a:
	cd ../Madgraph/PROC_SA_CPP_sm_4/src
	make
	cd -

$(THAD_MADLIBDIR)/lib$(THAD_MADMODELLIB).a:
	cd ../Madgraph/PROC_SA_CPP_sm_DECAY_tbwjj/src
	make
	cd -

$(TLEP_MADLIBDIR)/lib$(TLEP_MADMODELLIB).a:
	cd ../Madgraph/PROC_SA_CPP_sm_DECAY_tbwlnu/src
	make
	cd -

$(H2L2NU_MADLIBDIR)/lib$(H2L2NU_MADMODELLIB).a:
	cd ../Madgraph/PROC_SA_CPP_sm_DECAY_hw2l2nu/src
	make
	cd -

$(TTLL_MADLIBDIR)/lib$(TTLL_MADMODELLIB).a:
	cd ../Madgraph/PROC_SA_CPP_sm_4/src
	make
	cd -

#../gosam/ttH_ME/libgolem_olp.so:
#	cd ../../gosam/ttH_ME
#	make



