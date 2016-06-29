# Makefile

CC = g++ -g -Wno-writable-strings -Wno-write-strings -fsanitize=address -lTMVA

NTUPLEINC    = -I../NtupleProducer/include -I$(CMSSW_RELEASE_BASE)/src
NTUPLELIB    = -L../NtupleProducer -L$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc530 -lCondFormatsJetMETObjects \
              -lNtuple

KINFITINC    = -I../../KinFit
KINFITLIB    = -L../../KinFit \
              -lKinFit

SOVERSION               = $(shell echo $(VERSION) | \
                               sed 's/\(.*\)\.\(.*\)\.\(.*\)/\1:\2:0/')
ROOTCONF                = root-config
ROOTCINT                = $(shell $(ROOTCONF) --prefix)/bin/rootcint
ROOTLIBS                = $(shell $(ROOTCONF) --libs)
ROOTGLIBS               = $(shell $(ROOTCONF) --glibs)
ROOTCFLAGS              = $(shell $(ROOTCONF) --cflags)
INCLUDES                = $(shell $(ROOTCONF) --incdir) $(ROOTCFLAGS)
CINTINCLUDES            = $(shell $(ROOTCONF) --incdir)
LIBS                    = $(ROOTLIBS)

LDLIBS  = $(ROOTLIBS) $(NTUPLELIB) $(KINFITLIB) 

ADDINC = $(NTUPLEINC) $(KINFITINC) 

SRC = src/Analyzer.cxx  src/Hist.cxx \
src/TTbarDileptonAnalysis.cxx src/HistoManager.cxx \
src/TTbarHiggsMultileptonAnalysis.cxx src/TTbarHiggsTFAnalysis.cxx \
src/Lepton.cxx \
src/TTbarHiggsTriggerAnalysis.cxx


INC = include/Hist.h include/Ranges.h include/Analyzer.h \
      include/TTbarDileptonAnalysis.h \
      include/HistoManager.h include/TTbarHiggsMultileptonAnalysis.h \
      include/TTbarHiggsTFAnalysis.h \
      include/Lepton.h \
      include/TTbarHiggsTriggerAnalysis.h

OBJ = Analyzer.o Hist.o \
TTbarDileptonAnalysis.o HistoManager.o \
TTbarHiggsMultileptonAnalysis.o \
TTbarHiggsTFAnalysis.o \
Lepton.o \
TTbarHiggsTriggerAnalysis.o

EXE = Analyzer

all:  	     		$(EXE)
#			@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(OBJ):     		$(SRC) $(INC)
		        $(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)

clean: 		
			@rm -f $(OBJ) $(EXE)
