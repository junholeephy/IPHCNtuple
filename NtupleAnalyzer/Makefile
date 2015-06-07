# Makefile

CC = g++ -g

NTUPLEINC    = -I$(CMSSW_BASE)/src/NtupleProducer/include
NTUPLELIB    = -L$(CMSSW_BASE)/src/NtupleProducer \
              -lNtuple

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

LDLIBS  = $(ROOTLIBS) $(NTUPLELIB)

ADDINC = $(NTUPLEINC)

SRC = src/Analyzer.cxx src/Hist.cxx src/TransferFunc.cxx \
      src/TopHiggs.cxx

INC = include/Hist.h include/Ranges.h include/Analyzer.h \
      include/TransferFunc.h include/TopHiggs.h

OBJ = Analyzer.o Hist.o TransferFunc.o TopHiggs.o

EXE = Analyzer

all:  	     		$(EXE)
#			@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(OBJ):     		$(SRC) $(INC)
		        $(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)

clean: 		
			@rm -f $(OBJ) $(EXE)
