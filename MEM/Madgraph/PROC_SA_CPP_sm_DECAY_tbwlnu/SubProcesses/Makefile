LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm_tbwlnu
CXXFLAGS= -O -I$(INCDIR) -I.
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

main=check

all: $(main)

objects=CPPProcess_tbwlnu.o check_sa.o

$(main): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -o $@ $(objects) $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

