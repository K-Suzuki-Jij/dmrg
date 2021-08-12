SMLDIR  := sml
DMRGDIR  := dmrg
MEKLMDIR := model/EKLM
INCDIR   := ../include

CC = g++
CFLAGS = -Wall -Wextra -std=c++17 -Ofast -mtune=native -march=native -Xpreprocessor -fopenmp -I $(INCDIR)
LFLAGS = -llapack -lomp

.PHONY: all sml dmrg model clean

all: sml dmrg model

sml  : libsml.a
dmrg : libdmrg.a
model: libmeklm.a

libsml.a: $(addprefix $(SMLDIR)/,$(patsubst %.cpp,%.o,$(wildcard ../$(SMLDIR)/*.cpp)))
	ar r libsml.a $(SMLDIR)/*.o
	
libmeklm.a: $(addprefix $(MEKLMDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MEKLMDIR)/*.cpp))))
	ar r libmeklm.a $(MEKLMDIR)/*.o

libdmrg.a: $(addprefix $(DMRGDIR)/,$(patsubst %.cpp,%.o,$(wildcard ../$(DMRGDIR)/*.cpp)))
	ar r libdmrg.a $(DMRGDIR)/*.o

$(SMLDIR)/%.o: ../$(SMLDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(SMLDIR)/$(@F) -c ../$(SMLDIR)/$(<F)
	
$(MEKLMDIR)/%.o: ../$(MEKLMDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MEKLMDIR)/$(@F) -c ../$(MEKLMDIR)/$(<F)

$(DMRGDIR)/%.o: ../$(DMRGDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(DMRGDIR)/$(@F) -c ../$(DMRGDIR)/$(<F)
	
clean:
	rm -rf $(SMLDIR) model $(DMRGDIR)
	rm *.a

