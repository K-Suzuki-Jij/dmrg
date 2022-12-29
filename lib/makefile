SMLDIR      := sml
DMRGDIR     := dmrg
EDDIR       := exact_diag
INCDIR      := ../include
MAKLMDIR    := model/AKLM
MAKLTVFMDIR := model/AKLM_TVF
MEKLMDIR    := model/EKLM
MHUBDIR     := model/HUBBARD
MTAKLMDIR   := model/TAKLM
MXXZDIR     := model/XXZ

# MacOS clang
CC = g++
CFLAGS = -Wall -Wextra -std=c++17 -Ofast -Xpreprocessor -fopenmp -I $(INCDIR)

# Linux
#CC = g++
#CFLAGS = -Wall -Wextra -std=c++17 -Ofast -fopenmp -I $(INCDIR)

# Linux OneAPI
#CC = icx
#CFLAGS = -Wall -std=c++17 -Ofast -fiopenmp -I $(INCDIR)

.PHONY: all sml dmrg model exact clean

all: sml dmrg model exact

sml  : libsml.a
dmrg : libdmrg.a
exact: libed.a
model: libmaklm.a libmaklm_tvf.a libmeklm.a libmhubbard.a libmtaklm.a libmxxz.a

###SML
libsml.a: $(addprefix $(SMLDIR)/,$(patsubst %.cpp,%.o,$(wildcard ../$(SMLDIR)/*.cpp)))
	ar r $@ $^
	
###Model
libmaklm.a: $(addprefix $(MAKLMDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MAKLMDIR)/*.cpp))))
	ar r $@ $^
	
libmaklm_tvf.a: $(addprefix $(MAKLTVFMDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MAKLTVFMDIR)/*.cpp))))
	ar r $@ $^
	
libmeklm.a: $(addprefix $(MEKLMDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MEKLMDIR)/*.cpp))))
	ar r $@ $^

libmhubbard.a: $(addprefix $(MHUBDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MHUBDIR)/*.cpp))))
	ar r $@ $^
	
libmtaklm.a: $(addprefix $(MTAKLMDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MTAKLMDIR)/*.cpp))))
	ar r $@ $^
	
libmxxz.a: $(addprefix $(MXXZDIR)/,$(notdir $(patsubst %.cpp,%.o,$(wildcard ../$(MXXZDIR)/*.cpp))))
	ar r $@ $^

###DMRG
libdmrg.a: $(addprefix $(DMRGDIR)/,$(patsubst %.cpp,%.o,$(wildcard ../$(DMRGDIR)/*.cpp)))
	ar r $@ $^
	
libed.a: $(addprefix $(EDDIR)/,$(patsubst %.cpp,%.o,$(wildcard ../$(EDDIR)/*.cpp)))
	ar r $@ $^
	
$(SMLDIR)/%.o: ../$(SMLDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(SMLDIR)/$(@F) -c ../$(SMLDIR)/$(<F)
	
$(MAKLMDIR)/%.o: ../$(MAKLMDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MAKLMDIR)/$(@F) -c ../$(MAKLMDIR)/$(<F)
	
$(MAKLTVFMDIR)/%.o: ../$(MAKLTVFMDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MAKLTVFMDIR)/$(@F) -c ../$(MAKLTVFMDIR)/$(<F)
	
$(MEKLMDIR)/%.o: ../$(MEKLMDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MEKLMDIR)/$(@F) -c ../$(MEKLMDIR)/$(<F)
	
$(MHUBDIR)/%.o: ../$(MHUBDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MHUBDIR)/$(@F) -c ../$(MHUBDIR)/$(<F)
	
$(MTAKLMDIR)/%.o: ../$(MTAKLMDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MTAKLMDIR)/$(@F) -c ../$(MTAKLMDIR)/$(<F)
	
$(MXXZDIR)/%.o: ../$(MXXZDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(MXXZDIR)/$(@F) -c ../$(MXXZDIR)/$(<F)

$(DMRGDIR)/%.o: ../$(DMRGDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(DMRGDIR)/$(@F) -c ../$(DMRGDIR)/$(<F)
	
$(EDDIR)/%.o: ../$(EDDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(EDDIR)/$(@F) -c ../$(EDDIR)/$(<F)
	
clean:
	rm -rf $(SMLDIR) model $(DMRGDIR) $(EDDIR)
	rm *.a

