PROGNAME := 1D_AKLM_DMRG.out
SRCDIR   := src
OUTDIR   := build
INCDIR   := ../../../include
LIBDIR   := ../../../lib
TARGET   := $(OUTDIR)/$(PROGNAME)
SRCS     := $(wildcard $(SRCDIR)/*.cpp)
OBJS     := $(addprefix $(OUTDIR)/,$(patsubst %.cpp,%.o,$(SRCS)))

#CC = icpc
#CFLAGS = -Wall -std=c++17 -O3 -qopenmp -parallel -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I $(INCDIR)
#LFLAGS = -mkl

CC = g++
CFLAGS = -Wall -std=c++17 -O3 -mtune=native -march=native -Xpreprocessor -fopenmp -I $(INCDIR)
LFLAGS = -llapack -lomp -lsml -ldmrg -lmaklm -L $(LIBDIR)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(OUTDIR)/%.o: %.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR)
