PROGNAME := 1D_XXZ_ED.out
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

# MacOS clang
#CC = g++
#CFLAGS = -Wall -std=c++17 -O3 -mtune=native -march=native -Xpreprocessor -fopenmp -I $(INCDIR)
#LFLAGS = -lmxxz -led -lsml -llapack -lomp -L $(LIBDIR)

CC = g++
CFLAGS = -Wall -std=c++17 -O3 -fopenmp -I $(INCDIR)
LFLAGS = -lmxxz -led -lsml -llapack -L $(LIBDIR)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(OUTDIR)/%.o: %.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR)
