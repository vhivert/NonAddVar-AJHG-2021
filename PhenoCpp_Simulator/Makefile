
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)

	CFLAGS = -std=c++11 -Wall -Wextra -fopenmp -ggdb
	CC = g++
endif
ifeq ($(UNAME), Darwin)
        
	CFLAGS = -std=c++17 -Wall -Wextra -ggdb
	CC = clang++ -Xpreprocessor -fopenmp -lomp -I"/usr/local/opt/libomp/include" -L"/usr/local/opt/libomp/lib"
	LIBS = /usr/local/lib/libomp.dylib
endif

SRCDIR := src
BUILDDIR := build
TARGET := bin/PhenoSim



$(BUILDDIR)/rand.o: $(SRCDIR)/rand.cpp $(SRCDIR)/rand.h
	mkdir -p $(BUILDDIR)
	$(CC) -o $(BUILDDIR)/rand.o -c $(SRCDIR)/rand.cpp $(LIBS) $(CFLAGS) $(OPTS)
	
$(BUILDDIR)/read.o: $(SRCDIR)/read.cpp $(SRCDIR)/read.h
	$(CC) -o $(BUILDDIR)/read.o -c $(SRCDIR)/read.cpp $(LIBS) $(CFLAGS) $(OPTS)
	
$(BUILDDIR)/simul.o: $(SRCDIR)/simul.cpp $(SRCDIR)/simul.h
	$(CC) -o $(BUILDDIR)/simul.o -c $(SRCDIR)/simul.cpp $(LIBS) $(CFLAGS) $(OPTS)
	
$(BUILDDIR)/mem.o: $(SRCDIR)/mem.cpp $(SRCDIR)/mem.h
	$(CC) -o $(BUILDDIR)/mem.o -c $(SRCDIR)/mem.cpp $(LIBS) $(CFLAGS) $(OPTS)

$(BUILDDIR)/write.o: $(SRCDIR)/write.cpp $(SRCDIR)/write.h
	$(CC) -o $(BUILDDIR)/write.o -c $(SRCDIR)/write.cpp $(LIBS) $(CFLAGS) $(OPTS)

all: export OPTS=-O3

fulldebug: export OPTS=-g -O0

debug: export OPTS=-g -O2

fastdebug: export OPTS=-g -O3

profiling: export OPTS=-pg -g

all fulldebug debug fastdebug profiling: $(BUILDDIR)/rand.o $(BUILDDIR)/read.o $(BUILDDIR)/mem.o $(BUILDDIR)/write.o $(BUILDDIR)/simul.o
	$(CC) -o $(TARGET) $(SRCDIR)/PhenoCppSimulator.cpp $(BUILDDIR)/rand.o $(BUILDDIR)/read.o $(BUILDDIR)/mem.o $(BUILDDIR)/write.o $(BUILDDIR)/simul.o $(LIBS) $(CFLAGS) $(OPTS)
	
clean:
	rm -rf $(BUILDDIR) $(TARGET)*
