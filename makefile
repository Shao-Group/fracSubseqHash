CC=gcc
CPP=g++
CFLAGS+= -m64 -g -Wall
LDFLAGS= -L$$GUROBI_HOME/lib -lgurobi91
LIBS=
INC= $$GUROBI_HOME/include/
ALLDEP:= $(patsubst %.h,%.o,$(wildcard *.h)) $(wildcard *.hpp) $(wildcard *.tpp)
ALLILP:= $(wildcard *_ILP.c)

.PHONY: all
all: $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.h,%.cpp,$(wildcard *.h)), $(wildcard *.cpp)))

.PHONY: product
product: CFLAGS = -O3
product: $(ALLDEP) $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.o,%.cpp,$(ALLDEP)) $(ALLILP), $(wildcard *.cpp)))

sampleFast%.out: sampleFast%.cpp
	$(CPP) $(CFLAGS) -std=c++17 -o $@ $^ $(LIBS)
genSubseqSeed%.out: genSubseqSeed%.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $^ $(LIBS) -pthread

%.out: %.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $(filter-out $(wildcard *.tpp), $^) $(LIBS)
%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.cpp %.h makefile
	$(CPP) $(CFLAGS) -MMD -c $< -o $@
%.o: %.c %.h makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

#gurobi make
%_ILP: %_ILP.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ -I$(INC) $(LDFLAGS) -lm

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     
