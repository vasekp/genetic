PROGS = regex binary
HEADERS = genetic.h

GCC = g++
OPTIONS = -Wall -fno-diagnostics-show-caret -std=c++11

all: $(PROGS)

$(PROGS): %: %.cpp $(HEADERS)
	$(GCC) $(OPTIONS) -O3 $@.cpp -o $@
