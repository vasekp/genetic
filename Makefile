PROGS = regex binary
HEADERS = genetic.h

CXXFLAGS += -Wall -fno-diagnostics-show-caret -std=c++11 -pthread

all: CXXFLAGS += -O3
debug: CXXFLAGS += -g
all debug: $(PROGS)

clean:
	rm $(PROGS)

.PHONY: all debug clean

$(PROGS): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@
