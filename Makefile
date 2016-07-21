PROGS = regex binary
HEADERS = genetic.h

CXXFLAGS += -Wall -fno-diagnostics-show-caret -std=c++11 -pthread

all: CXXFLAGS += -O3
bench: CXXFLAGS += -DBENCH -O3
debug: CXXFLAGS += -g
optdebug: CXXFLAGS += -O3 -g
profile: CXXFLAGS += -O1 -pg
optprofile: CXXFLAGS += -O3 -pg

all bench debug optdebug profile optprofile: $(PROGS)

clean:
	rm $(PROGS)

.PHONY: all debug clean

$(PROGS): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@
