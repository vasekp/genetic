PROGS = regex binary
HEADERS = genetic.h

CXXFLAGS += -Wall -fno-diagnostics-show-caret -std=c++11 -pthread

all: CXXFLAGS += -O3
bench: CXXFLAGS += -DBENCH -O3
debug: CXXFLAGS += -DBENCH -g
optdebug: CXXFLAGS += -O3 -DBENCH -fno-inline-functions -g
profile: CXXFLAGS += -O1 -DBENCH -pg
optprofile: CXXFLAGS += -O3 -DBENCH -fno-inline-functions -pg

all bench debug optdebug profile optprofile: $(PROGS)

clean:
	rm $(PROGS)

.PHONY: all debug clean

$(PROGS): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@
