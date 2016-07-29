PROGS = regex binary
HEADERS = genetic.h

CXXFLAGS += -pedantic -Wall -Wextra -Weffc++
CXXFLAGS += -fno-diagnostics-show-caret -pthread

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
	$(CXX) -std=c++11 $(CXXFLAGS) $@.cpp -o $@
