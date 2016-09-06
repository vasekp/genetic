PROGS = regex binary binary-mo
HEADERS = include/genetic include/genetic_bits/*.hpp

# This allows to distinguish between binaries made for debug, profiling, etc.
# The main purpose is that "make debug" doesn't think there's nothing to make 
# after just "make", but also separate targets may be specified, like binary.d.
SFXS = b d p od op
ALLTARGETS = $(PROGS) $(foreach SFX,$(SFXS),$(addsuffix .$(SFX),$(PROGS)))

# This allows to compile all programs e.g. for profiling.
# Suffixes need to use . so that they can be stripped using $(basename).
all: SFX =
bench: SFX = .b
debug: SFX = .d
optdebug: SFX = .od
profile: SFX = .p
optprofile: SFX = .op

# Common C++ flags
CXXFLAGS += -Iinclude
CXXFLAGS += -std=c++11 -march=native
CXXFLAGS += -pedantic -Wall -Wextra -Weffc++
CXXFLAGS += -fno-diagnostics-show-caret -fopenmp

# Flags specific for different targets
FLAGS += -O3
FLAGS.b += -DBENCH -DDEBUG -O3
FLAGS.d += -DBENCH -DDEBUG -g
FLAGS.od += -O3 -DDEBUG -fno-inline-functions -g
FLAGS.p += -O1 -DBENCH -DDEBUG -pg
FLAGS.op += -O3 -DDEBUG -fno-inline-functions -pg


# SECONDEXPANSION is necessary because otherwise $(SFX) wouldn't be known at 
# this point. Also allows $(basename %) below.
.SECONDEXPANSION:
all bench debug optdebug profile optprofile: $$(addsuffix $$(SFX),$$(PROGS))

clean:
	rm -f $(ALLTARGETS)

doc:	$(HEADERS) Doxyfile
	doxygen Doxyfile

.PHONY: all doc clean bench debug optdebug profile optprofile


# basename: binary.d -> binary
# $(patsubst $(basename $@)%,%,$@): binary.d -> .d
# $(FLAGS...): binary.d -> $(FLAGS.d)

$(ALLTARGETS): %: $$(basename %).cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(FLAGS$(patsubst $(basename $@)%,%,$@)) $< -o $@
