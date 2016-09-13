PROGS = regex binary binary-mo
HEADERS := $(addprefix ../,$(HEADERS))

# This allows to distinguish between binaries made for debug, profiling, etc.
# The main purpose is that "make debug" doesn't think there's nothing to make 
# after just "make", but also separate targets may be specified, like binary.d.
SFXS = b d p od op p2 gp gpb
ALLTARGETS = $(PROGS) $(foreach SFX,$(SFXS),$(addsuffix .$(SFX),$(PROGS)))

# This allows to compile all programs e.g. for profiling.
# Suffixes need to use . so that they can be stripped using $(basename).
all: SFX =
bench: SFX = .b
debug: SFX = .d
optdebug: SFX = .od
profile: SFX = .p
optprofile: SFX = .op
gpmulti: SFX = .gp
gpsingle: SFX = .gpb

# Common C++ flags
CXXFLAGS += -I../include
CXXFLAGS += -std=c++11 -march=native
CXXFLAGS += -pedantic -Wall -Wextra -Weffc++
CXXFLAGS += -fno-diagnostics-show-caret -fopenmp

# Flags specific for different targets
# DEBUG: smaller number of iterations, less inlining
# BENCH: single thread and repeatable RNG
FLAGS     += -O3
FLAGS.b   += -O3 -DDEBUG -DBENCH
FLAGS.d   += -O1 -DDEBUG -DBENCH -g
FLAGS.od  += -O3 -DDEBUG         -g -fno-inline-functions
FLAGS.p   += -O1 -DDEBUG -DBENCH -pg
FLAGS.op  += -O3 -DDEBUG         -pg -fno-inline-functions
FLAGS.gp  += -O3 -DDEBUG         -lprofiler
FLAGS.gpb += -O3 -DDEBUG -DBENCH -lprofiler


# SECONDEXPANSION is necessary because otherwise $(SFX) wouldn't be known at 
# this point. Also allows $(basename %) below.
.SECONDEXPANSION:
all bench debug optdebug profile optprofile: $$(addsuffix $$(SFX),$$(PROGS))

clean:
	rm -f $(ALLTARGETS)

.PHONY: all clean bench debug optdebug profile optprofile


# basename: binary.d -> binary
# $(patsubst $(basename $@)%,%,$@): binary.d -> .d
# $(FLAGS...): binary.d -> $(FLAGS.d)

$(ALLTARGETS): %: $$(basename %).cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(FLAGS$(patsubst $(basename $@)%,%,$@)) $< -o $@
