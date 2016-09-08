#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <deque>
#include <list>
#include <random>
#include <algorithm>
#include <type_traits>
#include <mutex>
#include <condition_variable>
#include <omp.h>

#ifndef NOINLINE
#define NOINLINE
#endif

#include "genetic_bits/internal.hpp"
#include "genetic_bits/rw_semaphore.hpp"
#include "genetic_bits/CTIterator.hpp"
#include "genetic_bits/globals.hpp"
#include "genetic_bits/Candidate.hpp"

#include "genetic_bits/CandidateTagged.hpp"
#include "genetic_bits/BasePopulation.hpp"
#include "genetic_bits/OrdPopulation.hpp"
#include "genetic_bits/FloatPopulation.hpp"
#include "genetic_bits/DomPopulation.hpp"
#include "genetic_bits/NSGAPopulation.hpp"
#include "genetic_bits/Population.hpp"

#endif
