#include <vector>
#include <random>
#include <algorithm>


/* Exactly one source file needs to define these quantities
 * (more can be declared there). */

namespace Context {
  extern std::mt19937 rng;
};

namespace Config {
  /* > 0. This tunes rank-based selection. Zero would mean no account on
   * fitness in the selection process whatsoever. The bigger the value the
   * more candidates with low fitness are likely to be selected. */
  extern const float selectBias;

  /* > 0. In this header file, this only serves the purpose of a default value
   * to Population::trim(). */
  extern const int popSize;
};

/* Example:
 * std::mt19937 Context::rng;
 * 
 * namespace Config {
 *   const float selectBias = 1.0;
 *   const float pCrossover = 0.7;
 *   const int popSize = 50;
 *   const int popSize2 = 30;
 *   const int nGen = 50;
 * }
 *
 * in main():
 * 
 * Context::rng = std::mt19937((std::random_device())());
 */


/* The Candidate template. Needs a typename Fitness which can be a simple type
 * or a class but needs to implement operator<. The virtual implementation
 * only keeps track of whether fitness has been computed, and provides the
 * precomputed value when available. The inner function to compute fitness
 * when required, computeFitness(), needs to be provided in specializations. */
template<typename Fitness>
class ICandidate {
  mutable Fitness _fitness;
  mutable bool fitnessValid = false;

  public:
  Fitness fitness() const {
    if(!fitnessValid) {
      _fitness = computeFitness();
      fitnessValid = true;
    }
    return _fitness;
  }

  friend bool operator< (const ICandidate& c1, const ICandidate& c2) {
    return c1.fitness() < c2.fitness();
  }

  private:
  virtual Fitness computeFitness() const = 0;
};


/* The Population template. Requires a Candidate class which is expected to be
 * derived from ICandidate although no exposed properties except operator <
 * are accessed. Default constructor creates an empty population. */
template<class Candidate>
class Population : private std::vector<Candidate> {
  bool sorted = false;

  public:
  /* Creates an empty population. */
  Population() = default;

  /* Draws a population from a source function. */
  Population(size_t size, std::function<Candidate()> src) {
    add(size, src);
  }

  /* Pushes back a new candidate. */
  void add(const Candidate& c) {
    this->push_back(c);
    sorted = false;
  }

  /* Draws n candidates from a source function. */
  void add(size_t size, std::function<Candidate()> src) {
    for(size_t j = 0; j < size; j++)
      add(src());
  }


  using std::vector<Candidate>::begin;
  using std::vector<Candidate>::end;
  using std::vector<Candidate>::size;
  using std::vector<Candidate>::clear;

  /* Retrieves a candidate randomly chosen by rank-based selection,
   * Config::selectBias determines how much low-fitness solutions are
   * preferred, see discussion in Config above. */
  Candidate& rankSelect(float bias = Config::selectBias) {
    ensureSorted();
    float x = (std::uniform_real_distribution<float>(0, 1))(Context::rng);
    return (*this)[(int)(-log(1 - x + x*exp(-bias))/bias*this->size())];
  }

  /* Returns unconditionally the best candidate of population. If more
   * candidates have equal best fitness the returned reference may be any of
   * them. */
  Candidate& best() {
    ensureSorted();
    return *this->begin();
  }

  /* Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample. */
  void trim(size_t size = Config::popSize) {
    std::sort(this->begin(), this->end());
    sorted = true;
    if(this->size() > size)
      this->resize(size);
  }

  private:
  inline void ensureSorted() {
    if(!sorted) {
      std::sort(this->begin(), this->end());
      sorted = true;
    }
  }
};


/* This derived class provides a new function stat() for candidate classes
 * whose fitness is a simple floating point type or allows an implicit
 * convertion to one. */
template<class Candidate, typename Float>
class SPopulation : public Population<Candidate> {
  public:
    struct Stat {
      Float mean;
      Float stdev;
    };
   
    /* Calculates the mean and standard deviation of the sample's fitnesses. */
    Stat stat() {
      Float sf = 0, sf2 = 0;
      for(Candidate c : *this) {
        Float f = c.fitness();
        sf += f;
        sf2 += f*f;
      }
      size_t sz = this->size();
      Float dev2 = sf2/sz - sf/sz*sf/sz;
      return Stat{sf/sz, dev2 >= 0 ? (Float)sqrt(dev2) : 0};
    }
};
