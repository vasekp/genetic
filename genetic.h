#include <vector>
#include <random>
#include <algorithm>


/* Exactly one source file needs to define these quantities
 * (more can be declared there). */

namespace Context {
  extern std::mt19937 rng;
};

namespace Config {
  /* This tunes rank-based selection. Zero would mean no account on
   * fitness in the selection process whatsoever. The bigger the value the
   * more candidates with low fitness are likely to be selected. A negative
   * value can be used if high fitness candidates are to be sought. */
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
class Population : public std::vector<Candidate> {
  bool sorted = false;

  public:
  /* Pushes back a new candidate. */
  void add(const Candidate& c) {
    this->push_back(c);
    sorted = false;
  }

  /* Retrieves a candidate randomly chosen by rank-based selection,
   * Config::selectBias determines how much low-fitness solutions are
   * preferred, see discussion in Config above. */
  Candidate& rankSelect() {
    float x = (std::uniform_real_distribution<float>(0, 1))(Context::rng);
    return (*this)[(int)(-log(1 - x + x*exp(-Config::selectBias))/Config::selectBias*this->size())];
  }

  /* Returns unconditionally the best candidate of population. If more
   * candidates have equal best fitness the returned reference may be any of
   * them. */
  Candidate& best() {
    if(!sorted) {
      std::sort(this->begin(), this->end());
      sorted = true;
    }
    return *this->begin();
  }

  /* Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample. */
  Population& trim(size_t size = Config::popSize) {
    std::sort(this->begin(), this->end());
    sorted = true;
    if(this->size() > size)
      this->resize(size);
    return *this;
  }
};
