#include <vector>
#include <random>
#include <algorithm>
#include <type_traits>


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
  typedef Fitness _FitnessType;

  /* Returns this Candidate's fitness, calculating it on request if not known
   * from before. */
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
  /* Every Candidate class must implement this routine. */
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

  /* Create an empty population but preallocate space for count candidates. */
  Population(size_t count) {
    this->reserve(count);
  }

  /* Draws a population from a source function. */
  template<class F>
  Population(size_t count, F src) {
    add(count, src);
  }

  /* Pushes back a new candidate. */
  void add(const Candidate& c) {
    this->push_back(c);
    sorted = false;
  }

  /* Pushes back a new candidate using the move semantics. */
  void add(Candidate&& c) {
    this->push_back(std::forward<Candidate>(c));
    sorted = false;
  }

  /* Draws n candidates from a source function
   * F can be:
   * - std::function<Candidate>: returning by copy
   * - std::function<const Candidate&>: returning by reference
   * - lambda returning either
   * The template allows for optimizations (inlining) in the latter case. */
  template<class F>
  void add(size_t count, F src) {
    this->reserve(size() + count);
    for(size_t j = 0; j < count; j++)
      this->push_back(src());
    sorted = false;
  }

  /* Draws n candidates from a source function (returning reference). */
  void add(size_t count, std::function<const Candidate&()> src) {
    this->reserve(size() + count);
    for(size_t j = 0; j < count; j++)
      this->push_back(src());
    sorted = false;
  }

  /* Takes all candidates from another population. */
  void add(Population<Candidate>& pop) {
    this->reserve(size() + pop.size());
    this->insert(end(), pop.begin(), pop.end());
    sorted = false;
  }

  /* Like add(Population&) but moving the contents of the argument. */
  void merge(Population<Candidate>& pop) {
    this->reserve(size() + pop.size());
    this->insert(end(), std::make_move_iterator(pop.begin()), std::make_move_iterator(pop.end()));
    pop.clear();
    sorted = false;
  }

  using std::vector<Candidate>::begin;
  using std::vector<Candidate>::end;
  using std::vector<Candidate>::size;
  using std::vector<Candidate>::clear;
  using std::vector<Candidate>::operator[];

  /* Retrieves a candidate randomly chosen by rank-based selection.
   * selectBias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected. */
  template<class Rng>
  const Candidate& rankSelect(Rng& rng, float bias) {
    static thread_local std::uniform_real_distribution<float> rDist(0, 1);
    ensureSorted();
    float x = rDist(rng);
    if(x == 1)
      return this->back();
    else
      return (*this)[(int)(-log(1 - x + x*exp(-bias))/bias*size())];
  }

  /* Returns unconditionally the best candidate of population. If more
   * candidates have equal best fitness the returned reference may be any of
   * them. */
  const Candidate& best() {
    ensureSorted();
    return this->front();
  }

  /* Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample. */
  void trim(size_t newSize) {
    std::sort(begin(), end());
    sorted = true;
    if(size() > newSize)
      this->resize(newSize);
  }


  struct Stat {
    double mean;
    double stdev;
  };

  /* Returns the mean fitness of the population and the standard deviation.
   * Conditional member (using SFINAE) for candidate classes whose fitness is
   * a simple floating point type or allows an implicit convertion to one. */
  template<typename FT = typename Candidate::_FitnessType>
  auto stat() -> typename std::enable_if<std::is_convertible<FT, double>::value, Stat>::type {
    double f, sf = 0, sf2 = 0;
    for(Candidate &c : *this) {
      f = c.fitness();
      sf += f;
      sf2 += f*f;
    }
    size_t sz = size();
    double dev2 = sf2/sz - sf/sz*sf/sz;
    return Stat{sf/sz, dev2 >= 0 ? sqrt(dev2) : 0};
  }

  inline void precomputeFitnesses() {
    for(Candidate &c : *this)
      c.fitness();
  }

  inline void ensureSorted() {
    if(!sorted) {
      std::sort(begin(), end());
      sorted = true;
    }
  }
};
