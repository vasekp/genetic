#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <random>
#include <algorithm>
#include <type_traits>
#include <mutex>

namespace gen {

namespace internal {
  /** \brief The default random number generator for Population::rankSelect(). */
  static thread_local std::ranlux48_base rng{std::random_device()()};
}

/** \brief The `Candidate` template.
 *
 * Needs a typename `Fitness` which can be a simple type or a class but needs
 * to implement zero-argument construction and `operator<()`.  The virtual
 * implementation only keeps track of whether fitness has been computed, and
 * provides the precomputed value when available. The inner function to
 * compute fitness when required, computeFitness(), needs to be provided in
 * specializations. */
template<typename Fitness>
class ICandidate {
  mutable Fitness _fitness{};
  mutable bool fitnessValid = false;

  public:
  /** \brief The Fitness type provided for this template specialization. */
  typedef Fitness _FitnessType;

  /** \brief Returns this `Candidate`'s fitness, calculating it on request if not
   * known from before. */
  Fitness fitness() const {
    if(!fitnessValid) {
      _fitness = computeFitness();
      fitnessValid = true;
    }
    return _fitness;
  }

  /** \brief Compares two `Candidate`s by the Fitness's `operator<`. */
  friend bool operator< (const ICandidate<Fitness>& c1, const ICandidate<Fitness>& c2) {
    return c1.fitness() < c2.fitness();
  }

  virtual ~ICandidate() { }

  protected:
  /** \brief The internal fitness computation, called the first time this
   * `Candidate`'s fitness() is queried. Every specialization must implement
   * this routine. */
  virtual Fitness computeFitness() const = 0;
}; // class ICandidate


/** \brief The Population template.
 *
 * Requires a Candidate class which is expected to be derived from ICandidate
 * although no exposed properties except `operator<` are accessed. */
template<class Candidate>
class Population : private std::vector<Candidate> {
  bool sorted = false;
  std::mutex mtx{};

  public:
  /** \brief Creates an empty population. */
  Population() = default;

  /** \brief Creates an empty population but preallocate space for count
   * candidates. */
  Population(size_t count) {
    this->reserve(count);
  }

  /** \brief Creates a population of size `count` whose candidates are results
   * of calls to the source function `src`. For discussion about the latter
   * parameter see add(size_t, Source).
   *
   * \see add(size_t, Source) */
  template<class Source>
  Population(size_t count, Source src) {
    add(count, src);
  }

  /* The Big Four: trivial but we need them because the mutex can't be 
   * default copied or moved */

  /** \brief The copy constructor. */
  Population(const Population& _p): std::vector<Candidate>(_p), sorted(_p.sorted) { }

  /** \brief The move constructor. */
  Population(Population&& _p): std::vector<Candidate>(std::move(_p)), sorted(_p.sorted) { }

  /** \brief Copy assignment operator. */
  Population& operator=(const Population& _p) {
    std::lock_guard<std::mutex> lock(mtx);
    std::vector<Candidate>::operator=(_p);
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Move assignment operator. */
  Population& operator=(Population&& _p) {
    std::lock_guard<std::mutex> lock(mtx);
    std::vector<Candidate>::operator=(std::move(_p));
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Adds a new candidate. */
  void add(const Candidate& c) {
    std::lock_guard<std::mutex> lock(mtx);
    this->push_back(c);
    sorted = false;
  }

  /** \brief Pushes back a new candidate using the move semantics. */
  void add(Candidate&& c) {
    std::lock_guard<std::mutex> lock(mtx);
    this->push_back(std::forward<Candidate>(c));
    sorted = false;
  }

  /** \brief Draws `count` candidates from a source function `src`.
   *
   * Source can be:
   * - `std::function<Candidate>`: returning by copy,
   * - `std::function<const Candidate&>`: returning by reference,
   * - a lambda function returning either.
   *
   * The template allows for optimizations (inlining) in the latter case. */
  template<class Source>
  void add(size_t count, Source src) {
    std::lock_guard<std::mutex> lock(mtx);
    this->reserve(size() + count);
    for(size_t j = 0; j < count; j++)
      this->push_back(src());
    sorted = false;
  }

  /** \brief Takes all candidates from another population. */
  void add(Population<Candidate>& pop) {
    std::lock_guard<std::mutex> lock(mtx);
    this->reserve(size() + pop.size());
    this->insert(end(), pop.begin(), pop.end());
    sorted = false;
  }

  /** \brief Like add(Population&) but moving the contents of the argument. */
  void merge(Population<Candidate>& pop) {
    std::lock_guard<std::mutex> lock(mtx);
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

  /** \brief Retrieves a candidate randomly chosen by rank-based selection.
   *
   * bias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected. */
  template<class Rng = decltype(internal::rng)>
  const Candidate& rankSelect(float bias, Rng& rng = internal::rng) {
    static thread_local std::uniform_real_distribution<float> rDist(0, 1);
    ensureSorted();
    float x = rDist(rng);
    if(x == 1)
      return this->back();
    else
      return (*this)[(int)(-log(1 - x + x*exp(-bias))/bias*size())];
  }

  /** \brief Returns the best candidate of population.
   *
   * If more candidates have equal best fitness the returned reference may be
   * any of them. */
  const Candidate& best() {
    ensureSorted();
    return this->front();
  }

  /** \brief Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample. */
  Population<Candidate>& trim(size_t newSize) {
    ensureSorted();
    std::lock_guard<std::mutex> lock(mtx);
    if(size() > newSize)
      this->resize(newSize);
    return *this;
  }


  /** \brief The return type of stat(). */
  struct Stat {
    double mean;  ///< The mean fitness of the Population.
    double stdev; ///< The standard deviation of fitness in the Population.
  };

  /** \brief Returns the mean fitness of the population and the standard deviation.
   *
   * Conditional member (using SFINAE) for candidate classes whose fitness is
   * a simple floating point type or allows an implicit convertion to one. The
   * method does not appear in specializations for which this condition is not
   * satisfied. */
#ifdef DOXYGEN
  Stat stat() {
#else
  template<typename FT = typename Candidate::_FitnessType>
  auto stat() -> typename std::enable_if<std::is_convertible<FT, double>::value, Stat>::type {
#endif
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

  private:
  inline void ensureSorted() {
    std::lock_guard<std::mutex> lock(mtx);
    if(!sorted) {
      std::sort(begin(), end());
      sorted = true;
    }
  }
}; // class Population

} // namespace gen

#endif
