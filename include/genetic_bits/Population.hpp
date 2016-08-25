namespace gen {

/** \brief The Population template.
 *
 * Requires a `Candidate` class derived from ICandidate and implementing a
 * default constructor. */
template<class Candidate>
class Population : private std::vector<Candidate> {
  bool sorted = false;
  typedef std::shared_timed_mutex mutex_t;
  mutable mutex_t mtx{};

  static_assert(std::is_default_constructible<Candidate>::value,
      "The Candidate type needs to provide a default constructor.");

  typedef decltype(internal::detectFT<Candidate>(nullptr)) _FitnessType;

  static_assert(internal::hasFT<Candidate>(nullptr) &&
      std::is_base_of<ICandidate<_FitnessType>, Candidate>::value,
      "The Candidate type needs to be derived from ICandidate.");

  public:
  /** \brief Creates an empty population. */
  Population() = default;

  /** \brief Creates an empty population but preallocate space for count
   * candidates. */
  Population(size_t count) {
    this->reserve(count);
  }

  /** \brief Creates a population of size `count` whose candidates are results
   * of calls to the source function `src`.
   *
   * For discussion about the latter parameter see add(size_t, Source).
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
    std::lock_guard<mutex_t> lock(mtx);
    std::vector<Candidate>::operator=(_p);
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Move assignment operator. */
  Population& operator=(Population&& _p) {
    std::lock_guard<mutex_t> lock(mtx);
    std::vector<Candidate>::operator=(std::move(_p));
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Adds a new candidate. */
  void add(const Candidate& c) {
    std::lock_guard<mutex_t> lock(mtx);
    this->push_back(c);
    sorted = false;
  }

  /** \brief Pushes back a new candidate using the move semantics. */
  void add(Candidate&& c) {
    std::lock_guard<mutex_t> lock(mtx);
    this->push_back(std::forward<Candidate>(c));
    sorted = false;
  }

  /** \brief Draws `count` candidates from a source function `src`.
   *
   * Source can be:
   * - `std::function<(const) Candidate>`: returning by copy,
   * - `std::function<(const) Candidate&>`: returning by reference,
   * - a pointer to function returning `Candidate` or `Candidate&`,
   * - a lambda function returning either.
   *
   * The template allows for optimizations (inlining) in the latter case. */
  template<class Source>
  void NOINLINE add(size_t count, Source src) {
    std::lock_guard<mutex_t> lock(mtx);
    this->reserve(size() + count);
    for(size_t j = 0; j < count; j++)
      this->push_back(src());
    sorted = false;
  }

  /** \brief Copies all candidates from a vector of `Candidate`s. */
  void NOINLINE add(const std::vector<Candidate>& vec) {
    add(vec.begin(), vec.end());
  }

  /** \brief Copies an iterator range from a container of `Candidate`s. */
  template<class InputIt>
  void add(InputIt first, InputIt last) {
    std::lock_guard<mutex_t> lock(mtx);
    this->reserve(size() + std::distance(first, last));
    this->insert(end(), first, last);
    sorted = false;
  }

  /** \brief Moves all candidates from a vector of `Candidate`s. */
  void NOINLINE add(std::vector<Candidate>&& vec) {
    std::lock_guard<mutex_t> lock(mtx);
    this->reserve(size() + vec.size());
    this->insert(end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));
    vec.clear();
    sorted = false;
  }

  /** \brief Copies all candidates from another Population. */
  void add(const Population<Candidate>& pop) {
    return add(static_cast<const std::vector<Candidate>&>(pop));
  }

  /** \brief Moves all candidates from another Population. */
  void add(Population<Candidate>&& pop) {
    return add(static_cast<std::vector<Candidate>&&>(pop));
  }

  using std::vector<Candidate>::begin;
  using std::vector<Candidate>::end;
  using std::vector<Candidate>::size;
  using std::vector<Candidate>::clear;
  using std::vector<Candidate>::operator[];
  using std::vector<Candidate>::erase; // TODO temporary

  /** \brief Retrieves a candidate randomly chosen by rank-based selection.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double)` that will receive arguments linearly spaced between
   * `1/size * bias` and `bias` for candidates ranked `1` through `size` and
   * its return value will be interpreted as inverse probability, and as such,
   * is expected to be positive and strictly increasing in its argument. This
   * function will be built in at compile time, eliminating a function pointer
   * lookup.  The default value is `std::exp`, for which an equivalent fast
   * replacement algorithm is provided.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates an error at compile time in
   * specializations for which this condition is not satisfied.
   *
   * \param bias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   * 
   * \returns a constant reference to a randomly chosen `Candidate`. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate& NOINLINE rankSelect(float bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp(bias, rng);
    else
      return rankSelect<
        static_cast<double(*)(double, double)>(&internal::eval_in_product<fun>)
        > (bias);
  }

  /** \brief Retrieves a candidate chosen using uniform random selection. */
  template<class Rng = decltype(rng)>
  const Candidate& NOINLINE randomSelect(Rng& rng = rng) {
    std::shared_lock<mutex_t> lock(mtx);
    std::uniform_int_distribution<size_t> dist{0, size() - 1};
    return (*this)[dist(rng)];
  }

  /** \brief Retrieves a candidate randomly chosen by rank-based selection.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double, double)` that will receive its first argument linearly
   * spaced between `1/size` and `1` for candidates ranked `1` through `size`
   * and second argument equal to `bias` and its return value will be
   * interpreted as inverse probability. As such, is expected to be positive
   * and strictly increasing in its argument. This function will be built in
   * at compile time, eliminating a function pointer lookup. A usual choice
   * for `fun` is `std::pow`.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates an error at compile time in
   * specializations for which this condition is not satisfied.
   *
   * \param bias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   * 
   * \returns a constant reference to a randomly chosen `Candidate`. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate& NOINLINE rankSelect(double bias, Rng& rng = rng) {
    static thread_local std::discrete_distribution<size_t> iDist{};
    static thread_local size_t last_sz = 0;
    static thread_local std::vector<double> probs{};
    std::lock_guard<mutex_t> lock(mtx);
    size_t sz = size();
    if(sz != last_sz) {
      probs.clear();
      probs.reserve(sz);
      for(size_t i = 0; i < sz; i++)
        probs.push_back(1 / fun((double)(i+1) / sz, bias));
      iDist = std::discrete_distribution<size_t>(probs.begin(), probs.end());
      last_sz = sz;
    }
    ensureSorted();
    return (*this)[iDist(rng)];
  }

  private:
  template<class Rng>
  const Candidate& rankSelect_exp(double bias, Rng& rng = rng) {
    static thread_local std::uniform_real_distribution<double> rDist(0, 1);
    double x = rDist(rng);
    std::lock_guard<mutex_t> lock(mtx);
    ensureSorted();
    if(x == 1)
      return this->back();
    else
      return (*this)[(int)(-log(1 - x + x*exp(-bias))/bias*size())];
  }

  public:

  /** \brief Returns the best candidate of population.
   *
   * If more candidates have equal best fitness the returned reference may be
   * any of them.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates an error at compile time in
   * specializations for which this condition is not satisfied. */
  const Candidate& best() {
    std::lock_guard<mutex_t> lock(mtx);
    ensureSorted();
    return std::vector<Candidate>::front();
  }

  /** \brief Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates an error at compile time in
   * specializations for which this condition is not satisfied. */
  void trim(size_t newSize) {
    if(size() <= newSize)
      return;
    std::lock_guard<mutex_t> lock(mtx);
    ensureSorted();
    this->resize(newSize);
  }

  /** \brief Reduces the population to a maximum size given by the argument,
   * using random selection if the latter is smaller. */
  template<class Rng = decltype(rng)>
  void randomTrim(size_t newSize, Rng& rng = rng) {
    if(size() <= newSize)
      return;
    std::lock_guard<mutex_t> lock(mtx);
    std::shuffle(begin(), end(), rng);
    this->resize(newSize);
  }

  /** \brief Returns the number of `Candidate`s in this population dominated by
   * a given `Candidate`. */
  friend size_t operator<< (const Candidate& c, const Population<Candidate>& pop) {
    size_t cnt = 0;
    std::shared_lock<mutex_t> lock(pop.mtx);
    for(auto& cmp : pop)
      if(c << cmp)
        cnt++;
    return cnt;
  }

  /** \brief Returns the number of `Candidate`s in this population that
   * dominate a given `Candidate`. */
  friend size_t operator<< (const Population<Candidate>& pop, const Candidate& c) {
    size_t cnt = 0;
    std::shared_lock<mutex_t> lock(pop.mtx);
    for(auto& cmp : pop)
      if(cmp << c)
        cnt++;
    return cnt;
  }

  /** \brief Returns a nondominated subset of this population. */
  Population<Candidate> front() {
    Population<Candidate> ret{};
    std::shared_lock<mutex_t> lock(mtx);
    for(auto& c : *this)
      if(*this << c == 0)
        ret.add(c);
    return ret;
  }


  /** \brief The return type of stat(). */
  struct Stat {
    double mean;  ///< The mean fitness of the Population.
    double stdev; ///< The standard deviation of fitness in the Population.
  };

  /** \brief Returns the mean fitness of the population and the standard
   * deviation.
   *
   * Applicable only to candidate classes whose fitness is a simple floating
   * point type or allows an implicit convertion to one. This method
   * generates an error at compile time in specializations for which this
   * condition is not satisfied. */
  Stat stat() {
    static_assert(std::is_convertible<_FitnessType, double>::value,
        "This method requires the fitness type to be convertible to double.");
    double f, sf = 0, sf2 = 0;
    std::shared_lock<mutex_t> lock(mtx);
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
  void ensureSorted() {
    static_assert(internal::comparable<_FitnessType>(0),
        "This method requires the fitness type to implement an operator<.");
    if(!sorted) {
      std::sort(begin(), end());
      sorted = true;
    }
  }
}; // class Population

} // namespace gen
