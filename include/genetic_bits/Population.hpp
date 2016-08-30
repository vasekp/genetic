namespace gen {

/** \brief The Population template.
 *
 * Requires a `Candidate` class derived from gen::Candidate. */
template<class Candidate>
class Population : private std::vector<Candidate> {
  bool sorted = false;
  mutable internal::rw_semaphore smp{};

  /* If this is a reference population of another, Candidate is actually
   * std::reference_wrapper<const _Candidate>. This removes the wrapper and
   * the const, so _Candidate refers to the template parameter of the original
   * class. Otherwise _Candidate = Candidate. This is later used to detect
   * whether we're a reference or not. It's guaranteed that Candidate& can be
   * cast to _Candidate& in any case. */
  typedef decltype(internal::unwrap(std::declval<Candidate>())) _Candidate;
  // Let that not confuse poor Doxygen.
#ifdef DOXYGEN
#define _Candidate Candidate
#endif

  typedef decltype(internal::detectFT<Candidate>(nullptr)) _FitnessType;

  typedef std::vector<Candidate> Base;

  static_assert(std::is_convertible<Candidate&, const gen::Candidate<_FitnessType>&>::value,
      "The Candidate type needs to be derived from gen::Candidate.");


public:
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::clear;
  using Base::operator[];

  /** \copybrief gen::RefPopulation
   *
   * Population::Ref can be used as short for gen::RefPopulation<Candidate> if
   * Population is defined elsewhere to be `gen::Population<Candidate>`.
   *
   * This is the return type of functions that return a subset of an existing
   * Population by reference. For convenience, objects of this type can be assigned to a
   * `Population<Candidate>`. Also, it's guaranteed that Population::Ref::Ref is
   * identical to Population::Ref, which makes it convenient to chain
   * selection function, e.g. \link randomSelect(size_t, Rng&) const `pop.randomSelect(5)`
   * \endlink`.front().randomSelect()` for
   * a simple tournament selection.
   *
   * \see RefPopulation for more discussion. */
  typedef Population<std::reference_wrapper<const _Candidate>> Ref;

  friend class Population<std::reference_wrapper<const _Candidate>>;
  friend class Population<_Candidate>;

  /** \brief Creates an empty population. */
  Population() = default;

  /** \brief Creates an empty population but preallocate space for count
   * candidates. */
  explicit Population(size_t count) {
    Base::reserve(count);
  }

  /** \brief Creates a population of size `count` whose candidates are results
   * of calls to the source function `src`.
   *
   * \see add(size_t, Source, bool, bool) for discussion about the parameters. */
  template<class Source>
  explicit Population(size_t count, Source src, bool precompute = false, bool parallel = true) {
    add(count, src, precompute, parallel);
  }

  /* The copy and move constructors: trivial but we need explicit definition
   * because the semaphore can't be default copied or moved */

  /** \brief The copy constructor. */
  Population(const Population& _p): Base(_p), sorted(_p.sorted) { }

  /** \brief The move constructor. */
  Population(Population&& _p): Base(std::move(_p)), sorted(_p.sorted) { }

  /** \brief Initializes a population as a copy of a gen::RefPopulation. */
#ifdef DOXYGEN
  explicit Population(const RefPopulation<Candidate>& _p) {
#else
  template<bool is_ref = !std::is_same<_Candidate, Candidate>::value>
  explicit Population(typename std::enable_if<!is_ref, const Ref&>::type _p) {
#endif
    Base::insert(end(), _p.begin(), _p.end());
    sorted = _p.sorted;
  }

  /** \brief Initializes a reference Population from a whole Population.
   * Exposed only for instances of gen::RefPopulation. */
#ifdef DOXYGEN
  explicit Population(const Population<Candidate>&) {
#else
  template<bool is_ref = !std::is_same<_Candidate, Candidate>::value>
  explicit Population(typename std::enable_if<is_ref, const Population<_Candidate>&>::type _p) {
#endif
    Base::insert(end(), _p.begin(), _p.end());
    sorted = _p.sorted;
  }

  /* The copy and move assignments: trivial but we need explicit definition
   * semaphore can't be default copied or moved */

  /** \brief Copy assignment operator. */
  Population& operator=(const Population& _p) {
    internal::write_lock lock(smp);
    Base::operator=(_p);
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Move assignment operator. */
  Population& operator=(Population&& _p) {
    internal::write_lock lock(smp);
    Base::operator=(std::move(_p));
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Copy assignment from a gen::RefPopulation. */
#ifdef DOXYGEN
  Population& operator=(const RefPopulation<Candidate>& _p) {
#else
  template<bool is_ref = !std::is_same<_Candidate, Candidate>::value>
  typename std::enable_if<!is_ref, Population&>::type operator=(const Ref& _p) {
#endif
    internal::write_lock lock(smp);
    Base::clear();
    Base::insert(end(), _p.begin(), _p.end());
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Adds a new candidate. */
  void add(const _Candidate& c) {
    internal::write_lock lock(smp);
    Base::push_back(c);
    sorted = false;
  }

  /** \brief Pushes back a new candidate using the move semantics. */
  void add(Candidate&& c) {
    internal::write_lock lock(smp);
    Base::push_back(c);
    sorted = false;
  }

  /** \brief Draws `count` candidates from a source function `src`.
   *
   * \param count the number of candidates to generate
   * \param src source function; one of:
   * - `std::function<(const) Candidate>`: returning by copy,
   * - `std::function<(const) Candidate&>`: returning by reference,
   * - a pointer to function returning `Candidate` or `Candidate&`,
   * - a lambda function returning either.
   * The template allows for optimizations (inlining) in the latter case.
   * \param precompute set to `true` to enable precomputing the candidates'
   * fitness (the evaluation is lazy by default)
   * \param parallel controls parallelization using OpenMP (on by default) */
  template<class Source>
  void NOINLINE add(size_t count, Source src, bool precompute = false, bool parallel = true) {
    internal::write_lock lock(smp);
    Base::reserve(size() + count);
    #pragma omp parallel if(parallel)
    {
      std::vector<Candidate> tmp{};
      tmp.reserve(count*2/omp_get_num_threads());
      #pragma omp for schedule(dynamic)
      for(size_t j = 0; j < count; j++) {
        Candidate c{src()};
        if(precompute)
          c.fitness();
        tmp.push_back(c);
      }
      #pragma omp critical
      Base::insert(end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
    }
    sorted = false;
  }

  /** \brief Copies all candidates from a vector of `Candidate`s. */
  void NOINLINE add(const std::vector<Candidate>& vec) {
    add(vec.begin(), vec.end());
  }

  /** \brief Copies an iterator range from a container of `Candidate`s. */
  template<class InputIt>
  void add(InputIt first, InputIt last) {
    internal::write_lock lock(smp);
    Base::insert(end(), first, last);
    sorted = false;
  }

  /** \brief Moves all candidates from a vector of `Candidate`s. */
  void NOINLINE add(std::vector<Candidate>&& vec) {
    internal::write_lock lock(smp);
    Base::reserve(size() + vec.size());
    Base::insert(end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));
    vec.clear();
    sorted = false;
  }

  /** \brief Copies all candidates from another Population. */
  void add(const Population<Candidate>& pop) {
    add(static_cast<const Base&>(pop));
  }

  /** \brief Moves all candidates from another Population. */
  void add(Population<Candidate>&& pop) {
    add(static_cast<Base&&>(pop));
  }

  /** \brief Copies all candidates from a gen::RefPopulation. */
#ifdef DOXYGEN
  void add(const RefPopulation<Candidate>& pop) {
#else
  template<bool is_ref = !std::is_same<_Candidate, Candidate>::value>
  typename std::enable_if<!is_ref, void>::type add(const Ref& pop) {
#endif
    add(pop.begin(), pop.end());
  }

  /** \brief Takes reference to all candidates from another Population.
   * Exposed only for instances of gen::RefPopulation. */
#ifdef DOXYGEN
  void add(const Population<Candidate>&) {
#else
  template<bool is_ref = !std::is_same<_Candidate, Candidate>::value>
  typename std::enable_if<is_ref, void>::type add(const Population<_Candidate>& pop) {
#endif
    add(pop.begin(), pop.end());
  }

  /** \brief Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates an error at compile time in
   * specializations for which this condition is not satisfied.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged. */
  void rankTrim(size_t newSize) {
    internal::write_lock lock(smp);
    if(size() <= newSize)
      return;
    ensureSorted(lock);
    const Candidate& dummy = *begin();  // needed by resize() if size() < newSize which can't happen
    Base::resize(newSize, dummy);       // (otherwise we would need a default constructor)
  }

  /** \brief Reduces the population to a maximum size given by the argument,
   * using random selection if the latter is smaller.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged.
   * \param rng the random number generator, or gen::rng by default. */
  template<class Rng = decltype(rng)>
  void randomTrim(size_t newSize, Rng& rng = rng) {
    internal::write_lock lock(smp);
    if(size() <= newSize)
      return;
    shuffle(rng);
    const Candidate& dummy = *begin();  // see rankTrim()
    Base::resize(newSize, dummy);
  }

  /** \brief Reduces the population by selective removal of candidates.
   *
   * Candidates are tested against a given predicate. If a candidate `c`
   * satisfies `pred(c)` it is removed from the population. A minimum number
   * of candidates can be set; if so, the procedure stops when enough
   * candidates have been removed to satisfy this bound.
   *
   * \param test a boolean function accepting a `const Candidate` reference.
   * If the return value is `true` the candidate is removed from the
   * population.
   * \param minSize a minimum number of candidates to be kept if possible. If
   * zero (the default value), all candidates satisfying the predicate are
   * removed.
   * \param randomize whether to randomly shuffle the sample prior to pruning
   * (this is the default). If `false` then earlier appearing candidates are
   * preferred in survival.
   * \param rng the random number generator, or gen::rng by default. Ignored
   * if `randomize` is `false`. */
  template<class Rng = decltype(rng)>
  void prune(bool (*test)(const _Candidate&), size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::write_lock lock(smp);
    if(size() <= minSize)
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++)
      if(test((*this)[i])) {
        Base::erase(begin() + i);
        if(--sz <= minSize)
          return;
      }
  }

  /** \brief Reduces the population by selective removal of candidates.
   * 
   * Candidates are tested for similarity according to a provided crierion
   * function. If a pair of candidates `(a, b)` satisfies the test, only `a`
   * is kept. A minimum number of candidates can be set; if so, the procedure
   * stops when enough candidates have been removed to satisfy this bound.
   *
   * \param test a boolean function accepting two `const Candidate`
   * references. Should be symmetric in its arguments. If the return value is
   * `true` the latter candidate is removed from the population.
   * \param minSize a minimum number of candidates to be kept if possible. If
   * zero (the default value), all duplicates are removed.
   * \param randomize whether to randomly shuffle the sample prior to pruning
   * (this is the default). If `false` then earlier appearing candidates are
   * preferred in survival.
   * \param rng the random number generator, or gen::rng by default. Ignored
   * if `randomize` is `false`. */
  template<class Rng = decltype(rng)>
  void prune(bool (*test)(const _Candidate&, const _Candidate&), size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::write_lock lock(smp);
    if(size() <= minSize)
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++)
      for(size_t j = sz - 1; j > i; j--)
        if(test((*this)[i], (*this)[j])) {
          Base::erase(begin() + j);
          if(--sz <= minSize)
            return;
        }
  }

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
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use rankSelect_v(double, Rng&) instead.
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
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const _Candidate& rankSelect(double bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<const _Candidate&>(bias, rng);
    else
      return rankSelect_two<const _Candidate&, &internal::eval_in_product<fun>>(bias, rng);
  }

  /**
   * \brief Retrieves a candidate randomly chosen by rank-based selection.
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
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use rankSelect_v(double, Rng&) instead.
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
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const _Candidate& rankSelect(double bias, Rng& rng = rng) {
    return rankSelect_two<const _Candidate&, fun>(bias, rng);
  }

  /** \copybrief rankSelect(double, Rng&)
   *
   * Works like rankSelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  _Candidate rankSelect_v(double bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<_Candidate>(bias, rng);
    else
      return rankSelect_two<_Candidate, &internal::eval_in_product<fun>>(bias, rng);
  }

  /** \copybrief rankSelect(double, Rng&)
   *
   * Works like rankSelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  _Candidate rankSelect_v(double bias, Rng& rng = rng) {
    return rankSelect_two<_Candidate, fun>(bias, rng);
  }

private:
  template<class Ret, class Rng>
  Ret NOINLINE rankSelect_exp(double bias, Rng& rng = rng) {
    static thread_local std::uniform_real_distribution<double> rDist(0, 1);
    double x = rDist(rng);
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("rankSelect(): Population is empty.");
    ensureSorted(lock);
    if(x == 1)
      return Base::back();
    else
      return (*this)[(int)(-log(1 - x + x*exp(-bias))/bias*sz)];
  }

  template<class Ret, double (*fun)(double, double), class Rng = decltype(rng)>
  Ret NOINLINE rankSelect_two(double bias, Rng& rng = rng) {
    static thread_local std::discrete_distribution<size_t> iDist{};
    static thread_local size_t last_sz = 0;
    static thread_local std::vector<double> probs{};
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("rankSelect(): Population is empty.");
    if(sz != last_sz) {
      probs.clear();
      probs.reserve(sz);
      for(size_t i = 0; i < sz; i++)
        probs.push_back(1 / fun((double)(i+1) / sz, bias));
      iDist = std::discrete_distribution<size_t>(probs.begin(), probs.end());
      last_sz = sz;
    }
    ensureSorted(lock);
    return (*this)[iDist(rng)];
  }

public:
  /** \brief Retrieves a candidate chosen using uniform random selection.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use randomSelect_v(Rng&) const instead. */
#ifdef DOXYGEN
  template<class Rng = decltype(rng)>
  const _Candidate& NOINLINE randomSelect(Rng& rng = rng) const {
#else
  template<class Rng = decltype(rng), class Ret = const _Candidate&>
  Ret NOINLINE randomSelect(Rng& rng = rng) const {
#endif
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("randomSelect(): Population is empty.");
    std::uniform_int_distribution<size_t> dist{0, sz - 1};
    return (*this)[dist(rng)];
  }

  /** \copybrief randomSelect(Rng&) const
   * \brief Works like randomSelect(Rng&) const but returns by value. */
  template<class Rng = decltype(rng)>
  _Candidate NOINLINE randomSelect_v(Rng& rng = rng) const {
    return randomSelect<Rng, _Candidate>(rng);
  }

  /** \brief Randomly selects `k` different candidates. If `k < size()`, the
   * whole population is returned.
   *
   * The returned \link RefPopulation \endlink remains valid until the
   * original population is modified.  Therefore there is a risk of
   * invalidating it in a multi-threaded program if another thread
   * concurrently modifies the population. If your code allows this, use
   * randomSelect_v(size_t, Rng&) const instead. */
#ifdef DOXYGEN
  template<class Rng = decltype(rng)>
  RefPopulation<Candidate> NOINLINE randomSelect(size_t k, Rng& rng = rng) const {
#else
  template<class Rng = decltype(rng), class Ret = Ref>
  Ret NOINLINE randomSelect(size_t k, Rng& rng = rng) const {
#endif
    internal::read_lock lock(smp);
    size_t sz = size();
    k = std::min(k, sz);
    std::vector<size_t> idx(sz);
    /* Fisher-Yates intentionally without initialization! */
    for(size_t i = 0; i < k; i++) {
      size_t d = rng() % (sz - i), j = i+d;
      std::swap(idx[i], idx[j]);
      idx[j] -= d;
      idx[i] += i + d;
    }
    idx.resize(k);
    Ret ret(k);
    for(auto i : idx)
      ret.add((*this)[i]);
    return ret;
  }

  /** \copybrief randomSelect(size_t, Rng&) const
   * \brief Works like randomSelect(size_t, Rng&) const but returns an independent
   * Population. */
  template<class Rng = decltype(rng)>
  Population<_Candidate> NOINLINE randomSelect_v(size_t k, Rng& rng = rng) const {
    return randomSelect<Rng, Population<_Candidate>>(k, rng);
  }

  /** \brief Returns the best candidate of population.
   *
   * If more candidates have equal best fitness the returned reference may be
   * any of them.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use best_v() instead.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates an error at compile time in
   * specializations for which this condition is not satisfied. */
#ifdef DOXYGEN
  const Candidate& best() {
#else
  template<class Ret = const _Candidate&>
  Ret best() {
#endif
    static_assert(internal::comparable<_FitnessType>(0),
        "This method requires the fitness type to implement an operator<.");
    internal::read_lock lock(smp);
    if(this->empty())
      throw std::out_of_range("best(): Population is empty.");
    if(sorted)
      return Base::front();
    else
      return *std::min_element(begin(), end());
  }

  /** \copybrief best()
   * \brief Works like best() but returns by value. */
  _Candidate best_v() {
    return best<_Candidate>();
  }

  /** \brief Returns the number of candidates in this population dominated by
   * a given candidate. */
  friend size_t operator<< (const _Candidate& c, const Population<Candidate>& pop) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(const _Candidate& cmp : pop)
      if(c << cmp)
        cnt++;
    return cnt;
  }

  /** \brief Returns the number of candidates in this population that
   * dominate a given candidate. */
  friend size_t operator<< (const Population<Candidate>& pop, const _Candidate& c) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(const _Candidate& cmp : pop)
      if(cmp << c)
        cnt++;
    return cnt;
  }

  /** \brief Returns a nondominated subset of this population.
   *
   * The returned \link RefPopulation \endlink remains valid until the
   * original population is modified.  Therefore there is a risk of
   * invalidating it in a multi-threaded program if another thread
   * concurrently modifies the population. If your code allows this, use
   * front_v() instead.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
#ifdef DOXYGEN
  RefPopulation<Candidate> NOINLINE front(bool parallel = true) const {
#else
  template<class Ret = Ref>
  Ret front(bool parallel = true) const {
#endif
    Ret ret{};
    internal::read_lock lock(smp);
    size_t sz = size();
    std::vector<char> dom(sz, 0);
    #pragma omp parallel for if(parallel)
    for(size_t i = 0; i < sz; i++) {
      for(size_t j = 0; j < sz; j++)
        if(!dom[j] && (*this)[j] << (*this)[i]) {
          dom[i] = 1;
          break;
        }
      if(!dom[i])
        ret.add((*this)[i]);
    }
    return ret;
  }

  /** \copybrief front()
   * \brief Works like front() but returns an independent Population. */
  Population<_Candidate> NOINLINE front_v(bool parallel = true) const {
    return front<Population<_Candidate>>(parallel);
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
   * condition is not satisfied.
   *
   * \see Stat
   */
  Stat stat() const {
    static_assert(std::is_convertible<_FitnessType, double>::value,
        "This method requires the fitness type to be convertible to double.");
    if(this->empty())
      throw std::out_of_range("stat(): Population is empty.");
    double f, sf = 0, sf2 = 0;
    internal::read_lock lock(smp);
    for(const _Candidate &c : *this) {
      f = c.fitness();
      sf += f;
      sf2 += f*f;
    }
    size_t sz = size();
    double dev2 = sf2/sz - sf/sz*sf/sz;
    return Stat{sf/sz, dev2 >= 0 ? sqrt(dev2) : 0};
  }

  private:
  void ensureSorted(internal::rw_lock& lock) {
    static_assert(internal::comparable<_FitnessType>(0),
        "This method requires the fitness type to implement an operator<.");
    if(!sorted) {
      internal::upgrade_lock up(lock);
      std::sort(begin(), end());
      sorted = true;
    }
  }

  template<class Rng = decltype(rng)>
  void shuffle(Rng& rng) {
    std::shuffle(begin(), end(), rng);
    sorted = false;
  }
}; // class Population


/** \brief A "reference population", a helper type for functions returning a
 * selection from an existing Population.
 *
 * This effectively allows to make a Population of references to Candidate,
 * e.g., to store the return value of Population::randomSelect(size_t, Rng&) const.
 * `RefPopulation<Candidate>` retains the capabilities of `Population<Candidate>`
 * like functions dependent on the existence of Candidate::operator< or
 * Candidate::operator<<. Most functions can be called on a `RefPopulation` and
 * it can be converted from and to normal population (taking references and
 * creating copies of all elements, respectively).
 *
 * Note that the elements of a \link RefPopulation \endlink are in fact
 * [`std::reference_wrapper`](http://en.cppreference.com/w/cpp/utility/functional/
 * reference_wrapper)`<const Candidate>`. This can implicitly be converted
 * to `const Candidate&` but some other uses may need to be modified. Most
 * importantly, the dot operator can not be overloaded within C++11 so
 * `x.fitness()` won't work for `x` taken from a \link RefPopulation \endlink.
 * If this is needed, convert explicitly to `const Candidate&` or use
 * [`std::reference_wrapper::get()`](http://en.cppreference.com/w/cpp/utility/
 * functional/reference_wrapper/get).
 *
 * \see Population::Ref */
template<class Candidate>
using RefPopulation = typename Population<Candidate>::Ref;

} // namespace gen
