namespace gen {

/** \brief The Population template.
 *
 * \param Candidate the class describing individual members of this
 * population. Must be derived from gen::Candidate.
 * \param Tag an optional supplement class to accompany each candidate. Used
 * for internal purposes. This does not enter iterations over this population
 * and is not duplicated when copies or references are taken.
 * \param ref if set to `true`, this is a reference population. See
 * Population::Ref for more details. */
template<class Candidate, class Tag = internal::empty, bool ref = false>
class Population : private std::vector<internal::Tagged<Candidate, Tag, ref>> {
  bool sorted = false;
  mutable internal::rw_semaphore smp{};

  typedef decltype(internal::detectFT<Candidate>(nullptr)) _FitnessType;

  typedef std::vector<internal::Tagged<Candidate, Tag, ref>> Base;

  typedef internal::cast_iterator<const Candidate, typename Base::const_iterator> CastIter;

  static_assert(std::is_convertible<Candidate&, const gen::Candidate<_FitnessType>&>::value,
      "The Candidate type needs to be derived from gen::Candidate.");

public:
  using Base::size;
  using Base::clear;
  using Base::operator[];

#ifndef DOXYGEN
  CastIter begin() const {
    return CastIter(Base::cbegin());
  }

  CastIter end() const {
    return CastIter(Base::cend());
  }
#endif

  /** \brief A corresponding "reference population", a helper type for
   * functions returning a selection from an existing Population.
   *
   * This is the return type of functions that return a subset of an existing
   * Population by reference as it holds its members by reference rather than
   * by value. Internally it is a different Population specialization with a
   * matching Candidate type, so the same functions including adding or
   * erasing can be called on it despite its temporary nature and assignments
   * between Population and Population::Ref objects are allowed and behave
   * as expected, see add(Container&).
   *
   * It is guaranteed that Population::Ref::Ref is identical to
   * Population::Ref, which makes it convenient to chain selection function,
   * e.g. \link randomSelect(size_t, Rng&) const `pop.randomSelect(5)`\endlink
   * `.front().randomSelect()` for a simple tournament selection.
   *
   * It is the user's responsibility not to use the references beyond their
   * scope. They are invalidated by operations which modify the original
   * Population, namely all operations adding, removing, and reordering its
   * elements.  This includes operations which need to sort the contents
   * first, like rankSelect(). */
  typedef Population<Candidate, Tag, true> Ref;

  /* Befriend all compatible Populations */
  template<class Tag_, bool ref_>
  friend class Population<Candidate, Tag, ref>;

  /** \brief Creates an empty population. */
  Population() = default;

  /** \brief Creates an empty population but preallocate space for count
   * candidates. */
  explicit Population(size_t count) {
    Base::reserve(count);
  }

  /** \brief Creates a population of size `count` whose candidates are results
   * of calls to the source function `src`.
   * \copydetails add(size_t, Source, bool, bool) */
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

  /** \brief Initializes this population from a compatible Population.
   * \copydetails add(Container&) */
  template<class Tag_, bool ref_>
  Population(const Population<Candidate, Tag_, ref_>& _p) {
    add(_p);
    sorted = _p.sorted;
  }

  /** \brief Initializes this population from a compatible Population using
   * move semantics. */
  template<class Tag_, bool ref_>
  Population(Population<Candidate, Tag_, ref_>&& _p) {
    add(_p);
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

  /** \brief Copy assignment of a compatible Population.
   * \copydetails add(Container&) */
  template<class Tag_, bool ref_>
  Population& operator=(const Population<Candidate, Tag_, ref_>& _p) {
    internal::write_lock lock(smp);
    Base::clear();
    Base::insert(Base::end(), _p.begin(), _p.end());
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Move assignment of a compatible Population. */
  template<class Tag_, bool ref_>
  Population& operator=(Population<Candidate, Tag_, ref_>&& _p) {
    internal::write_lock lock(smp);
    Base::clear();
    Base::insert(Base::end(),
        internal::move_iterator<decltype(_p.begin())>(_p.begin()),
        internal::move_iterator<decltype(_p.end())>(_p.end()));
    sorted = _p.sorted;
    _p.clear();
    return *this;
  }

  /** \brief Adds a new candidate. */
  void add(const Candidate& c) {
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
      Base::insert(Base::end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
    }
    sorted = false;
  }

  /** \brief Copies an iterator range from a container of `Candidate`s. */
  template<class InputIt>
  void add(InputIt first, InputIt last) {
    internal::write_lock lock(smp);
    Base::insert(Base::end(), first, last);
    sorted = false;
  }

  /** \brief Copies all candidates from a container of `Candidate`s.
   *
   * If this population is a reference population, references to all members
   * of the argument are taken, copies are made otherwise. */
  template<class Container>
  void NOINLINE add(Container& vec) {
    add(vec.begin(), vec.end());
  }

#ifndef DOXYGEN
  /** \brief Copies all candidates from a container of `Candidate`s. */
  template<class Container>
  void NOINLINE add(const Container& vec) {
    add(vec.begin(), vec.end());
  }
#endif

  /** \brief Moves all candidates from a container of `Candidate`s. */
  template<class Container>
  void NOINLINE add(Container&& vec) {
    add(internal::move_iterator<decltype(vec.begin())>(vec.begin()),
        internal::move_iterator<decltype(vec.end())>(vec.end()));
    vec.clear();
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
  void prune(bool (*test)(const Candidate&), size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::write_lock lock(smp);
    if(size() <= minSize)
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++)
      if(test((*this)[i])) {
        Base::erase(Base::begin() + i);
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
  void prune(bool (*test)(const Candidate&, const Candidate&), size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::write_lock lock(smp);
    if(size() <= minSize)
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++)
      for(size_t j = sz - 1; j > i; j--)
        if(test((*this)[i], (*this)[j])) {
          Base::erase(Base::begin() + j);
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
  const Candidate& rankSelect(double bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<const Candidate&>(bias, rng);
    else
      return rankSelect_two<const Candidate&, &internal::eval_in_product<fun>>(bias, rng);
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
  const Candidate& rankSelect(double bias, Rng& rng = rng) {
    return rankSelect_two<const Candidate&, fun>(bias, rng);
  }

  /** \copybrief rankSelect(double, Rng&)
   *
   * Works like rankSelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate rankSelect_v(double bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<Candidate>(bias, rng);
    else
      return rankSelect_two<Candidate, &internal::eval_in_product<fun>>(bias, rng);
  }

  /** \copybrief rankSelect(double, Rng&)
   *
   * Works like rankSelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate rankSelect_v(double bias, Rng& rng = rng) {
    return rankSelect_two<Candidate, fun>(bias, rng);
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
  const Candidate& randomSelect(Rng& rng = rng) const {
#else
  template<class Rng = decltype(rng), class Ret = const Candidate&>
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
  Candidate NOINLINE randomSelect_v(Rng& rng = rng) const {
    return randomSelect<Rng, Candidate>(rng);
  }

  /** \brief Randomly selects `k` different candidates. If `k < size()`, the
   * whole population is returned.
   *
   * The returned Ref remains valid until the original population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use randomSelect_v(size_t, Rng&) const instead. */
#ifdef DOXYGEN
  template<class Rng = decltype(rng)>
  Ref randomSelect(size_t k, Rng& rng = rng) const {
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
      ret.add((const Candidate&)(*this)[i]);
    return ret;
  }

  /** \copybrief randomSelect(size_t, Rng&) const
   * \brief Works like randomSelect(size_t, Rng&) const but returns an independent
   * Population. */
  template<class Rng = decltype(rng)>
  Population NOINLINE randomSelect_v(size_t k, Rng& rng = rng) const {
    return randomSelect<Rng, Population>(k, rng);
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
  template<class Ret = const Candidate&>
  Ret NOINLINE best() {
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
  Candidate best_v() {
    return best<Candidate>();
  }

  /** \brief Returns the number of candidates in this population dominated by
   * a given candidate. */
  friend size_t operator<< (const Candidate& c, const Population& pop) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(const Candidate& cmp : pop)
      if(c << cmp)
        cnt++;
    return cnt;
  }

  /** \brief Returns the number of candidates in this population that
   * dominate a given candidate. */
  friend size_t operator<< (const Population& pop, const Candidate& c) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(const Candidate& cmp : pop)
      if(cmp << c)
        cnt++;
    return cnt;
  }

  /** \brief Returns a nondominated subset of this population.
   *
   * The returned Ref remains valid until the original population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use front_v() instead.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
#ifdef DOXYGEN
  Ref front(bool parallel = true) const {
#else
  template<class Ret = Ref>
  Ret NOINLINE front(bool parallel = true) const {
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
        ret.add((const Candidate&)(*this)[i]);
    }
    return ret;
  }

  /** \copybrief front()
   * \brief Works like front() but returns an independent Population. */
  Population NOINLINE front_v(bool parallel = true) const {
    return front<Population>(parallel);
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
    for(const Candidate &c : *this) {
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
      std::sort(Base::begin(), Base::end());
      sorted = true;
    }
  }

  template<class Rng = decltype(rng)>
  void shuffle(Rng& rng) {
    std::shuffle(Base::begin(), Base::end(), rng);
    sorted = false;
  }
}; // class Population

} // namespace gen
