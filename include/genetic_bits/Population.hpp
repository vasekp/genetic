namespace gen {

namespace internal {

template<class CBase, class Tag, bool is_ref>
using PBase = std::vector<internal::CandidateTagged<CBase, Tag, is_ref>>;

} // namespace internal


/** \brief The Population template.
 *
 * From the following template parameters, only `CBase` is mandatory. For most
 * purposes only `Population<CBase>` should be needed externally. The other
 * parameters are used mainly for return types of temporary objects. For any
 * values of `Tag` and `is_ref`, `Population<CBase, Tag, is_ref>` can be
 * implicitly converted to a `Population<CBase>`.
 *
 * Population can be used as a container of `Candidate`s with read-only
 * access. (See Candidate for discussion about the relation between a
 * Candidate and `CBase`.) The functions begin() and end() are exposed,
 * returning random access iterators dereferencable to `const
 * Candidate<CBase>&` and allowing the iteration patterns
 * ```
 * for(auto& c : pop) { ... }
 * ```
 * and
 * ```
 * for(auto c : pop) { ... }
 * ```
 * Also, read-only element accessors at() and operator[]() are available for
 * direct access to candidates.
 *
 * \tparam CBase the base class of the member candidates of this population.
 * See Candidate for details.
 * \tparam Tag an optional supplement class or literal type to accompany each
 * candidate. Used for internal purposes. This does not enter iterations over
 * this population and is not duplicated when copies or references are taken.
 * When supplied, must represent a default-constructible type.
 * \tparam is_ref if set to `true`, this is a reference population. See
 * Population::Ref for more details. */
template<class CBase, class Tag = internal::empty, bool is_ref = false>
class Population : private internal::PBase<CBase, Tag, is_ref> {
  bool sorted = false;
  mutable internal::rw_semaphore smp{};

  typedef internal::PBase<CBase, Tag, is_ref> Base;

  typedef internal::cast_iterator<const Candidate<CBase>&, typename Base::iterator> iterator;
  typedef internal::cast_iterator<const Candidate<CBase>&, typename Base::const_iterator> const_iterator;

public:
  /** \brief A corresponding "reference population", a helper type for
   * functions returning a selection from an existing Population.
   *
   * This is the return type of functions that return a subset of an existing
   * Population by reference as it holds its members by reference rather than
   * by value. Internally it is a different Population specialization with a
   * matching CBase type, so the same functions including adding or erasing
   * can be called on it despite its temporary nature and assignments between
   * Population and Population::Ref objects are allowed and behave as
   * expected, see add(Container&).
   *
   * It is guaranteed that Population::Ref::Ref is identical to
   * Population::Ref, which makes it convenient to chain selection functions,
   * e.g.
   * ```
   * pop.randomSelect(5).front().randomSelect()
   * ```
   * for a simple multi-objective tournament selection.
   *
   * It is the user's responsibility not to use the references beyond their
   * scope. They are invalidated by operations which modify the original
   * Population, namely all operations adding, removing, and reordering its
   * elements.  This includes rankSelect(), which needs to sort the contents
   * for its operation, and reserve(), which may move the contents to a new
   * memory location. */
  typedef Population<CBase, Tag, true> Ref;

   /* Befriend all compatible Populations */
   template<class, class, bool>
   friend class Population;

  /** \brief Creates an empty population. */
  Population() = default;

  /* The copy and move constructors: trivial but we need explicit definition
   * because the semaphore can't be default copied or moved */

  /** \brief The copy constructor. */
  Population(const Population& _p): Base(_p), sorted(_p.sorted) { }

  /** \brief The move constructor. */
  Population(Population&& _p) noexcept: Base(std::move(_p)), sorted(_p.sorted) { }

  /** \brief Creates an empty population but preallocate space for `count`
   * candidates. */
  explicit Population(size_t count) {
    Base::reserve(count);
  }

  /** \brief Creates a population of size `count` whose candidates are results
   * of calls to the source function `src`.
   * \copydetails add(size_t, Source, bool) */
  template<class Source>
  explicit Population(size_t count, Source src, bool parallel = true) {
    add(count, src, parallel);
  }

  /** \brief Initializes this population from an iterator range from a
   * container of `Candidate`s or `CBase`s. */
  template<class InputIt>
  explicit Population(InputIt first, InputIt last) {
    add(first, last);
  }

  /** \brief Initializes this population from a container of `Candidate`s or
   * `CBase`s (e.g., a `std::vector` or another Population).
   * \copydetails add(const Container&) */
  template<class Container>
  explicit Population(const Container& vec) {
    add(std::forward<Container>(vec));
  }

  /** \brief Initializes this population from a container of `Candidate`s or
   * `CBase`s using move semantics, leaving the original container empty. */
  template<class Container>
  explicit Population(Container&& vec) {
    add(std::forward<Container>(vec));
  }

#ifndef DOXYGEN
  /* There's no added functionality w.r.t. the above as far as the user is
   * concerned, no need to document */
  template<class Tag_, bool ref_>
  Population(const Population<CBase, Tag_, ref_>& _p) {
    add(_p);
    sorted = _p.sorted;
  }

  template<class Tag_, bool ref_>
  Population(Population<CBase, Tag_, ref_>&& _p) {
    add(std::move(_p));
    sorted = _p.sorted;
  }
#endif

  /** \brief Copy assignment operator. */
  Population& operator=(const Population& _p) {
    internal::write_lock lock(smp);
    Base::operator=(_p);
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Move assignment operator. */
  Population& operator=(Population&& _p) noexcept {
    internal::write_lock lock(smp);
    Base::operator=(std::move(_p));
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Copy assignment of a compatible Population.
   * \copydetails add(const Container&) */
  template<class Tag_, bool ref_>
  Population& operator=(const Population<CBase, Tag_, ref_>& _p) {
    internal::write_lock lock(smp);
    Base::clear();
    Base::insert(Base::end(), _p.begin(), _p.end());
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Move assignment of a compatible Population. */
  template<class Tag_, bool ref_>
  Population& operator=(Population<CBase, Tag_, ref_>&& _p) {
    internal::write_lock lock(smp);
    Base::clear();
    move_add_unguarded(_p);
    sorted = _p.sorted;
    return *this;
  }

#ifdef DOXYGEN
  /** \brief Returns the current count of candidates. */
  size_t size();
#else
  using Base::size;
#endif

  /** \brief Empties the population. */
  void clear() {
    internal::write_lock lock(smp);
    Base::clear();
  }

  /** \brief Reserves space for `count` candidates.
   *
   * If `count` is larger than the actual size of the population, all
   * references may be invalidated. */
  void reserve(size_t count) {
    internal::write_lock lock(smp);
    Base::reserve(count);
  }

  /** \brief Returns an iterator to the beginning. */
  iterator begin() {
    return iterator(Base::begin());
  }

  /** \brief Returns an iterator to the end. */
  iterator end() {
    return iterator(Base::end());
  }

  /** \brief Returns a constant iterator to the beginning. */
  const_iterator begin() const {
    return const_iterator(Base::begin());
  }

  /** \brief Returns a constant iterator to the end. */
  const_iterator end() const {
    return const_iterator(Base::end());
  }

  /** \brief Read-only access to a specified element (no bounds checking). */
  const Candidate<CBase>& operator[](size_t pos) const {
    return static_cast<const Candidate<CBase>&>(Base::operator[](pos));
  }

  /** \brief Read-only access to a specified element (with bounds checking). */
  const Candidate<CBase>& at(size_t pos) const {
    return static_cast<const Candidate<CBase>&>(Base::at(pos));
  }

  /** \brief Adds a new candidate. */
  void add(const Candidate<CBase>& c) {
    internal::write_lock lock(smp);
    Base::push_back(c);
    sorted = false;
  }

  /** \brief Pushes back a new candidate using the move semantics. */
  void add(Candidate<CBase>&& c) {
    internal::write_lock lock(smp);
    Base::push_back(std::move(c));
    sorted = false;
  }

  /** \brief Draws `count` candidates from a source function `src`.
   *
   * \param count the number of candidates to generate
   * \param src source function; can be any callable object (e.g., a
   * `std::function`, a function pointer, or a lambda function) returning
   * either of `Candidate<CBase>` or `CBase` and either by value or by
   * reference (a copy will be taken). In many cases the function
   * call can be inlined by the optimizer if known at compile time.
   * \param parallel controls parallelization using OpenMP (on by default) */
  template<class Source>
  void NOINLINE add(size_t count, Source src, bool parallel = true) {
    internal::write_lock lock(smp);
    Base::reserve(size() + count);
    #pragma omp parallel if(parallel)
    {
      std::vector<Candidate<CBase>> tmp{};
      tmp.reserve(count*2/omp_get_num_threads());
      #pragma omp for schedule(dynamic)
      for(size_t j = 0; j < count; j++) {
        Candidate<CBase> c{src()};
        tmp.push_back(c);
      }
      #pragma omp critical
      Base::insert(Base::end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
    }
    sorted = false;
  }

  /** \brief Copies an iterator range from a container of candidates. */
  template<class InputIt>
  void add(InputIt first, InputIt last) {
    assert_iterator<InputIt>();
    internal::write_lock lock(smp);
    Base::insert(Base::end(), first, last);
    sorted = false;
  }

  /** \brief Copies all candidates from a container of candidates.
   *
   * If this population is a reference population, references to all members
   * of the argument are taken, copies are made otherwise. */
  template<class Container>
#ifdef DOXYGEN
  void
#else
  typename std::enable_if<internal::is_container<Container>(0), void>::type
#endif
  NOINLINE add(const Container& vec) {
    add(vec.begin(), vec.end());
  }

  /** \brief Moves all candidates from a container of candidates.
   *
   * Moves between `Population`s are only supported if source and destination
   * are both non-reference or both reference and are of the same `CBase` type. */
  template<class Container>
#ifdef DOXYGEN
  void
#else
  typename std::enable_if<
      internal::is_container<Container>(0) &&
      std::is_rvalue_reference<Container&&>::value,
    void>::type
#endif
  NOINLINE add(Container&& vec) {
    internal::write_lock lock(smp);
    move_add_unguarded(std::forward<Container>(vec));
    sorted = false;
  }

private:
  template<class Container>
  void NOINLINE move_add_unguarded(Container&& vec) {
    assert_iterator<decltype(vec.begin())>();
    Base::insert(Base::end(),
        internal::move_iterator<decltype(vec.begin())>(vec.begin()),
        internal::move_iterator<decltype(vec.end())>(vec.end()));
    vec.clear();
  }

public:
  /** \brief Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample.
   *
   * Applicable only if the type returned by `CBase::fitness()` allows total
   * ordering using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged. */
  void rankTrim(size_t newSize) {
    internal::read_lock lock(smp);
    if(!lock.upgrade_if([newSize,this]() -> bool { return size() > newSize; }))
      return;
    ensureSorted(lock);
    auto& dummy = Base::front();  // needed by resize() if size() < newSize which can't happen
    Base::resize(newSize, dummy); // (otherwise we would need a default constructor)
  }

  /** \brief Reduces the population to a maximum size given by the argument,
   * using random selection if the latter is smaller.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged.
   * \param rng the random number generator, or gen::rng by default. */
  template<class Rng = decltype(rng)>
  void randomTrim(size_t newSize, Rng& rng = rng) {
    internal::read_lock lock(smp);
    if(!lock.upgrade_if([newSize,this]() -> bool { return size() > newSize; }))
      return;
    shuffle(rng);
    auto& dummy = Base::front(); // see rankTrim()
    Base::resize(newSize, dummy);
  }

  /** \brief Reduces the population by selective removal of candidates.
   *
   * Candidates are tested against a given predicate. If a candidate `c`
   * satisfies `pred(c)` it is removed from the population. A minimum number
   * of candidates can be set; if so, the procedure stops when enough
   * candidates have been removed to satisfy this bound.
   *
   * \param test a boolean function accepting a constant Candidate<CBase>
   * reference. If the return value is `true` the candidate is removed from
   * the population.
   * \param minSize a minimum number of candidates to be kept if possible. If
   * zero (the default value), all candidates satisfying the predicate are
   * removed.
   * \param randomize whether to randomly shuffle the sample prior to pruning
   * (this is the default). If `false` then earlier appearing candidates are
   * preferred in survival.
   * \param rng the random number generator, or gen::rng by default. Ignored
   * if `randomize` is `false`. */
  template<class Rng = decltype(rng)>
  void prune(bool (*test)(const Candidate<CBase>&), size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::read_lock lock(smp);
    if(!lock.upgrade_if([minSize,this]() -> bool { return size() > minSize; }))
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++)
      if(test(operator[](i))) {
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
   * \param test a boolean function accepting two constant Candidate<CBase>
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
  void prune(bool (*test)(const Candidate<CBase>&, const Candidate<CBase>&), size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::read_lock lock(smp);
    if(!lock.upgrade_if([minSize,this]() -> bool { return size() > minSize; }))
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++)
      for(size_t j = sz - 1; j > i; j--)
        if(test(operator[](i), operator[](j))) {
          Base::erase(Base::begin() + j);
          if(--sz <= minSize)
            return;
        }
  }

#ifdef DOXYGEN

  /** \brief Retrieves a candidate randomly chosen by rank-based selection.
   *
   * In determining the probability of each candidate, its rank, rescaled to
   * the range (0, 1], is processed by a given function as described below,
   * and the returned value is interpreted as an inverse probability.
   *
   * The population is sorted from the lowest to the highest fitness value and
   * each candidate is assigned a rank between 1 and `max`. If two or more
   * candidates have an equal fitness, they will be assigned different ranks
   * in an undefined order. This rank is then divided by `max` before passing
   * it to the processing function.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use rankSelect_v() instead.
   *
   * Applicable only if the type returned by `CBase::fitness()` allows total
   * ordering using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * \tparam fun A `constexpr` pointer to a function of signature
   * either `double(*)(double)` or `double(*)(double, double)`. In the former
   * case, the argument is `x * bias`, in the latter case the arguments
   * are `x` and `bias`, where `x` denotes the rescaled rank of each candidate.
   * It must be positive and strictly increasing in `x` for `bias > 0`.
   * This function will be built in at compile time, eliminating a function
   * pointer lookup. The default is `std::exp`, for which an fast specialized
   * algorithm is provided, another usual choice is `std::pow`.
   * \param bias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& rankSelect(double bias, Rng& rng = rng);

  /** \copybrief rankSelect()
   *
   * Works like rankSelect() but returns by value.
   *
   * \returns a copy of the randomly chosen candidate. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> rankSelect_v(double bias, Rng& rng = rng);

#else

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& rankSelect(double max, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<const Candidate<CBase>&>(max, rng);
    else
      return rankSelect_two<const Candidate<CBase>&, &internal::eval_in_product<fun>>(max, rng);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate<CBase>& rankSelect(double bias, Rng& rng = rng) {
    return rankSelect_two<const Candidate<CBase>&, fun>(bias, rng);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> rankSelect_v(double bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<Candidate<CBase>>(bias, rng);
    else
      return rankSelect_two<Candidate<CBase>, &internal::eval_in_product<fun>>(bias, rng);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate<CBase> rankSelect_v(double bias, Rng& rng = rng) {
    return rankSelect_two<Candidate<CBase>, fun>(bias, rng);
  }

#endif

private:
  std::uniform_real_distribution<double> _uniform{0, 1};

  template<class Ret, class Rng>
  Ret NOINLINE rankSelect_exp(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    double x = _uniform(rng);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("rankSelect(): Population is empty.");
    ensureSorted(lock);
    if(x == 1)
      return Base::back();
    else
      return operator[]((int)(-log(1 - x + x*exp(-bias))/bias*sz));
  }

  std::discrete_distribution<size_t> _rankSelect_dist{};
  std::vector<double> _rankSelect_probs{};
  size_t _rankSelect_last_sz{};
  double _rankSelect_last_bias{};

  template<class Ret, double (*fun)(double, double), class Rng>
  Ret NOINLINE rankSelect_two(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("rankSelect(): Population is empty.");
    if(sz != _rankSelect_last_sz || bias != _rankSelect_last_bias) {
      lock.upgrade();
      _rankSelect_probs.clear();
      _rankSelect_probs.reserve(sz);
      for(size_t i = 0; i < sz; i++)
        _rankSelect_probs.push_back(1 / fun((double)(i+1) / sz, bias));
      _rankSelect_dist = std::discrete_distribution<size_t>(_rankSelect_probs.begin(), _rankSelect_probs.end());
      _rankSelect_last_sz = sz;
      _rankSelect_last_bias = bias;
    }
    ensureSorted(lock);
    return operator[](_rankSelect_dist(rng));
  }

public:

#ifdef DOXYGEN

  /** \brief Retrieves a candidate randomly chosen by fitness-based selection.
   *
   * In determining the probability of each candidate, its fitness is
   * processed by a given function as described below, and the returned value
   * is interpreted as an inverse probability.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use fitnessSelect_v() instead.
   *
   * Applicable only if the type returned by `CBase::fitness()` can be
   * converted to `double`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * \tparam fun A `constexpr` pointer to a function of signature
   * either `double(*)(double)` or `double(*)(double, double)`. In the former
   * case, the argument is `x * bias`, in the latter case the arguments
   * are `x` and `bias`, where `x` denotes the fitness of each candidate.
   * It must be positive and strictly increasing in `x` for `bias > 0`.
   * This function will be built in at compile time, eliminating a function
   * pointer lookup. The default is `std::exp`, another usual choice is
   * `std::pow`.
   * \param bias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& fitnessSelect(double bias, Rng& rng = rng);

  /** \copybrief fitnessSelect()
   *
   * Works like fitnessSelect() but returns by value.
   *
   * \returns a copy of the randomly chosen candidate. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> fitnessSelect_v(double bias, Rng& rng = rng);

#else

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& fitnessSelect(double mult, Rng& rng = rng) {
    return fitnessSelect_int<const Candidate<CBase>&, &internal::eval_in_product<fun>>(mult, rng);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate<CBase>& fitnessSelect(double bias, Rng& rng = rng) {
    return fitnessSelect_int<const Candidate<CBase>&, fun>(bias, rng);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> fitnessSelect_v(double bias, Rng& rng = rng) {
    return fitnessSelect_int<Candidate<CBase>, &internal::eval_in_product<fun>>(bias, rng);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate<CBase> fitnessSelect_v(double bias, Rng& rng = rng) {
    return fitnessSelect_int<Candidate<CBase>, fun>(bias, rng);
  }

#endif

private:
  std::discrete_distribution<size_t> _fitnessSelect_dist{};
  std::vector<double> _fitnessSelect_probs{};
  size_t _fitnessSelect_last_mod{(size_t)-1};
  double _fitnessSelect_last_bias{};

  template<class Ret, double (*fun)(double, double), class Rng>
  Ret NOINLINE fitnessSelect_int(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    size_t sz = size();
    assert_double();
    if(sz == 0)
      throw std::out_of_range("fitnessSelect(): Population is empty.");
    if(_fitnessSelect_last_mod != smp.get_mod_cnt() || bias != _fitnessSelect_last_bias) {
      lock.upgrade();
      _fitnessSelect_probs.clear();
      _fitnessSelect_probs.reserve(sz);
      for(auto& c : *this)
        _fitnessSelect_probs.push_back(1 / fun(c.fitness(), bias));
      _fitnessSelect_dist = std::discrete_distribution<size_t>(_fitnessSelect_probs.begin(), _fitnessSelect_probs.end());
      _fitnessSelect_last_mod = smp.get_mod_cnt();
      _fitnessSelect_last_bias = bias;
    }
    return operator[](_fitnessSelect_dist(rng));
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
  const Candidate<CBase>& randomSelect(Rng& rng = rng) const {
#else
  template<class Rng = decltype(rng), class Ret = const Candidate<CBase>&>
  Ret NOINLINE randomSelect(Rng& rng = rng) const {
#endif
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("randomSelect(): Population is empty.");
    std::uniform_int_distribution<size_t> dist{0, sz - 1};
    return operator[](dist(rng));
  }

  /** \copybrief randomSelect(Rng&) const
   *
   * Works like randomSelect(Rng&) const but returns by value. */
  template<class Rng = decltype(rng)>
  Candidate<CBase> NOINLINE randomSelect_v(Rng& rng = rng) const {
    return randomSelect<Rng, Candidate<CBase>>(rng);
  }

  /** \brief Randomly selects `k` different candidates. If `k ≥ size()`, the
   * entire population is returned.
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
    if(k >= sz)
      return Ret(*this);
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
      ret.add(operator[](i));
    return ret;
  }

  /** \copybrief randomSelect(size_t, Rng&) const
   *
   * Works like randomSelect(size_t, Rng&) const but returns an independent
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
   * Applicable only if the type returned by `CBase::fitness()` allows total
   * ordering using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied. */
#ifdef DOXYGEN
  const Candidate<CBase>& best() {
#else
  template<class Ret = const Candidate<CBase>&>
  Ret NOINLINE best() {
#endif
    assert_comparable();
    internal::read_lock lock(smp);
    if(this->empty())
      throw std::out_of_range("best(): Population is empty.");
    if(sorted)
      return Base::front();
    else
      return *std::min_element(begin(), end());
  }

  /** \copybrief best()
   *
   * Works like best() but returns by value. */
  Candidate<CBase> best_v() {
    return best<Candidate<CBase>>();
  }

  /** \brief Returns the number of candidates in this population dominated by
   * a given candidate. */
  friend size_t operator<< (const Candidate<CBase>& c, const Population& pop) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(auto& cmp : pop)
      if(c << cmp)
        cnt++;
    return cnt;
  }

  /** \brief Returns the number of candidates in this population that
   * dominate a given candidate. */
  friend size_t operator<< (const Population& pop, const Candidate<CBase>& c) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(auto& cmp : pop)
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
    /* If we're a NSGA population, the front may have been precomputed */
    if(front_nsga(ret))
      return ret;
    size_t sz = size();
    std::vector<char> dom(sz, 0);
    #pragma omp parallel for if(parallel)
    for(size_t i = 0; i < sz; i++) {
      for(size_t j = 0; j < sz; j++)
        if(!dom[j] && operator[](j) << operator[](i)) {
          dom[i] = 1;
          break;
        }
      if(!dom[i])
        ret.add(operator[](i)); // thread-safe
    }
    return ret;
  }

  /** \copybrief front()
   *
   * Works like front() but returns an independent Population. */
  Population NOINLINE front_v(bool parallel = true) const {
    return front<Population>(parallel);
  }

private:
  template<class Ret, typename T = Tag>
  typename std::enable_if<std::is_same<T, size_t>::value, bool>::type
  front_nsga(Ret& ret) const {
    if(smp.get_mod_cnt() == _nsgaSelect_last_mod) {
      for(auto& tg : (Base&)(*this))
        if(tg.tag() == 0)
          ret.add(static_cast<Candidate<CBase>&>(tg));
      return true;
    } else {
      return false;
    }
  }

  template<class Ret, typename T = Tag>
  typename std::enable_if<!std::is_same<T, size_t>::value, bool>::type
  front_nsga(Ret&) const {
    return false;
  }

public:

#ifdef DOXYGEN

  /** \brief Retrieves a candidate randomly chosen by the NSGA algorithm.
   *
   * In determining the probability of each candidate, the length of the
   * longest chain of candidates successively dominating it is identified and
   * processed by a given function as described below, whose returned value
   * is interpreted as an inverse probability.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use NSGASelect_v() instead.
   *
   * Applicable only if the `Tag` parameter of this Population is `size_t`.
   * This method generates a compile-time error in specializations for
   * which this condition is not satisfied.
   *
   * \tparam fun A `constexpr` pointer to a function of signature
   * either `double(*)(double)` or `double(*)(double, double)`. In the former
   * case, the argument is `x * bias`, in the latter case the arguments
   * are `x` and `bias`, where `x` denotes the dominance rank of the candidate.
   * It must be positive and strictly increasing in `x` for `bias > 0`.
   * This function will be built in at compile time, eliminating a function
   * pointer lookup. The default is `std::exp`, another usual choice is
   * `std::pow`.
   * \param bias > 0 determines how much low-dominated solutions are preferred.
   * Zero would mean no account on dominance rank in the selection process
   * whatsoever. The bigger the value the more low-dominated candidates are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& NSGASelect(double bias, Rng& rng = rng);

  /** \copybrief NSGASelect()
   *
   * Works like NSGASelect() but returns by value.
   *
   * \returns a copy of the randomly chosen candidate. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> NSGASelect_v(double bias, Rng& rng = rng);

#else

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& NSGASelect(double mult, Rng& rng = rng) {
    return NSGASelect_int<const Candidate<CBase>&, &internal::eval_in_product<fun>>(mult, rng);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate<CBase>& NSGASelect(double bias, Rng& rng = rng) {
    return NSGASelect_int<const Candidate<CBase>&, fun>(bias, rng);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> NSGASelect_v(double bias, Rng& rng = rng) {
    return NSGASelect_int<Candidate<CBase>, &internal::eval_in_product<fun>>(bias, rng);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate<CBase> NSGASelect_v(double bias, Rng& rng = rng) {
    return NSGASelect_int<Candidate<CBase>, fun>(bias, rng);
  }

#endif

private:
  struct _nsga_struct {
    const Candidate<CBase>& r; // reference to a candidate
    size_t& rank;       // reference to its Tag in the original population
    size_t dom_cnt;     // number of candidates dominating this
    std::deque<_nsga_struct*> q;  // candidates dominated by this
  };

  /* The rating needs to be done with the particular Population's lock,
   * otherwise NSGASelect would think it's newer than the last record and
   * would recalculate it. This function helps in NSGAref() where only the
   * source Population's lock is used. */
  void _nsga_rate(bool parallel = false) {
    internal::write_lock lock(smp);
    _nsga_rate(lock, parallel);
    _nsgaSelect_last_mod = smp.get_mod_cnt();
  }

  void NOINLINE _nsga_rate(internal::rw_lock& lock, bool parallel = false) {
    internal::upgrade_lock up(lock);
    std::list<_nsga_struct> ref{};
    for(auto& tg : static_cast<Base&>(*this))
      ref.push_back(_nsga_struct{static_cast<const Candidate<CBase>&>(tg), tg.tag(), 0, {}});
    #pragma omp parallel if(parallel)
    {
      #pragma omp single
      for(auto& r : ref) {
        _nsga_struct* rr = &r;
        #pragma omp task firstprivate(rr)
        for(auto& s : ref)
          if(rr->r << s.r) {
            rr->q.push_back(&s);
            #pragma omp atomic
            s.dom_cnt++;
          }
      }
    }
    size_t cur_rank = 0;
    /* ref contains the candidates with yet unassigned rank. When this becomes
     * empty, we're done. */
    while(ref.size() > 0) {
      std::vector<_nsga_struct> cur_front{};
      {
        auto it = ref.begin();
        auto end = ref.end();
        while(it != end) {
          /* It *it is nondominated, move it to cur_front and remove from the
           * list. std::list was chosen so that both operations are O(1). */
          if(it->dom_cnt == 0) {
            cur_front.push_back(std::move(*it));
            ref.erase(it++);
          } else
            it++;
        }
      }
      /* Break arrows from members of cur_front and assign final rank to them */
      for(auto& r : cur_front) {
        for(auto q : r.q)
          q->dom_cnt--;
        r.rank = cur_rank;
      }
      cur_rank++;
    }
  }

  std::discrete_distribution<size_t> _nsgaSelect_dist{};
  std::vector<double> _nsgaSelect_probs{};
  size_t _nsgaSelect_last_mod{(size_t)-1};
  double _nsgaSelect_last_bias{};

  template<class Ret, double (*fun)(double, double), class Rng>
  Ret NOINLINE NSGASelect_int(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    assert_sztag();
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("NSGASelect(): Population is empty.");
    if(smp.get_mod_cnt() != _nsgaSelect_last_mod || bias != _nsgaSelect_last_bias) {
      lock.upgrade();
      if(smp.get_mod_cnt() != _nsgaSelect_last_mod + 1)
        _nsga_rate(lock);
      _nsgaSelect_probs.clear();
      _nsgaSelect_probs.reserve(sz);
      for(auto& tg : static_cast<Base&>(*this))
        _nsgaSelect_probs.push_back(1 / fun(tg.tag(), bias));
      _nsgaSelect_dist = std::discrete_distribution<size_t>(_nsgaSelect_probs.begin(), _nsgaSelect_probs.end());
      _nsgaSelect_last_mod = smp.get_mod_cnt();
      _nsgaSelect_last_bias = bias;
    }
    size_t ix = _nsgaSelect_dist(rng);
    return operator[](ix);
  }

public:
  /** \brief Returns the whole population by reference.
   *
   * Note that
   * ```
   * Population::Ref ref = pop.ref();
   * ```
   * is needless because the same effect is equivalently achieved with
   * ```
   * Population::Ref ref(pop);
   * ```
   * However, the availability of an explicit cast given by this function can
   * be used to benefit from C++11's `auto` keyword when Population has not
   * been `typedef`'d and the type name is getting too long. */
  Ref ref() const {
    return *this;
  }

  /** \brief Returns the whole population by reference supplemented by the
   * NSGA metainformation. The original Population does not need to be
   * specially prepared.
   *
   * This creates a reference population supplemented by the NSGA tag (`Tag =
   * sizeof_t`) and populates this by the ranks of each candidate. This is
   * ordinarily an expensive operation (<i>N</i><sup>2</sup> dominance
   * comparisons) that would otherwise be automatically invoked in the first
   * call to NSGASelect or NSGASelect_v. When in a multithreaded context, the
   * thread happening to make the first call would have to block others before
   * the assignment is done. This separate function allows precomputing it
   * using a parallelized algorithm before the thread split.
   *
   * Note that any modifications to the returned population invalidate the
   * NSGA information. If some modifications (e.g., pruning) of the references
   * are necessary, do them on an ordinary ref() and call NSGAref on that as
   * the last step.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
  Population<CBase, size_t, true> NSGAref(bool parallel = true) const {
    internal::read_lock lock(smp);
    Population<CBase, size_t, true> ref(*this);
    ref._nsga_rate(parallel);
    return ref;
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
   * generates a compile-time error in specializations for which this
   * condition is not satisfied.
   *
   * \see Stat
   */
  Stat stat() const {
    assert_double();
    if(this->empty())
      throw std::out_of_range("stat(): Population is empty.");
    double f, sf = 0, sf2 = 0;
    internal::read_lock lock(smp);
    for(auto& c : *this) {
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
    assert_comparable();
    if(!sorted) {
      internal::upgrade_lock up(lock);
      if(sorted) return;
      std::sort(Base::begin(), Base::end());
      sorted = true;
    }
  }

  template<class Rng = decltype(rng)>
  void shuffle(Rng& rng) {
    std::shuffle(Base::begin(), Base::end(), rng);
    sorted = false;
  }

  static void assert_double() {
    static_assert(Candidate<CBase>::Traits::is_float,
        "This method requires the fitness type to be convertible to double.");
  }

  static void assert_comparable() {
    static_assert(Candidate<CBase>::Traits::is_comparable,
        "This method requires the fitness type to implement an operator<.");
  }

  template<class InputIt>
  static void assert_iterator() {
    static_assert(std::is_convertible<typename InputIt::reference, const Candidate<CBase>&>::value,
        "The provided iterator does not return Candidate<CBase>.");
  }

  static void assert_sztag() {
    static_assert(std::is_same<Tag, size_t>::value,
        "This method requires a NSGA rank tag (Tag = size_t).");
  }
}; // class Population

} // namespace gen
