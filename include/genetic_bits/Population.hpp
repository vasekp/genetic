namespace gen {

/** \brief The Population template.
 *
 * From the following parameters, only Candidate is mandatory. For most
 * purposes only `Population<Candidate>` should be needed externally. The
 * other parameters are used mainly for return types of temporary objects. For
 * any values of `Tag` and `is_ref`, `Population<Candidate, Tag, is_ref>` can be
 * implicitly converted to a `Population<Candidate>`.
 *
 * Population can be used as a container of `Candidate`s with read-only access.
 * The functions begin() and end() are exposed, returning random access
 * iterators dereferencable to `const Candidate&` and allowing the iteration
 * patterns
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
 * \tparam Candidate the class describing individual members of this
 * population. Must be derived from gen::Candidate.
 * \tparam Tag an optional supplement class or literal type to accompany each
 * candidate. Used for internal purposes. This does not enter iterations over
 * this population and is not duplicated when copies or references are taken.
 * When supplied, must represent a default-constructible type.
 * \tparam is_ref if set to `true`, this is a reference population. See
 * Population::Ref for more details. */
template<class Candidate, class Tag = internal::empty, bool is_ref = false>
class Population : private std::vector<internal::Tagged<Candidate, Tag, is_ref>> {
  bool sorted = false;
  mutable internal::rw_semaphore smp{};

  typedef decltype(internal::detectFT<Candidate>(nullptr)) _FitnessType;

  typedef std::vector<internal::Tagged<Candidate, Tag, is_ref>> Base;

  typedef internal::cast_iterator<const Candidate&, typename Base::iterator> iterator;
  typedef internal::cast_iterator<const Candidate&, typename Base::const_iterator> const_iterator;

  static_assert(std::is_convertible<Candidate&, const gen::Candidate<_FitnessType>&>::value,
      "The Candidate type needs to be derived from gen::Candidate.");

public:
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
  typedef Population<Candidate, Tag, true> Ref;

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
  Population(Population&& _p): Base(std::move(_p)), sorted(_p.sorted) { }

  /** \brief Creates an empty population but preallocate space for `count`
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

  /** \brief Initializes this population from an iterator range from a
   * container of `Candidate`s. */
  template<class InputIt>
  explicit Population(InputIt first, InputIt last) {
    add(first, last);
  }

  /** \brief Initializes this population from a container of `Candidate`s
   * (e.g., a `std::vector` or another Population).
   * \copydetails add(const Container&) */
  template<class Container>
  explicit Population(const Container& vec) {
    add(std::forward<Container>(vec));
  }

  /** \brief Initializes this population from a container of `Candidate`s
   * using move semantics, leaving the original container empty. */
  template<class Container>
  explicit Population(Container&& vec) {
    add(std::forward<Container>(vec));
  }

#ifndef DOXYGEN
  /* There's no added functionality w.r.t. the above as far as the user is
   * concerned, no need to document */
  template<class Tag_, bool ref_>
  Population(const Population<Candidate, Tag_, ref_>& _p) {
    add(_p);
    sorted = _p.sorted;
  }

  template<class Tag_, bool ref_>
  Population(Population<Candidate, Tag_, ref_>&& _p) {
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
  Population& operator=(Population&& _p) {
    internal::write_lock lock(smp);
    Base::operator=(std::move(_p));
    sorted = _p.sorted;
    return *this;
  }

  /** \brief Copy assignment of a compatible Population.
   * \copydetails add(const Container&) */
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
  const Candidate& operator[](size_t pos) const {
    return static_cast<const Candidate&>(Base::operator[](pos));
  }

  /** \brief Read-only access to a specified element (with bounds checking). */
  const Candidate& at(size_t pos) const {
    return static_cast<const Candidate&>(Base::at(pos));
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
    Base::push_back(std::move(c));
    sorted = false;
  }

  /** \brief Draws `count` candidates from a source function `src`.
   *
   * \param count the number of candidates to generate
   * \param src source function; one of:
   * - `std::function<(const) Candidate ()>`: returning by copy,
   * - `std::function<(const) Candidate& ()>`: returning by reference,
   * - a pointer to function returning `Candidate` or `Candidate&`,
   * - a lambda function returning either.
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
    assert_iterator<InputIt>();
    internal::write_lock lock(smp);
    Base::insert(Base::end(), first, last);
    sorted = false;
  }

  /** \brief Copies all candidates from a container of `Candidate`s.
   *
   * If this population is a reference population, references to all members
   * of the argument are taken, copies are made otherwise. */
  template<class Container>
  void NOINLINE add(const Container& vec) {
    add(vec.begin(), vec.end());
  }

  /** \brief Moves all candidates from a container of `Candidate`s.
   *
   * Moves between `Population`s are only supported if source and destination
   * are both non-reference or both reference and are of the same Candidate type. */
  template<class Container>
#ifdef DOXYGEN
  void
#else
  typename std::enable_if<std::is_rvalue_reference<Container&&>::value, void>::type
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
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged. */
  void rankTrim(size_t newSize) {
    internal::read_lock lock(smp);
    if(size() <= newSize)
      return;
    lock.upgrade();
    ensureSorted(lock);
    const Candidate& dummy = Base::front(); // needed by resize() if size() < newSize which can't happen
    Base::resize(newSize, dummy);           // (otherwise we would need a default constructor)
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
    if(size() <= newSize)
      return;
    lock.upgrade();
    shuffle(rng);
    const Candidate& dummy = Base::front(); // see rankTrim()
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
    internal::read_lock lock(smp);
    if(size() <= minSize)
      return;
    lock.upgrade();
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
    internal::read_lock lock(smp);
    if(size() <= minSize)
      return;
    lock.upgrade();
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

  /** \brief Retrieves a candidate randomly chosen by rank-based selection.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double)` that will receive arguments linearly spaced between
   * `1/size * max` and `max` for candidates ranked `1` through `size` and
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
   * using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * This function may, as a side effect, reorder candidates within this
   * population and thus invalidate references.
   *
   * \param max > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate& rankSelect(double max, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp<const Candidate&>(max, rng);
    else
      return rankSelect_two<const Candidate&, &internal::eval_in_product<fun>>(max, rng);
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
   * <!-- NOTE: don't make this the default value, the two functions could not
   * be distinguished then. -->
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use rankSelect_v(double, Rng&) instead.
   *
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * This function may, as a side effect, reorder candidates within this
   * population and thus invalidate references.
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

  /** \brief Retrieves a candidate randomly chosen by fitness-based selection.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double)` that will receive the fitness of each candidate
   * multiplied by a given constant and its return value will be interpreted
   * as inverse probability, and as such, is expected to be positive and
   * strictly increasing in its argument. This function will be built in at
   * compile time, eliminating a function pointer lookup.  The default value
   * is `std::exp`.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use fitnessSelect_v(double, Rng&) instead.
   *
   * Applicable only if the fitness type of `Candidate` can be converted to
   * `double`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied.
   *
   * \param mult > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate& fitnessSelect(double mult, Rng& rng = rng) {
    return fitnessSelect_int<const Candidate&, &internal::eval_in_product<fun>>(mult, rng);
  }

  /** \brief Retrieves a candidate randomly chosen by fitness-based selection.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double, double)` that will receive the fitness of each candidate
   * as its first argument and the constant `bias` as the second and its
   * return value will be interpreted as inverse probability. As such, is
   * expected to be positive and strictly increasing in its argument. This
   * function will be built in at compile time, eliminating a function pointer
   * lookup. A usual choice `fun` is `std::pow`.
   * <!-- NOTE: don't make this the default value, the two functions could not
   * be distinguished then. -->
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use fitnessSelect_v(double, Rng&) instead.
   *
   * Applicable only if the fitness type of `Candidate` can be converted to
   * `double`. This method generates a compile-time error in
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
  const Candidate& fitnessSelect(double bias, Rng& rng = rng) {
    return fitnessSelect_int<const Candidate&, fun>(bias, rng);
  }

  /** \copybrief fitnessSelect(double, Rng&)
   *
   * Works like fitnessSelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate fitnessSelect_v(double bias, Rng& rng = rng) {
    return fitnessSelect_int<Candidate, &internal::eval_in_product<fun>>(bias, rng);
  }

  /** \copybrief fitnessSelect(double, Rng&)
   *
   * Works like fitnessSelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate fitnessSelect_v(double bias, Rng& rng = rng) {
    return fitnessSelect_int<Candidate, fun>(bias, rng);
  }

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
      for(const Candidate& c : *this)
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
    return operator[](dist(rng));
  }

  /** \copybrief randomSelect(Rng&) const
   *
   * Works like randomSelect(Rng&) const but returns by value. */
  template<class Rng = decltype(rng)>
  Candidate NOINLINE randomSelect_v(Rng& rng = rng) const {
    return randomSelect<Rng, Candidate>(rng);
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
   * Applicable only if the fitness type of `Candidate` allows total ordering
   * using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied. */
#ifdef DOXYGEN
  const Candidate& best() {
#else
  template<class Ret = const Candidate&>
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
          ret.add(static_cast<Candidate&>(tg));
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
  /** \brief Retrieves a candidate randomly chosen by the NSGA algorithm.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double)` that will receive the length of the longest chain of
   * successively dominating candidates of each candidate, multiplied by a
   * given constant, and its return value will be interpreted as inverse
   * probability, and as such, is expected to be positive and strictly
   * increasing in its argument. This function will be built in at compile
   * time, eliminating a function pointer lookup.  The default value is
   * `std::exp`.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use NSGASelect_v(double, Rng&) instead.
   *
   * Applicable only if the `Tag` parameter of this Population is `size_t`.
   * This method generates a compile-time error in specializations for
   * which this condition is not satisfied.
   *
   * \param mult > 0 determines how much nondominated solutions are preferred.
   * Zero would mean no account on dominance rank in the selection process
   * whatsoever. The bigger the value the more nondominated candidates are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate& NSGASelect(double mult, Rng& rng = rng) {
    return NSGASelect_int<const Candidate&, &internal::eval_in_product<fun>>(mult, rng);
  }

  /** \brief Retrieves a candidate randomly chosen by the NSGA algorithm.
   *
   * This method accepts as a template parameter a name of a function
   * `double(double, double)` that will receive the length of the longest
   * chain of successively dominating candidates of each candidate as its
   * first argument and the constant `bias` as the second and its return value
   * will be interpreted as inverse probability. As such, is expected to be
   * positive and strictly increasing in its argument. This function will be
   * built in at compile time, eliminating a function pointer lookup. A usual
   * choice `fun` is `std::pow`.
   * <!-- NOTE: don't make this the default value, the two functions could not
   * be distinguished then. -->
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use NSGASelect_v(double, Rng&) instead.
   *
   * Applicable only if the `Tag` parameter of this Population is `size_t`.
   * This method generates a compile-time error in specializations for
   * which this condition is not satisfied.
   *
   * \param bias > 0 determines how much nondominated solutions are preferred.
   * Zero would mean no account on dominance rank in the selection process
   * whatsoever. The bigger the value the more nondominated candidates are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate& NSGASelect(double bias, Rng& rng = rng) {
    return NSGASelect_int<const Candidate&, fun>(bias, rng);
  }

  /** \copybrief NSGASelect(double, Rng&)
   *
   * Works like NSGASelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate NSGASelect_v(double bias, Rng& rng = rng) {
    return NSGASelect_int<Candidate, &internal::eval_in_product<fun>>(bias, rng);
  }

  /** \copybrief NSGASelect(double, Rng&)
   *
   * Works like NSGASelect(double, Rng&) but returns by value.
   * <!-- FIXME: no way of distinguishing between the two functions in Doxygen ->
   *
   * \returns a copy of a randomly chosen candidate. */
  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate NSGASelect_v(double bias, Rng& rng = rng) {
    return NSGASelect_int<Candidate, fun>(bias, rng);
  }

private:
  struct _nsga_struct {
    const Candidate& r; // reference to a candidate
    size_t& rank;       // reference to its Tag in the original population
    size_t dom_cnt;     // number of candidates dominating this
    std::deque<_nsga_struct*> q;  // candidates dominated by this
  };

  void _nsga_rate(bool parallel = false) {
    internal::write_lock lock(smp);
    _nsga_rate(lock, parallel);
    _nsgaSelect_last_mod = smp.get_mod_cnt();
  }

  void NOINLINE _nsga_rate(internal::rw_lock& lock, bool parallel = false) {
    internal::upgrade_lock up(lock);
    std::list<_nsga_struct> ref{};
    for(auto& tg : static_cast<Base&>(*this))
      ref.push_back(_nsga_struct{static_cast<const Candidate&>(tg), tg.tag(), 0, {}});
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
   * been `typedef`'d and the type name is getting clumsy. */
  Ref ref() const {
    return *this;
  }

  /** \brief Returns the whole population by reference supplemented by the
   * NSGA metainformation. The original Population does not need to be
   * specially prepared.
   *
   * This creates a reference population supplemented by the NSGA tag (`Tag =
   * sizeof_t`) and populates this by the ranks of each candidate. This is
   * ordinarily an expensive operation (\f N^2 \f$ dominance comparisons) that
   * would otherwise be automatically invoked in the first call to NSGASelect or
   * NSGASelect_v. When in a multithreaded context, the thread happening to
   * make the first call would have to block others before the assignment is
   * done. This separate function allows precomputing it using a parallelized
   * algorithm before the thread split.
   *
   * Note that any modifications to the returned population invalidate the
   * NSGA information. If some modifications (e.g., pruning) of the references
   * are necessary, do them on an ordinary ref() and call NSGAref on that as
   * the last step.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
  Population<Candidate, size_t, true> NSGAref(bool parallel = true) const {
    internal::read_lock lock(smp);
    Population<Candidate, size_t, true> ref(*this);
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
    assert_comparable();
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

  static void assert_double() {
    static_assert(std::is_convertible<_FitnessType, double>::value,
        "This method requires the fitness type to be convertible to double.");
  }

  static void assert_comparable() {
    static_assert(internal::comparable<_FitnessType>(0),
        "This method requires the fitness type to implement an operator<.");
  }

  template<class InputIt>
  static void assert_iterator() {
    static_assert(std::is_convertible<typename InputIt::reference, const Candidate&>::value,
        "The provided iterator does not return Candidate.");
  }

  static void assert_sztag() {
    static_assert(std::is_same<Tag, size_t>::value,
        "This method requires a NSGA rank tag (Tag = size_t).");
  }
}; // class Population

} // namespace gen
