namespace gen {

namespace internal {

template<class CBase, bool is_ref, class Tag>
using PBase = std::vector<internal::CandidateTagged<CBase, is_ref, Tag>>;

} // namespace internal


/** \brief The BasePopulation template, covering functionality common to all
 * the derived Population classes.
 *
 * \tparam CBase the base class of the member candidates of this population.
 * See Candidate for details.
 * \tparam is_ref if set to \b true, this is a reference population. See
 * BasePopulation::Ref for more details.
 * \tparam Tag a class or literal type to accompany each candidate, used for
 * internal purposes. This does not enter iterations over this population and
 * is not duplicated when copies or references are taken. Must represent a
 * default-constructible type.
 * \tparam Population the outer Population type this should behave like.
 * Controls the return type of selection functions. */
template<class CBase, bool is_ref, class Tag, template<class, bool> class Population>
class BasePopulation: protected internal::PBase<CBase, is_ref, Tag> {

  typedef internal::PBase<CBase, is_ref, Tag> Base;

  typedef internal::CTIterator<typename Base::iterator> iterator;
  typedef internal::CTIterator<typename Base::const_iterator> const_iterator;

protected:

#ifndef DOXYGEN
  mutable internal::rw_semaphore smp{};
#endif

public:

  /** \brief A corresponding "reference population", a helper type for
   * functions returning a selection from an existing population.
   *
   * This is the return type of functions that return a subset of an existing
   * population by reference as it holds its members by reference rather than
   * by value. Internally it is another Population specialization with a
   * matching CBase type, so the same functions including adding or erasing
   * can be called on it despite its temporary nature and assignments between
   * Population and \b Population::Ref objects are allowed and behave as
   * expected, see add(Container&).
   *
   * It is guaranteed that \b Population::Ref::Ref is identical to
   * \b Population::Ref, which makes it convenient to chain selection functions,
   * e.g.
   * ```
   * pop.randomSelect(5).front().randomSelect()
   * ```
   * for a simple multi-objective tournament selection.
   *
   * It is the user's responsibility not to use the references beyond their
   * scope. They are invalidated by operations which modify the original
   * population, namely all operations adding, removing, and reordering its
   * elements.  This includes \link gen::OrdPopulation::rankSelect()
   * rankSelect() \endlink, which needs to sort the contents for its
   * operation, and reserve(), which may move the contents to a new memory
   * location. */
  typedef Population<CBase, true> Ref;

  /** \brief A corresponding "value population", a helper type for functions
   * returning a selection from an existing population by value.
   *
   * In the case of a \link gen::Population Population<CBase> \endlink this
   * will be identical to the calling class. However, when requested from a
   * \link Ref reference population \endlink, this will remove the reference
   * character thereof, so that elements stored in this type may outlive their
   * original source.
   *
   * It is guaranteed that \b Population::Val::Val = \b Population::Ref::Val =
   * \b Population::Val. */
  typedef Population<CBase, false> Val;

  /** \brief Creates an empty population. */
  BasePopulation() = default;

#ifndef DOXYGEN
  /* The copy and move constructors: trivial but we need explicit definition
   * because the semaphore can't be default copied or moved. They are defined
   * implicitly in any other *Population so let's not draw extra attention to
   * these in particular. */

  BasePopulation(const BasePopulation& _p): Base(_p) { }

  BasePopulation(BasePopulation&& _p) noexcept: Base(std::move(_p)) { }
#endif

  /** \brief Creates an empty population but preallocates space for `count`
   * candidates. */
  explicit BasePopulation(size_t count) {
    Base::reserve(count);
  }

  /** \brief Creates a population of size `count` whose candidates are results
   * of calls to the source function `src`.
   * \copydetails add(size_t, Source, bool) */
  template<class Source>
  explicit BasePopulation(size_t count, Source src, bool parallel = true) {
    add(count, src, parallel);
  }

  /** \brief Initializes this population from an iterator range from a
   * container of `Candidate`s or `CBase`s. */
  template<class InputIt>
  explicit BasePopulation(InputIt first, InputIt last) {
    add(first, last);
  }

  /** \brief Initializes this population from a container of `Candidate`s or
   * `CBase`s (e.g., a `std::vector` or another BasePopulation).
   * \copydetails add(const Container&) */
  template<class Container>
  explicit BasePopulation(const Container& vec) {
    add(vec);
  }

  /** \brief Initializes this population from a container of `Candidate`s or
   * `CBase`s using move semantics, leaving the original container empty. */
  template<class Container>
  explicit BasePopulation(Container&& vec) {
    add(std::forward<Container>(vec));
  }

  /* NB: the difference in the following two constructors is that they are
   * implicit. */

  /** \brief Initializes this population from a compatible population.
   * \copydetails add(const Container&) */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation(const BasePopulation<CBase, ref_, Tag_, Pop_>& _p) {
    add(_p);
  }

  /** \brief Initializes this population from a compatible population
   * using move semantics, leaving the original container empty. */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation(BasePopulation<CBase, ref_, Tag_, Pop_>&& _p) {
    add(std::move(_p));
  }

#ifndef DOXYGEN
  /* Copy and move assignment operators. Ditto as c&m constructors. */

  BasePopulation& operator=(const BasePopulation& _p) {
    internal::write_lock lock(smp);
    Base::operator=(_p);
    return *this;
  }

  BasePopulation& operator=(BasePopulation&& _p) noexcept {
    internal::write_lock lock(smp);
    Base::operator=(std::move(_p));
    return *this;
  }
#endif

  /** \brief Copy assignment of a compatible population.
   * \copydetails add(const Container&) */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation& operator=(const BasePopulation<CBase, ref_, Tag_, Pop_>& _p) {
    internal::write_lock lock(smp);
    Base::clear();
    Base::insert(Base::end(), _p.begin(), _p.end());
    return *this;
  }

  /** \brief Move assignment of a compatible population. */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation& operator=(BasePopulation<CBase, ref_, Tag_, Pop_>&& _p) {
    internal::write_lock lock(smp);
    Base::clear();
    move_add_unguarded(_p);
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
  }

  /** \brief Pushes back a new candidate using the move semantics. */
  void add(Candidate<CBase>&& c) {
    internal::write_lock lock(smp);
    Base::push_back(std::move(c));
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
  }

  /** \brief Copies an iterator range from a container of candidates. */
  template<class InputIt>
  void add(InputIt first, InputIt last) {
    internal::write_lock lock(smp);
    Base::insert(Base::end(), first, last);
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
   * Moves between populations are only supported if source and destination
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
  }

private:

  template<class Container>
  void NOINLINE move_add_unguarded(Container&& vec) {
    Base::insert(Base::end(),
        internal::move_iterator<decltype(vec.begin())>(vec.begin()),
        internal::move_iterator<decltype(vec.end())>(vec.end()));
    vec.clear();
  }

public:

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
      throw std::out_of_range("randomSelect(): BasePopulation is empty.");
    std::uniform_int_distribution<size_t> dist{0, sz - 1};
    return operator[](dist(rng));
  }

  /** \copybrief randomSelect(Rng&) const
   *
   * Works like randomSelect(Rng&) const but returns by value. */
  template<class Rng = decltype(rng)>
  Candidate<CBase> NOINLINE randomSelect_v(Rng& rng = rng) const {
    return randomSelect<Candidate<CBase>>(rng);
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
   * BasePopulation. */
  template<class Rng = decltype(rng)>
  Val NOINLINE randomSelect_v(size_t k, Rng& rng = rng) const {
    return randomSelect<Val>(k, rng);
  }

private:

  template<class Rng = decltype(rng)>
  void shuffle(Rng& rng) {
    std::shuffle(Base::begin(), Base::end(), rng);
  }

}; // class BasePopulation

} // namespace gen
