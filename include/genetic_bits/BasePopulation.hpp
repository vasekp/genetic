namespace gen {

namespace internal {

template<class CBase, bool is_ref, class Tag>
using PBase = std::vector<internal::CandidateTagged<CBase, is_ref, Tag>>;

} // namespace internal


/** \brief The BasePopulation template, covering functionality common to all
 * the derived Population classes.
 *
 * This is an inner class of the framework, not suitable to be used directly
 * by applications.
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
template<class CBase, bool is_ref, class Tag,
  template<class, bool> class Population>
#ifndef DOXYGEN
class BasePopulation: protected internal::PBase<CBase, is_ref, Tag> {
#else
class BasePopulation {
#endif

  using Base = internal::PBase<CBase, is_ref, Tag>;

  /* For access to begin(), end(), size() etc. via CRTP */
  friend class OrdPopulation<CBase, is_ref, Tag, Population>;
  friend class FloatPopulation<CBase, is_ref, Tag, Population>;
  friend class DomPopulation<CBase, is_ref, Tag, Population>;

protected:

#ifndef DOXYGEN
  using iterator = internal::CTIterator<typename Base::iterator>;
  using const_iterator = internal::CTIterator<typename Base::const_iterator>;

  mutable internal::rw_semaphore smp{};

  /* Allow access to smp */
  friend class PopulationLock;
#endif

public:

  /** \brief A corresponding "reference population", a helper type for
   * functions returning a selection from an existing population.
   *
   * This is the return type of functions that return a subset of an existing
   * population by reference as it holds its members by reference rather than
   * by value. Internally it is another Population specialization with a
   * matching \b CBase type, so the same functions including adding or erasing
   * can be called on it despite its dependent nature and assignments between
   * Population and \b Population::Ref objects are allowed and behave as
   * expected, see add(const Container&) and add(Container&&).
   *
   * It is guaranteed that \b Population::Ref::Ref is identical to
   * \b Population::Ref, which makes it convenient to chain selection
   * functions, e.g.
   * ```
   * pop.randomSelect(5).front().randomSelect()
   * ```
   * for a simple multi-objective tournament selection.
   *
   * It is the user's responsibility not to use the references beyond their
   * scope. They are invalidated by operations which modify the original
   * population, namely all operations adding, removing, and reordering its
   * elements.  This includes functions like \link
   * gen::OrdPopulation::rankSelect() rankSelect() \endlink, which needs to
   * sort the contents for its operation, and reserve(), which may move the
   * contents to a new memory location. */
  using Ref = Population<CBase, true>;

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
  using Val = Population<CBase, false>;

  /** \brief Creates an empty population. */
  BasePopulation() = default;

#ifndef DOXYGEN
  /* The copy and move constructors: trivial but we need explicit definition
   * because the semaphore can't be default copied or moved. They are defined
   * implicitly in any other *Population so let's not draw extra attention to
   * these in particular. */

  BasePopulation(const BasePopulation& p): Base(p) { }

  BasePopulation(BasePopulation&& p) noexcept: Base(std::move(p)) { }
#endif

  /** \brief Creates an empty population but preallocates space for \b count
   * candidates. */
  explicit BasePopulation(size_t count) {
    Base::reserve(count);
  }

  /** \brief Creates a population of size \b count whose candidates are results
   * of calls to the source function \b src.
   * \copydetails add(size_t, Source, bool) */
  template<class Source>
  explicit BasePopulation(size_t count, Source src, bool parallel = true) {
    add(count, src, parallel);
  }

  /** \brief Initializes this population from an iterator range from a
   * container of <b>Candidate</b>s or <b>CBase</b>s. */
  template<class InputIt>
  explicit BasePopulation(InputIt first, InputIt last) {
    add(first, last);
  }

  /** \brief Initializes this population from a container of <b>Candidate</b>s
   * or <b>CBase</b>s (e.g., a [**std::vector**] (http://en.cppreference.com/w/
   * cpp/container/vector) or another Population).
   *
   * \copydetails add(const Container&) */
  template<class Container>
  explicit BasePopulation(const Container& vec) {
    add(vec);
  }

  /** \brief Initializes this population from a container of <b>Candidate</b>s
   * or <b>CBase</b>s using move semantics, leaving the original container
   * empty. */
  template<class Container>
  explicit BasePopulation(Container&& vec) {
    add(std::forward<Container>(vec));
  }

  /* NB: the difference in the following two constructors is that they are
   * implicit. */

  /** \brief Initializes this population from a compatible population.
   * \copydetails add(const Container&) */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation(const BasePopulation<CBase, ref_, Tag_, Pop_>& p) {
    add(p);
  }

  /** \brief Initializes this population from a compatible population
   * using move semantics, leaving the original container empty. */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation(BasePopulation<CBase, ref_, Tag_, Pop_>&& p) {
    add(std::move(p));
  }

#ifndef DOXYGEN
  /* Copy and move assignment operators. Ditto as c&m constructors. */

  BasePopulation& operator=(const BasePopulation& p) {
    internal::write_lock lock{smp};
    Base::operator=(p);
    return *this;
  }

  BasePopulation& operator=(BasePopulation&& p) noexcept {
    internal::write_lock lock{smp};
    Base::operator=(std::move(p));
    return *this;
  }
#endif

  /** \brief Copy assignment of a compatible population.
   * \copydetails add(const Container&) */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation& operator=(const BasePopulation<CBase, ref_, Tag_, Pop_>& p) {
    internal::write_lock lock{smp};
    Base::clear();
    Base::insert(Base::end(), p.begin(), p.end());
    return *this;
  }

  /** \brief Move assignment of a compatible population. */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  BasePopulation& operator=(BasePopulation<CBase, ref_, Tag_, Pop_>&& p) {
    internal::write_lock lock{smp};
    Base::clear();
    move_add_unguarded(p);
    return *this;
  }

#ifdef DOXYGEN
  /** \brief Returns the current count of candidates. */
  size_t size() const;
#else
  using Base::size;
#endif

  /** \brief Empties the population. */
  void clear() {
    internal::write_lock lock{smp};
    Base::clear();
  }

  /** \brief Reserves space for \b count candidates.
   *
   * If \b count is larger than the actual size of the population, all
   * references may be invalidated. */
  void reserve(size_t count) {
    // Does not count as a modification
    // Invalidates iterators, but they may not be kept between locks anyway.
    internal::write_lock lock{smp, false};
    Base::reserve(count);
  }

  /** \brief Returns an iterator to the beginning.
   *
   * This iterator dereferences to a \link gen::Candidate const
   * Candidate<CBase>&\endlink.
   *
   * Note that using iterators to access the individual candidates in a
   * Population circumvents the memory locking mechanisms. In any case
   * references to members of a population are invalidated, so are any
   * currently stored iterators.
   *
   * \see PopulationLock */
  iterator begin() {
    return {Base::begin()};
  }

  /** \brief Returns the past-the-end iterator. */
  iterator end() {
    return {Base::end()};
  }

#ifndef DOXYGEN
  /* There is no user-observable difference between an iterator and a
   * const_iterator. Both dereference to a const Candidate&. The difference
   * arises internally when a Population is moved because iterator can be
   * turned into a move_iterator and const_iterator can not. */

  const_iterator begin() const {
    return {Base::begin()};
  }

  const_iterator end() const {
    return {Base::end()};
  }

protected:
  /* We also don't document protected members.
   * This is accessed from OrdPopulation, other kinds of Population are
   * insensitive to ordering. */

  std::reverse_iterator<iterator> rbegin() {
    return {Base::rbegin()};
  }

  std::reverse_iterator<iterator> rend() {
    return {Base::rend()};
  }

  std::reverse_iterator<const_iterator> rbegin() const {
    return {Base::rbegin()};
  }

  std::reverse_iterator<const_iterator> rend() const {
    return {Base::rend()};
  }
#endif

public:

  /** \brief Read-only access to a specified element by reference
   * (no bounds checking).
   *
   * Does not lock the population for read access. If another thread
   * simultaneously modifies the population the results are undefined. */
  const Candidate<CBase>& operator[](size_t pos) const {
    return static_cast<const Candidate<CBase>&>(Base::operator[](pos));
  }

  /** \brief Read-only access to a specified element by reference
   * (with bounds checking). */
  const Candidate<CBase>& at(size_t pos) const {
    internal::read_lock lock{smp};
    return static_cast<const Candidate<CBase>&>(Base::at(pos));
  }

  /** \brief Returns a specified element by value (with bounds checking) */
  Candidate<CBase> at_v(size_t pos) const {
    internal::read_lock lock{smp};
    return static_cast<const Candidate<CBase>>(Base::at(pos));
  }

protected:

#ifndef DOXYGEN
  /* Used for iterating directly over std::vector<CandidateTagged>. */
  Base& as_vec() {
    return static_cast<Base&>(*this);
  }

  const Base& as_vec() const {
    return static_cast<const Base&>(*this);
  }

  const Candidate<CBase>& first() {
    return Base::front();
  }

  const Candidate<CBase>& last() {
    return Base::back();
  }
#endif

public:

  /** \brief Adds a new candidate. */
  void add(const Candidate<CBase>& c) {
    internal::write_lock lock{smp};
    Base::push_back(c);
  }

  /** \brief Pushes back a new candidate using the move semantics. */
  void add(Candidate<CBase>&& c) {
    internal::write_lock lock{smp};
    Base::push_back(std::move(c));
  }

  /** \brief Draws \b count candidates from a source function \b src.
   *
   * \param count the number of candidates to generate
   * \param src source function; can be any callable object (e.g., a
   * [**std::function**] (http://en.cppreference.com/w/cpp/utility/functional/
   * function), a function pointer, or a lambda function) returning
   * either of \link gen::Candidate Candidate<CBase> \endlink or \b CBase and
   * either by value or by reference (a copy will be taken). In many cases the
   * function call can be inlined by the optimizer if known at compile time.
   * \param parallel controls parallelization using OpenMP (on by default) */
  template<class Source>
  NOINLINE void add(size_t count, Source src, bool parallel = true) {
    internal::write_lock lock{smp};
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
      Base::insert(Base::end(),
          std::make_move_iterator(tmp.begin()),
          std::make_move_iterator(tmp.end()));
    }
  }

  /** \brief Copies an iterator range from a container of <b>Candidate</b>s. */
  template<class InputIt>
  void add(InputIt first, InputIt last) {
    internal::write_lock lock{smp};
    Base::insert(Base::end(), first, last);
  }

  /** \brief Copies all candidates from a container of <b>Candidate</b>s.
   *
   * If this population is a reference population, references to all members
   * of the argument are taken, copies are made otherwise. */
  template<class Container>
#ifdef DOXYGEN
  void
#else
  typename std::enable_if<internal::is_container<Container>(0), void>::type
#endif
  add(const Container& vec) {
    add(vec.begin(), vec.end());
  }

  /** \brief Moves all candidates from a container of <b>Candidate</b>s.
   *
   * Moves between populations are only supported if source and destination
   * are both non-reference or both reference and are of the same \b CBase
   * type. */
  template<class Container>
#ifdef DOXYGEN
  void
#else
  typename std::enable_if<
      internal::is_container<Container>(0) &&
      std::is_rvalue_reference<Container&&>::value,
    void>::type
#endif
  add(Container&& vec) {
    internal::write_lock lock{smp};
    move_add_unguarded(std::forward<Container>(vec));
  }

private:

  template<class Container>
  void move_add_unguarded(Container&& vec) {
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
    internal::read_lock lock{smp};
    if(!lock.upgrade_if([newSize,this]() -> bool { return size() > newSize; }))
      return;
    if(size() == 0)
      return;
    shuffle(rng);
    /* A reference is needed for resize() to fill up new space if newSize >
     * size(). (There's an overload without this argument but that relies on a
     * default constructor which we don't require). This is prevented by the
     * lock so we can pass any reference at hand. Front is guaranteed to exist
     * after excluding an empty population above. */
    auto& dummy = Base::front();
    Base::resize(newSize, dummy);
  }

  /** \brief Reduces the population by selective removal of candidates.
   *
   * Candidates are tested against a given predicate. If a candidate \b c
   * satisfies \b pred(c) it is removed from the population. A minimum number
   * of candidates can be set; if so, the procedure stops when enough
   * candidates have been removed to satisfy this bound.
   *
   * \param test a boolean function accepting a constant \link gen::Candidate
   * Candidate<CBase>\endlink reference. If the return value is \b true the
   * candidate is removed from the population.
   * \param minSize a minimum number of candidates to be kept if possible. If
   * zero (the default value), all candidates satisfying the predicate are
   * removed.
   * \param rng the random number generator, or gen::rng by default. Unused
   * if \b randomize is \b false. */
  template<class Rng = decltype(rng)>
  NOINLINE void prune(bool (*test)(const Candidate<CBase>&),
      size_t minSize = 0, Rng& rng = rng) {
    internal::read_lock lock{smp};
    if(!lock.upgrade_if([minSize,this]() -> bool { return size() > minSize; }))
      return;
    size_t sz = size();
    for(size_t i = 0; i < sz; )
      if(test(operator[](i)))
        std::swap(Base::operator[](i), Base::operator[](--sz));
      else
        i++;
    Base::erase(Base::begin() + sz, Base::end());
  }

  /** \brief Reduces the population by selective removal of candidates.
   *
   * Candidates are tested for similarity according to a provided criterion
   * function. If a pair of candidates (\b a, \b b) satisfies the test, only
   * \b a is kept. A minimum number of candidates can be set; if so, the
   * procedure stops when enough candidates have been removed to satisfy this
   * bound.
   *
   * \param test a boolean function accepting two constant \link
   * gen::Candidate Candidate<CBase>\endlink references. Should be symmetric
   * in its arguments. If the return value is \b true the latter candidate is
   * removed from the population.
   * \param minSize a minimum number of candidates to be kept if possible. If
   * zero (the default value), all duplicates are removed.
   * \param randomize whether to randomly shuffle the sample prior to pruning
   * (this is the default). If \b false then earlier appearing candidates are
   * preferred in survival.
   * \param rng the random number generator, or gen::rng by default. Unused
   * if \b randomize is \b false. */
  template<class Rng = decltype(rng)>
  NOINLINE void prune(
      bool (*test)(const Candidate<CBase>&, const Candidate<CBase>&),
      size_t minSize = 0, bool randomize = true, Rng& rng = rng) {
    internal::read_lock lock{smp};
    if(!lock.upgrade_if([minSize,this]() -> bool { return size() > minSize; }))
      return;
    if(randomize)
      shuffle(rng);
    size_t sz = size();
    for(size_t i = 0; i < sz - 1; i++) {
      const Candidate<CBase>& op1 = operator[](i);
      for(size_t j = i + 1; j < sz; ) {
        if(test(op1, operator[](j)))
          std::swap(Base::operator[](j), Base::operator[](--sz));
        else
          j++;
      }
    }
    Base::erase(Base::begin() + sz, Base::end());
  }

  /** \brief Retrieves a candidate chosen using uniform random selection.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use randomSelect_v(Rng&) instead.
   *
   * \returns a constant reference to a randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<class Rng = decltype(rng)>
#ifdef DOXYGEN
  const Candidate<CBase>& randomSelect(Rng& rng = rng) {
#else
  typename std::enable_if<
    internal::is_URNG<Rng>(0),
    const Candidate<CBase>&
  >::type randomSelect(Rng& rng = rng) {
#endif
    internal::read_lock lock{smp};
    return *randomSelect_int(rng, lock, true);
  }

  /** \copybrief randomSelect(Rng&)
   *
   * Works like randomSelect(Rng&) but returns by value.
   *
   * \returns a copy of the randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<class Rng = decltype(rng)>
#ifdef DOXYGEN
  Candidate<CBase> randomSelect_v(Rng& rng = rng) {
#else
  typename std::enable_if<internal::is_URNG<Rng>(0), Candidate<CBase>>::type
  randomSelect_v(Rng& rng = rng) {
#endif
    internal::read_lock lock{smp};
    return *randomSelect_int(rng, lock, true);
  }

  /** \copybrief randomSelect(Rng&)
   *
   * Works like randomSelect(Rng&) but returns an iterator.
   *
   * This function relies on a read lock acquired externally for the
   * population via a PopulationLock. This lock will guard the validity of the
   * returned iterator.
   *
   * \returns an iterator pointing to the randomly selected candidate, end()
   * if the population is empty. */
  template<class Rng = decltype(rng)>
  iterator randomSelect_i(PopulationLock& lock, Rng& rng = rng) {
    return randomSelect_int(rng, lock.get(), false);
  }

private:

  template<class Rng>
  iterator randomSelect_int(Rng& rng, internal::rw_lock&, bool validate) {
    size_t sz = size();
    if(sz == 0) {
      if(validate)
        throw std::out_of_range("randomSelect(): Population is empty.");
      else
        return end();
    }
    std::uniform_int_distribution<size_t> dist{0, sz - 1};
    return begin() + dist(rng);
  }

public:

  /** \brief Randomly selects \b k different candidates. If <b>k ≥ size()</b>,
   * the entire population is returned.
   *
   * The returned \link Ref \endlink remains valid until the original
   * population is modified.  Therefore there is a risk of invalidating it in
   * a multi-threaded program if another thread concurrently modifies the
   * population. If your code allows this, use randomSelect_v(size_t, Rng&)
   * instead. */
#ifdef DOXYGEN
  template<class Rng = decltype(rng)>
  Ref randomSelect(size_t k, Rng& rng = rng) {
#else
  template<class Rng = decltype(rng), class Ret = Ref>
  NOINLINE Ret randomSelect(size_t k, Rng& rng = rng) {
#endif
    internal::read_lock lock{smp};
    size_t sz = size();
    if(k >= sz)
      return Ret{*this};
    std::vector<size_t> idx(sz);
    /* Fisher-Yates intentionally without initialization! */
    for(size_t i = 0; i < k; i++) {
      size_t d = rng() % (sz - i), j = i+d;
      std::swap(idx[i], idx[j]);
      idx[j] -= d;
      idx[i] += i + d;
    }
    idx.resize(k);
    Ret ret{k};
    for(auto i : idx)
      ret.add(operator[](i));
    return ret;
  }

  /** \copybrief randomSelect(size_t, Rng&)
   *
   * Works like randomSelect(size_t, Rng&) but returns an independent
   * population. */
  template<class Rng = decltype(rng)>
  Val randomSelect_v(size_t k, Rng& rng = rng) {
    return randomSelect<Val>(k, rng);
  }

private:

  template<class Rng = decltype(rng)>
  void shuffle(Rng& rng) {
    std::shuffle(Base::begin(), Base::end(), rng);
  }

}; // class BasePopulation<CBase, is_ref, Tag, Population>

} // namespace gen
