namespace gen {

/** \brief The OrdPopulation template, adding functionality dependent on total
 * ordering between candidates to a BasePopulation.
 * \copydetails gen::BasePopulation */
template<class CBase, bool is_ref, class Tag,
  template<class, bool> class Population>
class OrdPopulation: public BasePopulation<CBase, is_ref, Tag, Population> {

  using Base = BasePopulation<CBase, is_ref, Tag, Population>;

  /* Protects: last_sort_mod, rankSelect_* */
  /* Promise: to be only acquired from within a read lock on the Base. */
  mutable internal::rw_semaphore sort_smp{};

  size_t last_sort_mod{(size_t)(~0)};
  std::uniform_real_distribution<double> uniform{0, 1};
  std::discrete_distribution<size_t> rankSelect_dist{};
  std::vector<double> rankSelect_probs{};
  size_t rankSelect_last_sz{};
  double rankSelect_last_bias{};

public:

  using Base::Base;
  using Base::smp;
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::operator[];
  using typename Base::iterator;

  /** \brief Creates an empty population. */
  OrdPopulation() = default;

  /** \copydoc BasePopulation::reserve */
  void reserve(size_t count) {
    Base::reserve(count);
    /* The above command raised smp.mod_cnt by at least 1 but did not disturb
     * sorting. If our population was sorted we reflect that by manually
     * incrementing last_sort_mod. If we weren't level before, or if it rose by
     * more than 1, we won't have is_sorted() afterwards, as intended. */
    internal::write_lock sort_lock{sort_smp};
    ++last_sort_mod;
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
   * \returns a constant reference to the best-of-population.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  const Candidate<CBase>& best() {
    internal::read_lock lock{smp};
    return *best_int(lock, true);
  }

  /** \copybrief best()
   *
   * Works like best() but returns by value.
   *
   * \returns a copy of the best-of-population.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  Candidate<CBase> best_v() {
    internal::read_lock lock{smp};
    return *best_int(lock, true);
  }

  /** \copybrief best()
   *
   * Works like best() but returns an iterator.
   *
   * This function relies on a read lock acquired externally for the
   * population via a PopulationLock. This lock will guard the validity of the
   * returned iterator.
   *
   * \returns an iterator pointing to the best-of-population, end() if the
   * population is empty. */
  iterator best_i(PopulationLock& lock) {
    return best_int(lock.get(), false);
  }

private:

  iterator best_int(internal::rw_lock& lock, bool validate) {
    if(Base::empty()) {
      if(validate)
        throw std::out_of_range("best(): Population is empty.");
      else
        return end();
    }
    else if(is_sorted(lock))
      return {begin()};
    else
      return {std::min_element(begin(), end())};
  }

public:

#ifdef DOXYGEN

  /** \brief Retrieves a candidate randomly chosen by rank-based selection.
   *
   * In determining the probability of each candidate, its rank, rescaled to
   * the range (0, 1], is processed by a given function as described below,
   * and the returned value is interpreted as an inverse probability.
   *
   * The population is sorted from the lowest to the highest fitness value and
   * each candidate is assigned a rank between 1 and \b max. If two or more
   * candidates have an equal fitness, they will be assigned different ranks
   * in an undefined order. This rank is then divided by \b max before passing
   * it to the processing function.
   *
   * The returned reference remains valid until the population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use rankSelect_v() instead.
   *
   * \tparam fun A \b constexpr pointer to a function of signature either
   * <b>double(*)(double)</b> or <b>double(*)(double, double)</b>. In the
   * former case, the argument is <b>x * bias</b>, in the latter case the
   * arguments are \b x and \b bias, where \b x denotes the rescaled rank of
   * each candidate.  It must be positive and strictly increasing in \b x for
   * <b>bias > 0</b>.  This function will be built in at compile time,
   * eliminating a function pointer lookup. The default is \b std::exp, for
   * which an fast specialized algorithm is provided, another usual choice is
   * \b std::pow.
   * \param bias > 0 determines how much low-fitness solutions are preferred.
   * Zero would mean no account on fitness in the selection process
   * whatsoever. The bigger the value the more candidates with low fitness are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& rankSelect(double bias, Rng& rng = rng);

  /** \copybrief rankSelect()
   *
   * Works like rankSelect() but returns by value.
   *
   * \returns a copy of the randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> rankSelect_v(double bias, Rng& rng = rng);

  /** \copybrief rankSelect()
   *
   * Works like rankSelect() but returns an iterator.
   *
   * This function relies on a read lock acquired externally for the
   * population via a PopulationLock. This lock will guard the validity of the
   * returned iterator.
   *
   * This function may need to reorder candidates stored in the population.
   * This will temporarily upgrade the provided lock to a write lock (which
   * implies waiting for other threads to finish their running read
   * operations). This operation is guaranteed not to invalidate any iterators
   * but the contents they refer to changes. Within the calling thread, the
   * population observably changes if rankSelect_i() needed to sort the
   * population.
   *
   * \returns an iterator pointing to the randomly selected candidate, end()
   * if the population is empty. */
  template<class Rng = decltype(rng)>
  iterator rankSelect_i(PopulationLock& lock, double bias, Rng& rng = rng);

#else

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& rankSelect(double max, Rng& rng = rng) {
    internal::read_lock lock{smp};
    if(internal::is_exp<fun>::value)
      return *rankSelect_exp(max, rng, lock, true);
    else
      return *rankSelect_two<
        &internal::eval_in_product<fun>
      >(max, rng, lock, true);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate<CBase>& rankSelect(double bias, Rng& rng = rng) {
    internal::read_lock lock{smp};
    return *rankSelect_two<fun>(bias, rng, lock, true);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> rankSelect_v(double bias, Rng& rng = rng) {
    internal::read_lock lock{smp};
    if(internal::is_exp<fun>::value)
      return *rankSelect_exp(bias, rng, lock, true);
    else
      return *rankSelect_two<
        &internal::eval_in_product<fun>
      >(bias, rng, lock, true);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate<CBase> rankSelect_v(double bias, Rng& rng = rng) {
    internal::read_lock lock{smp};
    return *rankSelect_two<fun>(bias, rng, lock, true);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  iterator rankSelect_i(PopulationLock& lock, double bias, Rng& rng = rng) {
    if(internal::is_exp<fun>::value)
      return rankSelect_exp(bias, rng, lock.get(), false);
    else
      return rankSelect_two<
        &internal::eval_in_product<fun>
      >(bias, rng, lock.get(), false);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  iterator rankSelect_i(PopulationLock& lock, double bias, Rng& rng = rng) {
    return rankSelect_two<fun>(bias, rng, lock.get(), false);
  }

#endif

private:

  template<class Rng>
  NOINLINE iterator rankSelect_exp(double bias, Rng& rng,
      internal::rw_lock& lock, bool validate) {
    double x = uniform(rng);
    size_t sz = size();
    if(sz == 0) {
      if(validate)
        throw std::out_of_range("rankSelect(): Population is empty.");
      else
        return end();
    }
    ensure_sorted(lock);
    // Bug in GCC: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
    if(x == 1)
      return end() - 1;
    else
      return begin() + ((int)(-log(1 - x + x*exp(-bias))/bias*sz));
  }

  template<double (*fun)(double, double), class Rng>
  NOINLINE iterator rankSelect_two(double bias, Rng& rng,
      internal::rw_lock& lock, bool validate) {
    size_t sz = size();
    if(sz == 0) {
      if(validate)
        throw std::out_of_range("rankSelect(): Population is empty.");
      else
        return end();
    }
    internal::read_lock sort_lock{sort_smp};
    if(sz != rankSelect_last_sz || bias != rankSelect_last_bias) {
      sort_lock.upgrade();
      rankSelect_probs.clear();
      rankSelect_probs.reserve(sz);
      for(size_t i = 0; i < sz; i++)
        rankSelect_probs.push_back(1 / fun((double)(i+1) / sz, bias));
      rankSelect_dist = std::discrete_distribution<size_t>
        (rankSelect_probs.begin(), rankSelect_probs.end());
      rankSelect_last_sz = sz;
      rankSelect_last_bias = bias;
    }
    ensure_sorted(lock);
    return begin() + rankSelect_dist(rng);
  }

public:

  /** \brief Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged. */
  NOINLINE void rankTrim(size_t newSize) {
    internal::read_lock lock{smp};
    if(!lock.upgrade_if([newSize,this]() -> bool { return size() > newSize; }))
      return;
    if(size() == 0)
      return;
    ensure_sorted(lock);
    auto& dummy = Base::first(); // see BasePopulation::randomTrim()
    Base::as_vec().resize(newSize, dummy);
    // See reserve()
    // No sort_lock needed (have write_lock(smp))
    ++last_sort_mod;
  }

  /** \brief Sorts the population by fitness for faster response of future
   * rankSelect() calls.
   *
   * This happens by default whenever rankSelect() is called after the
   * population has been modified. However, since the sorting can be a rather
   * expensive operation, in a multithreaded setting this would mean all other
   * threads have to wait for the calling thread to finish the sorting before
   * proceeding. This method can be called before the work is split between
   * threads. Parallel sorting is currently not supported but this can still
   * be used to guarantee a better balanced workload for individual threads. */
  NOINLINE void sort() {
    internal::read_lock lock{smp};
    ensure_sorted(lock);
  }

  /** \copybrief sort()
   *
   * A variant of sort() for use in blocks protected by PopulationLock.
   *
   * This function may need to reorder candidates stored in the population.
   * This will temporarily upgrade the provided lock to a write lock (which
   * implies waiting for other threads to finish their running read
   * operations). This operation is guaranteed not to invalidate any iterators
   * but the contents they refer to changes. Within the calling thread, the
   * population observably changes if the order needed to be updated. */
  NOINLINE void sort(PopulationLock& lock) {
    ensure_sorted(lock.get());
  }

private:

  bool is_sorted(const internal::rw_lock&) const {
    internal::read_lock sort_lock{sort_smp};
    return smp.get_mod_cnt() == last_sort_mod;
  }

  void ensure_sorted(internal::rw_lock& lock) {
    if(!is_sorted(lock)) {
      internal::upgrade_lock up{lock};
      // No one else can be reading or modifying this now (⇐ promise)
      ++last_sort_mod;
      if(is_sorted(lock))
        return;
      std::sort(Base::as_vec().begin(), Base::as_vec().end());
      last_sort_mod = smp.get_mod_cnt();
    }
  }

}; // class OrdPopulation

} // namespace gen
