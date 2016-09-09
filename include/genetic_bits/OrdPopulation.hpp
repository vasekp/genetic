namespace gen {

/** \brief The OrdPopulation template, adding functionality dependent on total
 * ordering between candidates to a BasePopulation.
 * \copydetails gen::BasePopulation */
template<class CBase, bool is_ref, class Tag, template<class, bool> class Population>
class OrdPopulation: public BasePopulation<CBase, is_ref, Tag, Population> {

  using Base = BasePopulation<CBase, is_ref, Tag, Population>;
  using Base2 = internal::PBase<CBase, is_ref, Tag>;

  /* Protects: last_sort_mod, rankSelect_* */
  /* Promise: to be only acquired from within a read lock on the Base. */
  mutable internal::rw_semaphore sort_smp{};

  size_t last_sort_mod{(size_t)(~0)};
  std::uniform_real_distribution<double> uniform{0, 1};
  std::discrete_distribution<size_t> rankSelect_dist{};
  std::vector<double> rankSelect_probs{};
  size_t rankSelect_last_sz{};
  double rankSelect_last_bias{};

  /* Befriend all compatible OrdPopulations for access to their is_sorted() */
  template<class, bool, class, template<class, bool> class>
  friend class OrdPopulation;

public:

  using Base::Base;
  using Base::smp;
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::operator[];

  /** \brief Creates an empty population. */
  OrdPopulation() = default;

  /** \copydoc BasePopulation::reserve */
  void reserve(size_t count) {
    Base::reserve(count);
    /* The above command raised smp.mod_cnt by at least 1 but did not disturb
     * sorting. If our population was sorted we reflect that by manually
     * incrementing last_sort_mod. If we weren't level before, or if it rose by
     * more than 1, we won't have is_sorted() afterwards, as intended. */
    internal::write_lock sort_lock(sort_smp);
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
   * allows this, use best_v() instead. */
#ifdef DOXYGEN
  const Candidate<CBase>& best() {
#else
  template<class Ret = const Candidate<CBase>&>
  Ret NOINLINE best() {
#endif
    internal::read_lock lock(smp);
    if(this->empty())
      throw std::out_of_range("best(): Population is empty.");
    if(is_sorted(lock))
      return Base2::front();
    else
      return *std::min_element(begin(), end());
  }

  /** \copybrief best()
   *
   * Works like best() but returns by value. */
  Candidate<CBase> best_v() {
    return best<Candidate<CBase>>();
  }

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

  template<class Ret, class Rng>
  Ret NOINLINE rankSelect_exp(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    double x = uniform(rng);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("rankSelect(): Population is empty.");
    ensure_sorted(lock);
    if(x == 1)
      return Base2::back();
    else
      return operator[]((int)(-log(1 - x + x*exp(-bias))/bias*sz));
  }

  template<class Ret, double (*fun)(double, double), class Rng>
  Ret NOINLINE rankSelect_two(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("rankSelect(): Population is empty.");
    internal::read_lock sort_lock(sort_smp);
    if(sz != rankSelect_last_sz || bias != rankSelect_last_bias) {
      sort_lock.upgrade();
      rankSelect_probs.clear();
      rankSelect_probs.reserve(sz);
      for(size_t i = 0; i < sz; i++)
        rankSelect_probs.push_back(1 / fun((double)(i+1) / sz, bias));
      rankSelect_dist = std::discrete_distribution<size_t>(rankSelect_probs.begin(), rankSelect_probs.end());
      rankSelect_last_sz = sz;
      rankSelect_last_bias = bias;
    }
    ensure_sorted(lock);
    return operator[](rankSelect_dist(rng));
  }

public:

  /** \brief Reduces the population to a maximum size given by the argument,
   * dropping the worst part of the sample.
   *
   * \param newSize the maximum desired size of the population. If this bound
   * is satisfied, the population is unchanged. */
  void rankTrim(size_t newSize) {
    internal::read_lock lock(smp);
    bool was_sorted = is_sorted(lock);
    if(!lock.upgrade_if([newSize,this]() -> bool { return size() > newSize; }))
      return;
    if(size() == 0)
      return;
    ensure_sorted(lock);
    auto& dummy = Base2::front(); // see BasePopulation::randomTrim()
    Base2::resize(newSize, dummy);
    if(was_sorted)
      set_sorted(lock);
  }

private:

  bool is_sorted(const internal::rw_lock&) const {
    internal::read_lock sort_lock(sort_smp);
    return smp.get_mod_cnt() == last_sort_mod;
  }

  void set_sorted(const internal::rw_lock&) {
    internal::write_lock sort_lock(sort_smp);
    last_sort_mod = smp.get_mod_cnt();
  }

  void ensure_sorted(internal::rw_lock& lock) {
    if(!is_sorted(lock)) {
      internal::upgrade_lock up(lock);
      // No one else can be reading or modifying this now (⇐ promise)
      ++last_sort_mod;
      if(is_sorted(lock))
        return;
      std::sort(Base2::begin(), Base2::end());
      set_sorted(lock);
    }
  }

}; // class OrdPopulation

} // namespace gen
