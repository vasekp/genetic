namespace gen {

template<class CBase, bool is_ref, class Tag, template<class, bool> class Population>
class OrdPopulation : public BasePopulation<CBase, is_ref, Tag, Population> {

  typedef BasePopulation<CBase, is_ref, Tag, Population> Base;
  typedef internal::PBase<CBase, is_ref, Tag> Base2;

  size_t last_sort_mod{(size_t)(~0)};

public:

  using Base::Base;
  using Base::smp;
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::operator[];

  /* Befriend all compatible OrdPopulations */
  template<class, bool, class, template<class, bool> class>
  friend class OrdPopulation;

  /** \brief Creates an empty population. */
  OrdPopulation() = default;

  /** \brief The copy constructor. */
  OrdPopulation(const OrdPopulation& _p): Base(_p) {
    if(_p.is_sorted_g())
      set_sorted_g();
  }

  /** \brief The move constructor. */
  OrdPopulation(OrdPopulation&& _p) noexcept: Base(std::move(_p)) {
    if(_p.is_sorted_g())
      set_sorted_g();
  }

  /* TODO locks */
  /** \brief Copy assignment of a compatible OrdPopulation.
   * \copydetails add(const Container&) */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  OrdPopulation& operator=(const OrdPopulation<CBase, ref_, Tag_, Pop_>& _p) {
    Base::operator=(std::move(_p));
    if(_p.is_sorted_g())
      set_sorted_g();
    return *this;
  }

  /** \brief Move assignment of a compatible OrdPopulation. */
  template<bool ref_, class Tag_, template<class, bool> class Pop_>
  OrdPopulation& operator=(OrdPopulation<CBase, ref_, Tag_, Pop_>&& _p) {
    Base::operator=(std::move(_p));
    if(_p.is_sorted_g())
      set_sorted_g();
    return *this;
  }

  void reserve(size_t count) {
    Base::reserve(count);
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
   * Applicable only if the type returned by `CBase::fitness()` allows total
   * ordering using `operator<`. This method generates a compile-time error in
   * specializations for which this condition is not satisfied. */
#ifdef DOXYGEN
  const Candidate<CBase>& best() {
#else
  template<class Ret = const Candidate<CBase>&>
  Ret NOINLINE best() {
#endif
    internal::read_lock lock(smp);
    if(this->empty())
      throw std::out_of_range("best(): BasePopulation is empty.");
    if(is_sorted_ug())
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
      throw std::out_of_range("rankSelect(): BasePopulation is empty.");
    ensure_sorted(lock);
    if(x == 1)
      return Base2::back();
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
      throw std::out_of_range("rankSelect(): BasePopulation is empty.");
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
    ensure_sorted(lock);
    return operator[](_rankSelect_dist(rng));
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
    bool was_sorted = is_sorted_ug();
    if(!lock.upgrade_if([newSize,this]() -> bool { return size() > newSize; }))
      return;
    ensure_sorted(lock);
    auto& dummy = Base2::front();  // needed by resize() if size() < newSize which can't happen
    Base2::resize(newSize, dummy); // (otherwise we would need a default constructor)
    if(was_sorted)
      set_sorted_ug();
  }

private:

  bool is_sorted_ug() const {
    return smp.get_mod_cnt() == last_sort_mod;
  }

  bool is_sorted_g() const {
    internal::read_lock lock(smp);
    return is_sorted_ug();
  }

  void set_sorted_ug() {
    last_sort_mod = smp.get_mod_cnt();
  }

  void set_sorted_g() {
    internal::write_lock lock(smp);
    set_sorted_ug();
  }

  void ensure_sorted(internal::rw_lock& lock) {
    if(!is_sorted_ug()) {
      internal::upgrade_lock up(lock);
      ++last_sort_mod;
      if(is_sorted_ug())
        return;
      std::sort(Base2::begin(), Base2::end());
      set_sorted_ug();
    }
  }

}; // class OrdPopulation

} // namespace gen
