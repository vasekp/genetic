namespace gen {

/** \brief The FloatPopulation template, adding functionality dependent on the
 * candidates' fitness being convertible to a simple floating point type.
 * \copydetails gen::BasePopulation */
template<class CBase, bool is_ref, class Tag, template<class, bool> class Population>
class FloatPopulation : public OrdPopulation<CBase, is_ref, Tag, Population> {

  typedef OrdPopulation<CBase, is_ref, Tag, Population> Base;

public:

  using Base::Base;
  using Base::smp;
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::operator[];

  /** \brief Creates an empty population. */
  FloatPopulation() = default;

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
  size_t _fitnessSelect_last_mod{(size_t)(~0)};
  double _fitnessSelect_last_bias{};

  template<class Ret, double (*fun)(double, double), class Rng>
  Ret NOINLINE fitnessSelect_int(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("fitnessSelect(): BasePopulation is empty.");
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

  /** \brief The return type of FloatPopulation::stat(). */
  struct Stat {
    double mean;  ///< The mean fitness of the BasePopulation.
    double stdev; ///< The standard deviation of fitness in the BasePopulation.
  };

  /** \brief Returns the mean fitness of the population and the standard
   * deviation.
   *
   * Applicable only to candidate classes whose fitness is a simple floating
   * point type or allows an implicit convertion to one. This method
   * generates a compile-time error in specializations for which this
   * condition is not satisfied.
   *
   * \see Stat */
  Stat stat() const {
    if(this->empty())
      throw std::out_of_range("stat(): BasePopulation is empty.");
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

}; // class FloatPopulation

} // namespace gen
