namespace gen {

/** \brief The FloatPopulation template, adding functionality dependent on the
 * candidates' fitness being convertible to a simple floating point type.
 * \copydetails gen::BasePopulation */
template<class CBase, bool is_ref, class Tag,
  template<class, bool> class Population>
class FloatPopulation: public OrdPopulation<CBase, is_ref, Tag, Population> {

  using Base = BasePopulation<CBase, is_ref, Tag, Population>;
  using Derived = Population<CBase, is_ref>;

  Base& base() {
    return static_cast<Base&>(static_cast<Derived&>(*this));
  }

  const Base& base() const {
    return static_cast<const Base&>(static_cast<const Derived&>(*this));
  }

  using iterator = typename Base::iterator;

  std::discrete_distribution<size_t> fitnessSelect_dist{};
  std::vector<double> fitnessSelect_probs{};
  size_t fitnessSelect_last_mod{(size_t)(~0)};
  double fitnessSelect_last_bias{};

public:

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
   * \tparam fun A \b constexpr pointer to a function of signature either
   * <b>double(*)(double)</b> or <b>double(*)(double, double)</b>. In the
   * former case, the argument is <b>x * bias</b>, in the latter case the
   * arguments are \b x and \b bias, where \b x denotes the fitness of each
   * candidate.  It must be positive and strictly increasing in \b x for
   * <b>bias > 0</b>.  This function will be built in at compile time,
   * eliminating a function pointer lookup. The default is [**std::exp**]
   * (http://en.cppreference.com/w/cpp/numeric/math/exp), another usual choice
   * is [**std::pow**] (http://en.cppreference.com/w/cpp/numeric/math/pow).
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
  const Candidate<CBase>& fitnessSelect(double bias, Rng& rng = rng);

  /** \copybrief fitnessSelect()
   *
   * Works like fitnessSelect() but returns by value.
   *
   * \returns a copy of the randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  iterator fitnessSelect_v(double bias, Rng& rng = rng);

  /** \copybrief fitnessSelect()
   *
   * Works like fitnessSelect() but returns an iterator.
   *
   * This function relies on a read lock acquired externally for the
   * population via a PopulationLock. This lock will guard the validity of the
   * returned iterator.
   *
   * \returns an iterator pointing to the randomly selected candidate, end()
   * if the population is empty. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  iterator fitnessSelect_i(PopulationLock& lock, double bias, Rng& rng = rng);

#else

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& fitnessSelect(double mult, Rng& rng = rng) {
    internal::read_lock lock{base().smp};
    return *fitnessSelect_int<
      &internal::eval_in_product<fun>
    >(mult, rng, lock, true);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate<CBase>& fitnessSelect(double bias, Rng& rng = rng) {
    internal::read_lock lock{base().smp};
    return *fitnessSelect_int<fun>(bias, rng, lock, true);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> fitnessSelect_v(double bias, Rng& rng = rng) {
    internal::read_lock lock{base().smp};
    return *fitnessSelect_int<
      &internal::eval_in_product<fun>
    >(bias, rng, lock, true);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate<CBase> fitnessSelect_v(double bias, Rng& rng = rng) {
    internal::read_lock lock{base().smp};
    return *fitnessSelect_int<fun>(bias, rng, lock, true);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  iterator fitnessSelect_i(PopulationLock& lock, double bias, Rng& rng = rng) {
    return fitnessSelect_int<
      &internal::eval_in_product<fun>
    >(bias, rng, lock.get(), false);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  iterator fitnessSelect_i(PopulationLock& lock, double bias, Rng& rng = rng) {
    return fitnessSelect_int<fun>(bias, rng, lock.get(), false);
  }

#endif

private:

  template<double (*fun)(double, double), class Rng>
  NOINLINE iterator fitnessSelect_int(double bias, Rng& rng,
      internal::rw_lock& lock, bool validate) {
    size_t sz = base().size();
    if(sz == 0) {
      if(validate)
        throw std::out_of_range("fitnessSelect(): Population is empty.");
      else
        return base().end();
    }
    if(fitnessSelect_last_mod != base().smp.get_mod_cnt()
       || bias != fitnessSelect_last_bias) {
      lock.upgrade(false); // and keep upgraded
      fitnessSelect_probs.clear();
      fitnessSelect_probs.reserve(sz);
      for(auto& c : base())
        fitnessSelect_probs.push_back(1 / fun(c.fitness(), bias));
      fitnessSelect_dist = std::discrete_distribution<size_t>
        (fitnessSelect_probs.begin(), fitnessSelect_probs.end());
      fitnessSelect_last_mod = base().smp.get_mod_cnt();
      fitnessSelect_last_bias = bias;
    }
    return base().begin() + fitnessSelect_dist(rng);
  }

public:

  /** \brief The return type of FloatPopulation::stat(). */
  struct Stat {
    double mean;  ///< The mean fitness of the population.
    double stdev; ///< The standard deviation of fitness in the population.
  }; // struct FloatPopulation<>::Stat

  /** \brief Returns the mean fitness of the population and the standard
   * deviation.
   *
   * \see Stat
   *
   * \throws std::out_of_bounds if called on an empty population. */
  NOINLINE Stat stat() const {
    if(base().empty())
      throw std::out_of_range("stat(): Population is empty.");
    double f, sf = 0, sf2 = 0;
    internal::read_lock lock{base().smp};
    for(auto& c : base()) {
      f = c.fitness();
      sf += f;
      sf2 += f*f;
    }
    size_t sz = base().size();
    double dev2 = sf2/sz - sf/sz*sf/sz;
    return {sf/sz,
            dev2 >= 0 ? sqrt(dev2) : 0};
  }

}; // class FloatPopulation<CBase, is_ref, Tag, Population>

} // namespace gen
