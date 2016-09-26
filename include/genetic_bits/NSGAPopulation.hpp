namespace gen {

/** \brief A variant of Population enabling nondominance sorting and related
 * selection methods.
 *
 * \tparam CBase the base class of the member candidates of this population.
 * The type returned by \b CBase::fitness() must implement a strict partial
 * order (\b a *dominates* \b b) given by a <b>bool %operator<<()</b>,
 * otherwise a compile-time error is generated.
 * \tparam is_ref if set to \b true, this is a reference population. See \link
 * NSGAPopulation::Ref Ref \endlink for more details. */
template<class CBase, bool is_ref = false>
class NSGAPopulation:
  public BasePopulation<CBase, is_ref, size_t, NSGAPopulation>,
  public internal::PopulationChooser<CBase, is_ref, size_t, NSGAPopulation>
{

  static_assert(Candidate<CBase>::Traits::is_dominable,
      "The fitness type of CBase needs to support bool operator<<()!");

  using Base = BasePopulation<CBase, is_ref, size_t, NSGAPopulation>;
  using Dom = DomPopulation<CBase, is_ref, size_t, NSGAPopulation>;

  std::discrete_distribution<size_t> nsga_dist{};
  std::vector<double> nsga_probs{};
  size_t nsga_last_mod{(size_t)(~0)};
  double nsga_last_bias{};

public:

  using Base::Base;
  using Base::smp;
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::operator[];
  using typename Base::iterator;

  /** \brief Creates an empty population. */
  NSGAPopulation() = default;

  /** \copydoc DomPopulation::front() */
#ifdef DOXYGEN
  Ref front(bool parallel = true) const {
#else
  template<class Ret = typename Base::Ref>
  NOINLINE Ret front(bool parallel = true) const {
#endif
    {
      internal::read_lock lock{smp};
      if(smp.get_mod_cnt() == nsga_last_mod) {
        Ret ret{};
        for(auto& tg : Base::as_vec())
          if(tg.tag() == 0)
            ret.add(static_cast<const Candidate<CBase>&>(tg));
        return ret;
      }
    }
    return static_cast<const Dom&>(*this).front<Ret>(parallel);
  }

  /** \copydoc DomPopulation::front_v() */
  typename Base::Val front_v(bool parallel = true) const {
    return front<typename Base::Val>(parallel);
  }

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
   * \tparam fun A \b constexpr pointer to a function of signature either
   * <b>double(*)(double)</b> or <b>double(*)(double, double)</b>. In the
   * former case, the argument is <b>x * bias</b>, in the latter case the
   * arguments are \b x and \b bias, where \b x denotes the dominance rank of
   * the candidate.  It must be positive and strictly increasing in \b x for
   * <b>bias > 0</b>.  This function will be built in at compile time,
   * eliminating a function pointer lookup. The default is \b std::exp,
   * another usual choice is \b std::pow.
   * \param bias > 0 determines how much low-dominated solutions are preferred.
   * Zero would mean no account on dominance rank in the selection process
   * whatsoever. The bigger the value the more low-dominated candidates are
   * likely to be selected.
   * \param rng the random number generator, or gen::rng by default.
   *
   * \returns a constant reference to a randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& NSGASelect(double bias, Rng& rng = rng);

  /** \copybrief NSGASelect()
   *
   * Works like NSGASelect() but returns by value.
   *
   * \returns a copy of the randomly chosen candidate.
   *
   * \throws std::out_of_bounds if called on an empty population. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> NSGASelect_v(double bias, Rng& rng = rng);

  /** \copybrief NSGASelect()
   *
   * Works like NSGASelect() but returns an iterator.
   *
   * This function relies on a read lock acquired externally for the
   * population via a PopulationLock. This lock will guard the validity of the
   * returned iterator.
   *
   * \returns an iterator pointing to the randomly selected candidate, end()
   * if the population is empty. */
  template<double (*fun)(...) = std::exp, class Rng = decltype(rng)>
  iterator NSGASelect_i(PopulationLock& lock, double bias, Rng& rng = rng);

#else

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  const Candidate<CBase>& NSGASelect(double mult, Rng& rng = rng) {
    internal::read_lock lock{smp};
    return *NSGASelect_int<
      &internal::eval_in_product<fun>
    >(mult, rng, lock, true);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  const Candidate<CBase>& NSGASelect(double bias, Rng& rng = rng) {
    internal::read_lock lock{smp};
    return *NSGASelect_int<fun>(bias, rng, lock, true);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  Candidate<CBase> NSGASelect_v(double bias, Rng& rng = rng) {
    internal::read_lock lock{smp};
    return *NSGASelect_int<
      &internal::eval_in_product<fun>
    >(bias, rng, lock, true);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  Candidate<CBase> NSGASelect_v(double bias, Rng& rng = rng) {
    internal::read_lock lock{smp};
    return *NSGASelect_int<fun>(bias, rng, lock, true);
  }

  template<double (*fun)(double) = std::exp, class Rng = decltype(rng)>
  iterator NSGASelect_i(PopulationLock& lock, double bias, Rng& rng = rng) {
    return NSGASelect_int<
      &internal::eval_in_product<fun>
    >(bias, rng, lock, false);
  }

  template<double (*fun)(double, double), class Rng = decltype(rng)>
  iterator NSGASelect_i(PopulationLock& lock, double bias, Rng& rng = rng) {
    return NSGASelect_int<fun>(bias, rng, lock, false);
  }

#endif

private:

  struct nsga_struct {
    internal::CandidateTagged<CBase, is_ref, size_t>& rCT;
    size_t dom_cnt;                 // number of candidates dominating this
    std::vector<nsga_struct*> dom;  // candidates dominated by this

    nsga_struct(internal::CandidateTagged<CBase, is_ref, size_t>& ref):
      rCT(ref), dom_cnt(0), dom{} { };

    const Candidate<CBase>& rCand() {
      return static_cast<const Candidate<CBase>&>(rCT);
    }

    size_t& rRank() {
      return rCT.tag();
    }
  };

  NOINLINE void nsga_rate(bool parallel = false) {
    std::vector<nsga_struct> vec{};
    std::vector<nsga_struct*> ref{};
    size_t sz = size();
    /* vec: size() preallocated nsga_structs, each pointing to one Candidate
     * ref: size() pointers to vec, originally 1:1 in order
     * Other ways tried:
     * - std::unique_ptr<nsga_struct>: unnecessary individual allocation
     * - list embedded in a vector: slow due to special cases
     * - validity flag: slow */
    vec.reserve(sz);
    ref.reserve(sz);
    for(auto& tg : Base::as_vec()) {
      vec.push_back({tg});
      ref.push_back(&vec.back());
    }
    /* Calculate dominance counts and lists of dominated candidates */
    #pragma omp parallel for if(parallel)
    for(size_t i = 0; i < sz; i++) {
      nsga_struct& op1 = vec[i];
      for(nsga_struct& op2 : vec)
        if(op1.rCand() << op2.rCand()) {
          op1.dom.push_back(&op2);
          #pragma omp atomic
          op2.dom_cnt++;
        }
    }
    /* ref contains the candidates with yet unassigned rank. When it becomes
     * empty, we're done. */
    size_t cur_rank = 0;
    std::vector<nsga_struct> cur_front{};
    auto first = ref.begin();
    auto last = ref.end();
    while(first != last) {
      cur_front.clear();
      for(auto it = first; it != last; )
        if((*it)->dom_cnt == 0) {
          cur_front.push_back(std::move(**it));
          std::swap(*it, *(--last));
        } else
          it++;
      /* Break arrows from members of cur_front and assign final rank to
       * them */
      for(auto& r : cur_front) {
        for(auto pdom : r.dom)
          pdom->dom_cnt--;
        r.rRank() = cur_rank;
      }
      cur_rank++;
    }
  }

  template<double (*fun)(double, double), class Rng>
  NOINLINE iterator NSGASelect_int(double bias, Rng& rng,
      internal::rw_lock& lock, bool validate) {
    size_t sz = size();
    if(sz == 0) {
      if(validate)
        throw std::out_of_range("NSGASelect(): Population is empty.");
      else
        return end();
    }
    if(smp.get_mod_cnt() != nsga_last_mod || bias != nsga_last_bias) {
      lock.upgrade(false); // and keep upgraded
      if(smp.get_mod_cnt() != nsga_last_mod)
        nsga_rate(false);
      nsga_probs.clear();
      nsga_probs.reserve(sz);
      for(auto& tg : Base::as_vec())
        nsga_probs.push_back(1 / fun(tg.tag(), bias));
      nsga_dist = std::discrete_distribution<size_t>
        (nsga_probs.begin(), nsga_probs.end());
      nsga_last_mod = smp.get_mod_cnt();
      nsga_last_bias = bias;
    }
    return begin() + nsga_dist(rng);
  }

public:

  /** \brief Updates the metainformation about the NSGA rank of each
   * candidate.
   *
   * This happens by default whenever NSGASelect() is called after the
   * population has been modified. However, since the sorting is a rather
   * expensive operation, in a multithreaded setting this would mean all other
   * threads have to wait for the calling thread to finish the sorting before
   * proceeding. This method can be called before the work is split between
   * threads, and as such, can benefit from parallelization itself.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
  void precompute(bool parallel = true) {
    internal::read_lock lock{smp};
    if(smp.get_mod_cnt() != nsga_last_mod) {
      lock.upgrade(false);
      nsga_rate(parallel);
      nsga_last_mod = smp.get_mod_cnt();
    }
  }

}; // class NSGAPopulation

} // namespace gen
