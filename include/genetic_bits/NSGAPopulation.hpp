namespace gen {

/** \brief A variant of Population enabling nondominance sorting and related
 * selection methods.
 *
 * \tparam CBase the base class of the member candidates of this population.
 * Must implement a partial (dominance) order given by a <b>bool
 * %operator<<()</b>, otherwise a compile-time error is generated.  See
 * Candidate for further details.
 * \tparam is_ref if set to \b true, this is a reference population. See \link
 * NSGAPopulation::Ref Ref \endlink for more details. */
template<class CBase, bool is_ref = false>
class NSGAPopulation: public DomPopulation<CBase, is_ref, size_t, NSGAPopulation> {

  static_assert(Candidate<CBase>::Traits::is_dominable,
      "The fitness type of CBase needs to support bool operator<<()!");

  typedef DomPopulation<CBase, is_ref, size_t, NSGAPopulation> Base;
  typedef internal::PBase<CBase, is_ref, size_t> Base2;

  /* Protects: nsga_* */
  /* Promise: to be only acquired from within a read lock on the Base. */
  mutable internal::rw_semaphore nsga_smp{};

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

  /** \brief Creates an empty population. */
  NSGAPopulation() = default;

  /** \copydoc DomPopulation::front() */
#ifdef DOXYGEN
  Ref front(bool parallel = true) const {
#else
  template<class Ret = typename Base::Ref>
  Ret NOINLINE front(bool parallel = true) const {
#endif
    {
      internal::read_lock lock(smp);
      internal::read_lock nsga_lock(nsga_smp);
      if(smp.get_mod_cnt() == nsga_last_mod) {
        Ret ret{};
        for(auto& tg : (Base2&)(*this))
          if(tg.tag() == 0)
            ret.add(static_cast<Candidate<CBase>&>(tg));
        return ret;
      }
    }
    return Base::template front<Ret>(parallel);
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
    size_t& rank;       // reference to its tag in the original population
    size_t dom_cnt;     // number of candidates dominating this
    std::deque<_nsga_struct*> q;  // candidates dominated by this
  };

  void NOINLINE nsga_rate(bool parallel = false) {
    std::list<_nsga_struct> ref{};
    for(auto& tg : static_cast<Base2&>(*this))
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

  template<class Ret, double (*fun)(double, double), class Rng>
  Ret NOINLINE NSGASelect_int(double bias, Rng& rng) {
    internal::read_lock lock(smp);
    size_t sz = size();
    if(sz == 0)
      throw std::out_of_range("NSGASelect(): BasePopulation is empty.");
    internal::read_lock nsga_lock(nsga_smp);
    if(smp.get_mod_cnt() != nsga_last_mod || bias != nsga_last_bias) {
      nsga_lock.upgrade();
      if(smp.get_mod_cnt() != nsga_last_mod)
        nsga_rate(false);
      nsga_probs.clear();
      nsga_probs.reserve(sz);
      for(auto& tg : static_cast<Base2&>(*this))
        nsga_probs.push_back(1 / fun(tg.tag(), bias));
      nsga_dist = std::discrete_distribution<size_t>(nsga_probs.begin(), nsga_probs.end());
      nsga_last_mod = smp.get_mod_cnt();
      nsga_last_bias = bias;
    }
    size_t ix = nsga_dist(rng);
    return operator[](ix);
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
    internal::read_lock lock(smp);
    internal::read_lock nsga_lock(nsga_smp);
    if(smp.get_mod_cnt() != nsga_last_mod) {
      nsga_lock.upgrade();
      nsga_rate(parallel);
      nsga_last_mod = smp.get_mod_cnt();
    }
  }

}; // class NSGAPopulation

} // namespace gen
