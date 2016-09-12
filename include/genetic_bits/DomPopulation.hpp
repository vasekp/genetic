namespace gen {

/** \brief The DomPopulation template, adding functionality dependent on
 * partial (dominance) ordering between candidates to a BasePopulation.
 * \copydetails gen::BasePopulation */
template<class CBase, bool is_ref, class Tag,
  template<class, bool> class Population>
class DomPopulation: public BasePopulation<CBase, is_ref, Tag, Population> {

  using Base = BasePopulation<CBase, is_ref, Tag, Population>;

public:

  using Base::Base;
  using Base::smp;
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::operator[];

  /** \brief Creates an empty population. */
  DomPopulation() = default;

  /** \brief Returns the number of candidates in this population dominated by
   * a given candidate. */
  friend size_t operator<< (const Candidate<CBase>& c,
      const DomPopulation& pop) {
    size_t cnt = 0;
    internal::read_lock lock{pop.smp};
    for(auto& cmp : pop)
      if(c << cmp)
        cnt++;
    return cnt;
  }

  /** \brief Returns the number of candidates in this population that
   * dominate a given candidate. */
  friend size_t operator<< (const DomPopulation& pop,
      const Candidate<CBase>& c) {
    size_t cnt = 0;
    internal::read_lock lock{pop.smp};
    for(auto& cmp : pop)
      if(cmp << c)
        cnt++;
    return cnt;
  }

  /** \brief Returns a nondominated subset of this population.
   *
   * The returned \link Ref \endlink remains valid until the original
   * population is modified.  Therefore there is a risk of invalidating it in
   * a multi-threaded program if another thread concurrently modifies the
   * population. If your code allows this, use front_v() instead.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
#ifdef DOXYGEN
  Ref front(bool parallel = true) const {
#else
  template<class Ret = typename Base::Ref>
  Ret NOINLINE front(bool parallel = true) const {
#endif
    internal::read_lock lock{smp};
    size_t sz = size();
    std::vector<char> dom(sz, 0);
    #pragma omp parallel for if(parallel) schedule(dynamic)
    for(size_t i = 0; i < sz; i++)
      for(size_t j = 0; j < sz; j++)
        if(!dom[j] && operator[](j) << operator[](i)) {
          dom[i] = 1;
          break; // no omp collapse
        }
    Ret ret{};
    for(size_t i = 0; i < sz; i++)
      if(!dom[i])
        ret.add(operator[](i));
    return ret;
  }

  /** \copybrief front()
   *
   * Works like front() but returns an independent population. */
  typename Base::Val NOINLINE front_v(bool parallel = true) const {
    return front<typename Base::Val>(parallel);
  }

}; // class DomPupolation

} // namespace gen
