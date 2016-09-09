namespace gen {

/** \brief The DomPopulation template, adding functionality dependent on
 * partial (dominance) ordering between candidates to a BasePopulation.
 * \copydetails gen::BasePopulation */
template<class CBase, bool is_ref, class Tag, template<class, bool> class Population>
class DomPopulation : public BasePopulation<CBase, is_ref, Tag, Population> {

  typedef BasePopulation<CBase, is_ref, Tag, Population> Base;
  typedef internal::PBase<CBase, is_ref, Tag> Base2;

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
  friend size_t operator<< (const Candidate<CBase>& c, const DomPopulation& pop) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(auto& cmp : pop)
      if(c << cmp)
        cnt++;
    return cnt;
  }

  /** \brief Returns the number of candidates in this population that
   * dominate a given candidate. */
  friend size_t operator<< (const DomPopulation& pop, const Candidate<CBase>& c) {
    size_t cnt = 0;
    internal::read_lock lock(pop.smp);
    for(auto& cmp : pop)
      if(cmp << c)
        cnt++;
    return cnt;
  }

  /** \brief Returns a nondominated subset of this population.
   *
   * The returned Ref remains valid until the original population is modified.
   * Therefore there is a risk of invalidating it in a multi-threaded program
   * if another thread concurrently modifies the population. If your code
   * allows this, use front_v() instead.
   *
   * \param parallel controls parallelization using OpenMP (on by default) */
#ifdef DOXYGEN
  Ref front(bool parallel = true) const {
#else
  template<class Ret = typename Base::Ref>
  Ret NOINLINE front(bool parallel = true) const {
#endif
    Ret ret{};
    internal::read_lock lock(smp);
    size_t sz = size();
    std::vector<char> dom(sz, 0);
    #pragma omp parallel for if(parallel)
    for(size_t i = 0; i < sz; i++) {
      for(size_t j = 0; j < sz; j++)
        if(!dom[j] && operator[](j) << operator[](i)) {
          dom[i] = 1;
          break;
        }
      if(!dom[i])
        ret.add(operator[](i)); // thread-safe
    }
    return ret;
  }

  /** \copybrief front()
   *
   * Works like front() but returns an independent DomPopulation. */
  typename Base::Val NOINLINE front_v(bool parallel = true) const {
    return front<typename Base::Val>(parallel);
  }

}; // class DomPupolation

} // namespace gen
