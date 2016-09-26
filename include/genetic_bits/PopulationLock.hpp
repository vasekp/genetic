namespace gen {

/* Forward declarations */
template<class CBase, bool is_ref, class Tag,
  template<class, bool> class Population>
class BasePopulation;
/* End forward declarations */


/** \brief Read-locks a Population for the duration given by the definition
 * scope.
 *
 * Use a PopulationLock wherever unguarded direct read access to a population
 * is needed to prevent other threads for modifying it during such section.
 * (A PopulationLock is not needed if any member functions are called outside
 * a parallel region.) This could apply to, for example:
 * - copying, copy-assigning, or adding it to another population as a whole,
 * - using the population's iterators \link BasePopulation::begin()
 *   begin()\endlink, \link BasePopulation::end() end()\endlink or
 *   iterator-returning functions like \link OrdPopulation::rankSelect_i()
 *   rankSelect_i()
 *   \endlink,
 * - range-based for-loops.
 *
 * This mechanism allows to extend the thread-safe internal locking mechanism
 * to read-only external access to the population's members, most importantly
 * enumeration:
 * ```cpp
 * {
 *   gen::PopulationLock lock{pop}; // scope starts here
 *   for(auto& candidate : pop) {
 *     DoSomething(candidate);
 *   }
 * } // scope ends, lock automatically released
 * ```
 * Another important use case is using the container algorithms of the
 * Standard Template Library with a Population. For example,
 * ```cpp
 * {
 *   gen::PopulationLock lock{pop}; // scope starts here
 *   auto iterator = std::find_if(  // STL function
 *     pop.begin(),
 *     pop.end(),
 *     AnyPredicate // (const gen::Candidate<CBase>&) -> bool
 *   );
 *   if(iterator != pop.end()) {
 *     // iterator guaranteed to be valid even in a multithreaded scenario
 *     DoSomething(*iterator); // (const gen::Candidate<CBase>&) -> void
 *   }
 * } // scope ends, lock automatically released
 * ```
 * The iterators returned by Population do not allow to modify its contents so
 * only functions which do not attempt that can be used.
 *
 * **Important:** Most internal Population functions use their own read locks.
 * It is forbidden to call them while holding a PopulationLock as this will
 * lead to undefined behaviour (likely a dead-lock). The following functions
 * are guaranteed to be exempt from this rule and can be called:
 * - BasePopulation::begin()
 * - BasePopulation::end()
 * - BasePopulation::operator[]()
 * - BasePopulation::size() const
 * - any function with a name ending with \b _i.
 *
 * The invariance of the population under a PopulationLock may be broken by an
 * intentional call to a routine which modifies it in a compatible manner.
 * Currently this can be \link OrdPopulation::sort(PopulationLock&)
 * OrdPopulation::sort() \endlink or OrdPopulation::rankSelect_i() (which
 * calls the former). These functions reorder the members of a population.
 * They guarantee that no iterators are invalidated, but the candidates the
 * individual iterators point to after their call can be different than prior
 * to it. */
class PopulationLock {

  internal::read_lock lock;

#ifndef DOXYGEN
  /* for access to get() */
  template<class CBase, bool is_ref, class Tag,
    template<class, bool> class Population>
  friend class BasePopulation;
#endif

public:

  /** \brief Locks a given Population for read-only access.
   *
   * The lock is held as long as this object is in scope and released at
   * destruction, which can happen by leaving the syntactic block where it was
   * introduced or by throwing an exception.
   *
   * No other locking operations may be called in the invoking thread while
   * this object is being held. If a PopulationLock is requested while a write
   * operation is underway, the calling thread blocks until all modifications
   * are finished and an invariant population can be guaranteed.  However, if
   * any concurrent thread requests a write operation after a PopulationLock
   * is acquired, it is paused until the lock is released. */
  template<class CBase, bool is_ref, class Tag,
    template<class, bool> class Population>
  PopulationLock(BasePopulation<CBase, is_ref, Tag, Population>& pop):
    lock(pop.smp) { }

protected:

#ifndef DOXYGEN
  internal::read_lock& get() {
    return lock;
  }
#endif

}; // class PopulationLock

} // namespace gen
