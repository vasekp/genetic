namespace gen {

/** \brief Read-locks a Population for the duration given by the definition
 * scope.
 *
 * Use a PopulationLock wherever unguarded direct read access to a population
 * is needed to prevent other threads for modifying it during such section.
 * This could apply to, for example:
 * - copying, copy-assigning, or adding it to another population as a whole,
 * - using the population's iterators \link BasePopulation::begin()
 *   begin()\endlink and \link BasePopulation::end() end()\endlink,
 * - range-based for-loops.
 *
 * This mechanism allows to extend the thread-safe internal locking mechanism
 * to read-only external access to the population's members, most importantly
 * enumeration:
 * ```cpp
 * {
 *   gen::PopulationLock lock(pop); // scope starts here
 *   for(auto& candidate : pop) {
 *     DoSomething(candidate);
 *   }
 * } // scope ends, lock automatically released
 * ```
 * Another important use case is using the container algorithms of the
 * Standard Template Library with a Population. For example,
 * ```cpp
 * {
 *   gen::PopulationLock lock(pop); // scope starts here
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
 *
 * **Important:** Most internal Population functions use their own read locks.
 * It is forbidden to call them while holding a PopulationLock as this will
 * lead to undefined behaviour (likely a dead-lock). The following functions
 * are guaranteed to be exempt from this rule and can be called:
 * - BasePopulation::begin()
 * - BasePopulation::end()
 * - BasePopulation::begin() const
 * - BasePopulation::end() const
 * - BasePopulation::operator[]()
 * - BasePopulation::size() const
 *
 * Note that no PopulationLock is needed if any member functions are called
 * outside a parallel region. */
class PopulationLock {

  internal::read_lock lock;

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

}; // class PopulationLock

} // namespace gen
