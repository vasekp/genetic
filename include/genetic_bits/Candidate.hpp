namespace gen {

/** \brief The Candidate template.
 *
 * Takes a base class \b CBase which needs to implement a stateless function of
 * signature
 * ```
 * Fitness CBase::fitness() const
 * ```
 * Here, \b Fitness can be a simple type or a class, optionally supporting
 * <b>bool operator<()</b>, representing a total ordering, or <b>bool
 * operator<<()</b>, representing a partial ordering. If either of the
 * operators exist, they are extended to the Candidate type, comparing the
 * fitness of two candidates.  This is used to switch accessibility of various
 * Population methods.
 *
 * The value returned by Candidate::fitness() is the same as returned by its
 * <b>CBase::fitness()</b>, with the difference that the latter is called
 * exactly once, at the point of construction of this wrapper. This prevents
 * methods using the fitness from repeatedly requesting the computation which
 * is assumed to be expensive.
 *
 * All methods of \link Population Population<CBase>\endlink returning its
 * members or references to them return a \link Candidate<CBase>\endlink. This
 * can be implicitly converted to a \b CBase (or <b>const CBase&</b>), but
 * that loses the single-evaluation guarantee. It's recommended to use \b auto
 * or \b auto& in iterating over a Population or storing results of its
 * function calls. */
template<class CBase>
class Candidate: public CBase {

public:

#ifndef DOXYGEN
  struct Traits {

    using FitnessType = decltype(std::declval<CBase>().fitness());

    constexpr static bool is_comparable = internal::comparable<FitnessType>(0);

    constexpr static bool is_dominable = internal::dominable<FitnessType>(0);

    constexpr static bool is_float = std::is_convertible<FitnessType, double>::value;

  };
#endif

private:

  /* Making this const would mean deleting implicit copy and move assignments.
   * Let's suffice with it being private and all accessor functions returning
   * a const(&). */
  typename Traits::FitnessType _fitness;

public:

  /** \brief Copy-initialization from a \b CBase. */
  Candidate(const CBase& b): CBase(b), _fitness(CBase::fitness()) { }

  /** \brief Move-initialization from a \b CBase. */
  Candidate(CBase&& b): CBase(std::move(b)), _fitness(CBase::fitness()) { }

  /** \brief Returns this candidate's fitness, calculated at construction
   * time. */
  const typename Traits::FitnessType& fitness() const {
    return _fitness;
  }

  /** \brief Compares two <b>Candidate</b>s by the <b>Fitness</b>'s \b
   * %operator<().
   *
   * The \b %operator<() of \b Fitness is expected to be a total ordering for
   * the purposes of sorting the Population by fitness (used by OrdPopulation
   * and its derived classes). This method is not generated if
   * <b>Fitness::operator<()</b> does not exist or does not return a boolean
   * value. */
  friend inline bool
  operator< (const Candidate& c1, const Candidate& c2) {
    return c1.fitness() < c2.fitness();
  }

  /** \brief Compares two <b>Candidate</b>s by the <b>Fitness</b>'s \b
   * %operator<<().
   *
   * The \b %operator<<() of \b Fitness is expected to be a partial ordering
   * denoting dominance in multiobjective searches (used by DomPopulation
   * and its derived classes). This method is not generated if
   * <b>Fitness::operator<<()</b> does not exist or does not return a boolean
   * value. */
  friend inline bool
  operator<< (const Candidate& c1, const Candidate& c2) {
    return c1.fitness() << c2.fitness();
  }

}; // class Candidate

} // namespace gen
