namespace gen {

/** \brief The Candidate template.
 *
 * Takes a base class `CBase` which needs to implement a stateless function of
 * signature
 * ```
 * Fitness CBase::fitness() const
 * ```
 * Here, `Fitness` can be a simple type or a class, optionally supporting
 * `bool operator<()`, representing a total ordering, or `bool operator<<()`,
 * representing a partial ordering. If either of the operators exist, they are
 * extended to the Candidate type, comparing the fitness of two candidates.
 * This is used to switch accessibility of various Population methods.
 *
 * The value returned by Candidate::fitness() is the same as returned by its
 * `CBase::fitness()`, with the difference that the latter is called exactly
 * once, at the point of construction of this wrapper. This prevents methods
 * using the fitness from repeatedly requesting the computation which is
 * assumed to be expensive.
 *
 * All methods of `Population<CBase>` returning its members or references to
 * them return a `Candidate<CBase>`. This can be implicitly converted to a
 * `CBase` (or `const CBase&`), but that loses the single-evaluation
 * guarantee. It's recommended to use `auto` or `auto&` in iterating over a
 * Population or storing results of its function calls. */
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

  /** \brief Copy-initialization from a `CBase`. */
  Candidate(const CBase& b): CBase(b), _fitness(CBase::fitness()) { }

  /** \brief Move-initialization from a `CBase`. */
  Candidate(CBase&& b): CBase(std::move(b)), _fitness(CBase::fitness()) { }

  /** \brief Returns this candidate's fitness, calculated at construction
   * time. */
  const typename Traits::FitnessType& fitness() const {
    return _fitness;
  }

  /** \brief Compares two `Candidate`s by the `Fitness`'s `operator<`.
   *
   * The `operator<` of `Fitness` is expected to be a total ordering for the
   * purposes of sorting the Population by fitness. This method is not
   * generated if `Fitness::operator<()` does not exist or does not return a
   * boolean value. */
  friend inline bool
  operator< (const Candidate& c1, const Candidate& c2) {
    return c1.fitness() < c2.fitness();
  }

  /** \brief Compares two `Candidate`s by the `Fitness`'s `operator<<`.
   *
   * The `operator<<` of `Fitness` is expected to be a partial ordering
   * denoting dominance in multiobjective searches. This method is not
   * generated if `Fitness::operator<<()` does not exist or does not return a
   * boolean value. */
  friend inline bool
  operator<< (const Candidate& c1, const Candidate& c2) {
    return c1.fitness() << c2.fitness();
  }
}; // class Candidate

} // namespace gen
