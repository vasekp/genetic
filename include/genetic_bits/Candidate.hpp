namespace gen {

/** \brief The Candidate template.
 *
 * Takes a typename `Fitness` which can be a simple type or a class with a
 * a default constructor and, optionally, `operator<`. If the latter is
 * provided it is assumed that it defines a total ordering on the `Fitness`
 * type.
 *
 * The virtual implementation only keeps track of whether fitness has
 * been computed, and provides the precomputed value when available. The inner
 * function to compute fitness when required, computeFitness(), needs to be
 * provided in every derived class. Also, derived classes need to implement a
 * default constructor (or provide default values for all arguments in some
 * constructor) for the purposes of Population. */
template<typename Fitness>
class Candidate {
  mutable Fitness _fitness{};
  mutable bool fitnessValid = false;

  public:
  /** \brief Returns this candidate's fitness, calculating it on request if not
   * known from before. */
  Fitness fitness() const {
    if(!fitnessValid) {
      _fitness = computeFitness();
      fitnessValid = true;
    }
    return _fitness;
  }

  /** \brief Compares two `Candidate`s by the `Fitness`'s `operator<`.
   *
   * The `operator<` of `Fitness` is expected to be a total ordering for the
   * purposes of sorting the Population by fitness. This method is not
   * generated if `Fitness::operator<()` does not exist or does not return a
   * boolean value. */
  friend inline bool
  operator< (const Candidate<Fitness>& c1, const Candidate<Fitness>& c2) {
    return c1.fitness() < c2.fitness();
  }

  /** \brief Compares two `Candidate`s by the `Fitness`'s `operator<<`.
   *
   * The `operator<<` of `Fitness` is expected to be a partial ordering
   * denoting dominance in multiobjective searches. This method is not
   * generated if `Fitness::operator<<()` does not exist or does not return a
   * boolean value. */
  friend inline bool
  operator<< (const Candidate<Fitness>& c1, const Candidate<Fitness>& c2) {
    return c1.fitness() << c2.fitness();
  }

  virtual ~Candidate() { }

  protected:
  /** \brief The internal fitness computation, called the first time this
   * candidate's fitness() is queried. Every derived class must implement
   * this routine.
   *
   * In rare occasions it may happen that two threads ask to calculate a
   * Candidate's fitness simultaneously. This function must be designed
   * such that this causes no problem and must unconditionally return the
   * same value for the same input. */
  virtual Fitness computeFitness() const = 0;
}; // class Candidate

} // namespace gen
