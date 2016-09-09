namespace gen {

template<class CBase, bool is_ref = false>
  class Population;

namespace internal {

template<class CBase, bool is_ref = false>
using PopulationChooser = typename std::conditional<
  Candidate<CBase>::Traits::is_comparable,
  typename std::conditional<
    Candidate<CBase>::Traits::is_float,
    FloatPopulation<CBase, is_ref, internal::empty, gen::Population>,
    OrdPopulation<CBase, is_ref, internal::empty, gen::Population>
  >::type,
  typename std::conditional<
    Candidate<CBase>::Traits::is_dominable,
    DomPopulation<CBase, is_ref, internal::empty, gen::Population>,
    BasePopulation<CBase, is_ref, internal::empty, gen::Population>
  >::type
>::type;

} // namespace internal


/** \brief The main Population template, representing a collection of \link
 * Candidate Candidates \endlink. This is the entry point of most applications
 * of the framework.
 *
 * Depending on the properties of the candidate base class \b CBase,
 * \link gen::Population Population<CBase> \endlink becomes a synonyme for:
 * - OrdPopulation if \b CBase supports a <b>bool operator<()</b>,
 * - FloatPopulation if, moreover, its Fitness type is a simple floating type,
 * - DomPopulation if \b CBase supports a <b>bool operator<<()</b>,
 * - BasePopulation in all other cases (this is the subclass of the above
 *   three).
 *
 * See the documentation of the above classes for more detailed description. */
template<class CBase, bool is_ref>
class Population : public internal::PopulationChooser<CBase, is_ref> {

  typedef internal::PopulationChooser<CBase, is_ref> Base;

  public:

  using Base::Base;

};

} // namespace gen
