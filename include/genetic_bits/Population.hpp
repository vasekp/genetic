namespace gen {

namespace internal {

template<class CBase, class Tag = internal::empty, bool is_ref = false>
using PopulationChooser = typename std::conditional<
  Candidate<CBase>::Traits::is_comparable,
  typename std::conditional<
    Candidate<CBase>::Traits::is_float,
    FloatPopulation<CBase, Tag, is_ref>,
    OrdPopulation<CBase, Tag, is_ref>
  >::type,
  typename std::conditional<
    Candidate<CBase>::Traits::is_dominable,
    DomPopulation<CBase, Tag, is_ref>,
    BasePopulation<CBase, Tag, is_ref>
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
template<class CBase, class Tag = internal::empty, bool is_ref = false>
using Population = typename internal::PopulationChooser<CBase, Tag, is_ref>;

} // namespace gen
