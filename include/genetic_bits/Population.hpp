namespace gen {

template<class, bool>
class Population;

namespace internal {

template<class CBase, bool is_ref>
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
>::type; // type alias PopulationChooser

} // namespace internal


/** \brief The main Population template, representing a collection of
 * <b>Candidate</b>s. This is the entry point of most applications of the
 * framework.
 *
 * Depending on the properties of the candidate base class \b CBase,
 * \link gen::Population Population<CBase> \endlink becomes a synonyme for:
 * - OrdPopulation if \b CBase supports a <b>bool operator<()</b>,
 * - FloatPopulation if, moreover, its Fitness type is a simple floating type,
 * - DomPopulation if \b CBase supports a <b>bool operator<<()</b>,
 * - BasePopulation in all other cases (this is the subclass of the above
 *   three).
 *
 * Note that the existence of \b operator<() is checked first, so if both of
 * the operators are defined, \b operator<<() is ignored and the methods made
 * accessible in DomPopulation are not defined.
 *
 * A Population can be used as a container of <b>Candidate</b>s with read-only
 * access. (See Candidate for discussion about the relation between a
 * Candidate and \b CBase.) The functions \link BasePopulation::begin()
 * begin() \endlink and \link BasePopulation::end() end() \endlink are
 * exposed, returning random access iterators dereferencing to <b>const \link
 * Candidate Candidate<CBase>\endlink&</b> and allowing the iteration patterns
 * ```
 * for(auto& c : pop) { ... }
 * ```
 * and
 * ```
 * for(auto c : pop) { ... }
 * ```
 * Also, read-only element accessors \link BasePopulation::at() at() \endlink
 * and \link BasePopulation::operator[]() operator[]() \endlink are available
 * for direct access to candidates.
 *
 * See the documentation of the above classes for more detailed description.
 *
 * \tparam CBase the base class of the member candidates of this population.
 * See Candidate for details.
 * \tparam is_ref if set to \b true, this is a reference population. See
 * BasePopulation::Ref for more details.
 *
 * \see NSGAPopulation */
template<class CBase, bool is_ref = false>
class Population: public internal::PopulationChooser<CBase, is_ref> {

  using Base = internal::PopulationChooser<CBase, is_ref>;

public:

  using Base::Base;

}; // class Population

} // namespace gen
