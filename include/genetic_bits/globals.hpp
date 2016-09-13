/** \brief The main namespace of the framework.
 *
 * The entry point of most applications is either gen::Population or
 * gen::NSGAPopulation. Other <b>...%Population</b> classes are not intended
 * to be used directly. */
namespace gen {

/** \brief The default random number generator for this framework.
 *
 * A strong RNG is not a necessity for genetic applications, so speed was main
 * preference in choosing \b std::minstd_rand. Can be accessed freely by
 * applications of the framework. Similarly, alternative RNGs may be provided
 * to functions like \link BasePopulation::randomSelect()
 * randomSelect()\endlink. */
static thread_local std::minstd_rand rng{std::random_device{}()};

} // namespace gen
