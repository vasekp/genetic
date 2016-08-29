/** \brief The `gen` namespace. */
namespace gen {

/** \brief The default random number generator for Population::rankSelect().
 *
 * A strong RNG is not a necessity for genetic applications, so speed was main
 * preference in choosing `minstd_rand`. Can be accessed freely by applications
 * of the framework. Similarly, alternative RNGs may be provided to
 * Population::rankSelect(). */
static thread_local std::minstd_rand rng{std::random_device()()};

} // namespace gen
