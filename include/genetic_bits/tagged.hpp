namespace gen {

namespace internal {

/* Default if the user does not need any Tag */
struct empty { };


/* This is the element type of the std::vector underlying gen::Population. If
 * ref is false, it holds a Candidate and a Tag. If ref is true, it holds a
 * reference to a Candidate and an (independent) Tag. The user should never
 * encounter this mechanism directly as it transparently casts to a const
 * Candidate&. */
template<class Candidate, class Tag, bool ref>
class Tagged {
  typename std::conditional<ref,
             std::reference_wrapper<const Candidate>,
             Candidate
           >::type c;
  Tag t{};

  public:
  Tagged(const Candidate& _c): c(_c) { }
  operator const Candidate&() const { return c; }
  Tag& tag() { return t; }
};

} // namespace internal

} // namespace gen
