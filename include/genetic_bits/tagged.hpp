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
class Tagged :
public std::conditional<ref,
           std::reference_wrapper<const Candidate>,
           Candidate
         >::type,
private Tag {

  typedef typename std::conditional<ref,
    std::reference_wrapper<const Candidate>,
    Candidate
  >::type Base;

public:
  Tagged(const Candidate& _c): Base(_c) { }
  Tagged(Candidate&& _c): Base(std::move(_c)) { }
  Tag& tag() { return static_cast<Tag&>(*this); }
};

} // namespace internal

} // namespace gen
