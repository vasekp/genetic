namespace gen {

namespace internal {

/* Default if the user does not need any Tag */
struct empty { };


/* Wrapper for non-class Tags */
template<typename Tag>
struct TagWrap {
  Tag t;
  TagWrap(): t() { }
  TagWrap(Tag& _t): t(_t) { }
  operator Tag&() { return t; }
};

/* This specialization is an empty class, as a base it will take no extra
 * memory. (This would not hold if there was a named member like in the
 * above.) Also, tag() function is missing so the compilation fails if anyone
 * requests the tag. */
template<>
struct TagWrap<empty>: empty {
  TagWrap(...) { }
};


/* The base class for CandidateTagged. */
template<class CBase, bool ref>
using CTBase = typename std::conditional<ref,
                 std::reference_wrapper<const gen::Candidate<CBase>>,
                 gen::Candidate<CBase>
               >::type;


/* This is the element type of the std::vector underlying gen::Population. If
 * ref is false, it holds a Candidate and a Tag. If ref is true, it holds a
 * reference to a Candidate and an (independent) Tag. The user should never
 * encounter this mechanism directly as it transparently casts to a const
 * Candidate&. */
template<class CBase, typename Tag, bool ref>
struct CandidateTagged: public CTBase<CBase, ref>, private TagWrap<Tag> {
  CandidateTagged(const gen::Candidate<CBase>& _c): CTBase<CBase, ref>(_c) { }

  CandidateTagged(gen::Candidate<CBase>&& _c): CTBase<CBase, ref>(std::move(_c)) { }

  Tag& tag() { return static_cast<Tag&>(static_cast<TagWrap<Tag>&>(*this)); }
};

} // namespace internal

} // namespace gen
