namespace gen {

namespace internal {

/* Default if the user does not need any Tag */
struct empty { };


/* Wrapper for non-class Tags */
template<typename Tag>
struct TagWrap {

  Tag t;

  TagWrap(): t() { }

  TagWrap(Tag& t_): t(t_) { }

  operator Tag&() { return t; }

  operator const Tag&() const { return t; }

}; // struct TagWrap<Tag>


/* This specialization is an empty class, as a base it will take no extra
 * memory. (This would not hold if there was a named member like in the
 * above.) Also, operator Tag&() is missing so the compilation fails if anyone
 * requests the tag. */
template<>
struct TagWrap<empty> {

  TagWrap(...) { }

}; // struct TagWrap<empty>



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
template<class CBase, bool ref, typename Tag>
struct CandidateTagged : public CTBase<CBase, ref>, private TagWrap<Tag> {

  using reference = const gen::Candidate<CBase>&;
  using rv_reference = CTBase<CBase, ref>&&;

  CandidateTagged(const gen::Candidate<CBase>& c): CTBase<CBase, ref>(c) { }

  CandidateTagged(CTBase<CBase, ref>&& c): CTBase<CBase, ref>(std::move(c)) { }

  CandidateTagged& operator=(const gen::Candidate<CBase>& c) {
    CTBase<CBase, ref>::operator=(c);
    return *this;
  }

  CandidateTagged& operator=(CTBase<CBase, ref>&& c) {
    CTBase<CBase, ref>::operator=(std::move(c));
    return *this;
  }

  Tag& tag() {
    return static_cast<Tag&>(static_cast<TagWrap<Tag>&>(*this));
  }

  const Tag& tag() const {
    return static_cast<const Tag&>(static_cast<const TagWrap<Tag>&>(*this));
  }

  friend inline bool
  operator< (const CandidateTagged& c1, const CandidateTagged& c2) {
    return static_cast<reference>(c1) < static_cast<reference>(c2);
  }

  friend inline bool
  operator<< (const CandidateTagged& c1, const CandidateTagged& c2) {
    return static_cast<reference>(c1) << static_cast<reference>(c2);
  }

}; // struct CandidateTagged<Candidate, ref, Tag>



/* This class takes an iterator returning a CandidateTagged and modifies it so
 * that it returns const references to Candidate. (In case of move iterator,
 * Candidate&& or std::reference_wrapper<Candidate>&&.) */
template<class It, bool move = false>
class CTIterator: public It {

public:

  using Target = typename std::conditional<move,
    typename It::value_type::rv_reference,
    typename It::value_type::reference>::type;

  using iterator_category = typename It::iterator_category;
  using difference_type   = typename It::difference_type;
  using value_type        = Target;
  using reference         = Target;
  using pointer           = typename std::add_pointer<Target>::type;

  CTIterator(It it): It(it) { };

  reference operator*() const {
    return static_cast<Target>(It::operator*());
  }

  pointer operator->() const {
    return &operator*();
  }

  CTIterator operator+(difference_type n) {
    return {It::operator+(n)};
  }

  CTIterator operator-(difference_type n) {
    return {It::operator-(n)};
  }

}; // class CTIterator<It, move>


/* We need a specialization of std::move_iterator and std::reverse_iterator
 * for our CTIterator. */
template<class It>
class move_iterator: public std::move_iterator<It> {

public:

  explicit move_iterator(It it): std::move_iterator<It>(it) { }

}; // class move_iterator<It>


template<class It>
class move_iterator<CTIterator<It, false>>: public CTIterator<It, true> {

public:

  explicit move_iterator(CTIterator<It> it): CTIterator<It, true>(it) { }

}; // class move_iterator<CTIterator>

} // namespace internal

} // namespace gen


namespace std {

template<class It, bool move>
class reverse_iterator<gen::internal::CTIterator<It, move>> :
  public gen::internal::CTIterator<reverse_iterator<It>, move>
{

  using Base = gen::internal::CTIterator<reverse_iterator<It>, move>;

public:

  constexpr reverse_iterator(reverse_iterator<It> x): Base(x) { }

}; // class std::reverse_iterator<CTIterator>

} // namespace std
