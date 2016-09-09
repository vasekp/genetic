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
template<class CBase, bool ref, typename Tag>
struct CandidateTagged: public CTBase<CBase, ref>, private TagWrap<Tag> {

  using reference = const gen::Candidate<CBase>&;
  using rv_reference = CTBase<CBase, ref>&&;

  CandidateTagged(const gen::Candidate<CBase>& _c): CTBase<CBase, ref>(_c) { }

  CandidateTagged(CTBase<CBase, ref>&& _c): CTBase<CBase, ref>(std::move(_c)) { }

  Tag& tag() { return static_cast<Tag&>(static_cast<TagWrap<Tag>&>(*this)); }

  friend inline bool
  operator< (const CandidateTagged& c1, const CandidateTagged& c2) {
    return (const gen::Candidate<CBase>&)(c1) < (const gen::Candidate<CBase>&)(c2);
  }

  friend inline bool
  operator<< (const CandidateTagged& c1, const CandidateTagged& c2) {
    return (const gen::Candidate<CBase>&)(c1) << (const gen::Candidate<CBase>&)(c2);
  }
};



/* This class takes an iterator returning a CandidateTagged and modifies it so
 * that it returns const references to Candidate. (In case of move iterator,
 * Candidate&& or std::reference_wrapper<Candidate>&&.) */
template<class It, bool move = false>
class CTIterator: public It {
  public:
  typedef typename std::conditional<move,
    typename It::value_type::rv_reference,
    typename It::value_type::reference>::type Target;
  typedef typename It::iterator_category iterator_category;
  typedef typename It::difference_type difference_type;
  typedef Target value_type;
  typedef Target reference;
  typedef typename std::add_pointer<Target>::type pointer;

  CTIterator(It it): It(it) { };
  reference operator*() const { return static_cast<Target>(It::operator*()); }
  pointer operator->() const { return &operator*(); }
  CTIterator operator+(difference_type n) { return CTIterator(It::operator+(n)); }
  CTIterator operator-(difference_type n) { return CTIterator(It::operator-(n)); }
};


/* We need a specialization of std::move_iterator for our CTIterator. */
template<class It>
class move_iterator: public std::move_iterator<It> {
  public:
  explicit move_iterator(It it): std::move_iterator<It>(it) { }
};

template<class It>
class move_iterator<CTIterator<It, false>>:
public CTIterator<It, true> {
  public:
  explicit move_iterator(CTIterator<It> it):
      CTIterator<It, true>(it) { }
};

} // namespace internal

} // namespace gen
