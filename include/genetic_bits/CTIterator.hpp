namespace gen {

namespace internal {

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
