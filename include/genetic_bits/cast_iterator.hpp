namespace gen {

namespace internal {

/* This class takes an iterator and a class C and modifies the iterator so
 * that it returns const references to C. This is used to hide the fact that
 * the underlying class is different from C. */
template<class C, class It>
class cast_iterator: public It {
  public:
  typedef typename It::iterator_category iterator_category;
  typedef typename It::difference_type difference_type;
  typedef C value_type;
  typedef C reference;
  typedef typename std::add_pointer<C>::type pointer;

  cast_iterator(It it): It(it) { };
  reference operator*() const { return static_cast<C>(It::operator*()); }
  pointer operator->() const { return &operator*(); }
  cast_iterator operator+(difference_type n) { return cast_iterator(It::operator+(n)); }
  cast_iterator operator-(difference_type n) { return cast_iterator(It::operator-(n)); }
};


/* We need a specialization of std::move_iterator for our cast_iterator. */
template<class It>
class move_iterator: public std::move_iterator<It> {
  public:
  explicit move_iterator(It it): std::move_iterator<It>(it) { }
};

/* The move iterator needs to return a Candidate&& or a
 * std::reference_wrapper<Candidate>&&. */
template<class C, class It>
class move_iterator<cast_iterator<C, It>>:
public cast_iterator<typename It::value_type::rv_reference, It> {
  public:
  explicit move_iterator(cast_iterator<C, It> it):
      cast_iterator<typename It::value_type::rv_reference, It>(it) { }
};

} // namespace internal

} // namespace gen
