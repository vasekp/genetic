namespace gen {

namespace internal {

/* This class takes an iterator and a class C and modifies the iterator so
 * that it returns const references to C. This is used to hide the fact that
 * the underlying class is different from C. */
template<class C, class It>
class cast_iterator: public It {
  public:
  cast_iterator(It&& it): It(it) { };
  const C& operator*() const { return static_cast<C&>(It::operator*()); }
};


/* We need a specialization of std::move_iterator for our cast_iterator. */
template<class It>
class move_iterator: public std::move_iterator<It> {
  public:
  explicit move_iterator(It it): std::move_iterator<It>(it) { }
};

template<class C, class It>
class move_iterator<cast_iterator<C, It>>:
public cast_iterator<C&&, std::move_iterator<It>> {
  public:
  explicit move_iterator(cast_iterator<C, It> it):
    cast_iterator<C&&, std::move_iterator<It>>(std::move_iterator<It>(it)) { }
};

} // namespace internal

} // namespace gen
