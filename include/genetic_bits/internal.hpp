namespace gen {

/* Internal functions only to be used from other sources in this directory. */
namespace internal {

  /* Helpers for Population::rankSelect<std::exp> */
  template<double (*)(double)>
    struct is_exp: std::false_type { };

  template<>
    struct is_exp<std::exp>: std::true_type { };

  template<double (*f)(double)>
  double eval_in_product(double x, double y) {
    return f(x * y);
  }

  /* Helper for enabling functions dependent on total ordering */
  template<typename C>
  constexpr auto comparable(int) ->
    typename std::enable_if<
      std::is_convertible<decltype(std::declval<C>() < std::declval<C>()), bool>::value,
    bool>::type { return true; }

  template<typename C>
  constexpr bool comparable(...) { return false; }

  /* Helper for enabling functions dependent on dominance */
  template<typename C>
  constexpr auto dominable(int) ->
    typename std::enable_if<
      std::is_convertible<decltype(std::declval<C>() << std::declval<C>()), bool>::value,
    bool>::type { return true; }

  template<typename C>
  constexpr bool dominable(...) { return false; }
  
  /* Helpers for detecting if a Candidate is derived from ICandidate */
  template<class T>
  constexpr bool hasFT(typename T::_FitnessType*) { return true; }

  template<class T>
  constexpr bool hasFT(...) { return false; }

  template<class T>
  typename T::_FitnessType detectFT(T*);

  template<class T>
  void detectFT(...);

} // namespace internal

} // namespace gen
