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

  template<typename>
  constexpr bool comparable(...) { return false; }

  /* Helper for enabling functions dependent on dominance */
  template<typename C>
  constexpr auto dominable(int) ->
    typename std::enable_if<
      std::is_convertible<decltype(std::declval<C>() << std::declval<C>()), bool>::value,
    bool>::type { return true; }

  template<typename>
  constexpr bool dominable(...) { return false; }
  
  /* Helper for detecting if a Candidate is a class derived from gen::Candidate or a
   * reference_wrapper of some */
  template<typename T>
  T unwrap(std::reference_wrapper<T>);

  template<typename T>
  T unwrap(T);

  template<typename C>
  decltype(unwrap(std::declval<C>()).fitness()) detectFT(C*);

  template<typename C>
  // Returning C is a trick; void would break Candidate<_FitnessType> in
  // the static_assert of Population.hpp. This works as well, if C is not a
  // Candidate<?> then C& is not convertible to const Candidate<C>& (or should
  // not be).
  C detectFT(...);

} // namespace internal

} // namespace gen
