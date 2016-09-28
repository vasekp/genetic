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
      std::is_convertible<
        decltype(std::declval<C>() < std::declval<C>()), bool
      >::value, bool
    >::type { return true; }

  template<typename>
  constexpr bool comparable(...) { return false; }


  /* Helper for enabling functions dependent on dominance */

  template<typename C>
  constexpr auto dominable(int) ->
    typename std::enable_if<
      std::is_convertible<
        decltype(std::declval<C>() << std::declval<C>()), bool
      >::value, bool
    >::type { return true; }

  template<typename>
  constexpr bool dominable(...) { return false; }


  /* Helper for Population::add(Container&) and Population::add(Container&&) */

  template<class C>
  constexpr auto is_container(int) ->
    typename std::enable_if<
      std::is_same<
        decltype(std::declval<C>().begin()),
        decltype(std::declval<C>().end())
      >::value, bool
    >::type { return true; }

  template<class>
  constexpr bool is_container(...) { return false; }


  /* Helper for functions taking an optional random number generator:
   * not to be confused with another overload */

  template<typename C>
  constexpr auto is_URNG(int) ->
    typename std::enable_if<
      std::is_integral<typename C::result_type>::value,
      bool
    >::type { return true; }

  template<class>
  constexpr bool is_URNG(...) { return false; }


} // namespace internal

} // namespace gen
