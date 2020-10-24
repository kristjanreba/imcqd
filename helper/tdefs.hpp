#ifndef TDEFS_H
#define TDEFS_H
#include <type_traits>

template <class A> struct is_unique_pointer : std::false_type {};
template <class A, class B>
struct is_unique_pointer<std::unique_ptr<A, B>> : std::true_type {};

template <class A> struct is_const_unique_pointer : std::false_type {};
template <class A, class B>
struct is_const_unique_pointer<const std::unique_ptr<A, B>> : std::true_type {};

#endif
