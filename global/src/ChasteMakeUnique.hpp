#ifndef CHASTEMAKEUNIQUE_HPP_
#define CHASTEMAKEUNIQUE_HPP_

/** @file
 *
 * This file is not part of Chaste, per se.
 *
 * This file contains an implementation of make_unique adopted by the ISO C++ committee, which is now
 * in the C++14 standard.  The function was not added to the C++11 standard, but there are good reasons for
 * using it in place of using raw `new` (see, for example, Effective Modern C++ (Scott Meyers) Item 21).
 * We, therefore, include this implementation, which can be dropped once Chaste adopts C++14 proper.
 *
 * Details of the original file with this implementation is as follows:
 *
 * Document number: N3656
 * Date: 2013-04-18
 * Project: Programming Language C++, Library Working Group
 * Reply-to: Stephan T. Lavavej <stl@microsoft.com>
 *
 * and can be found at
 * https://isocpp.org/files/papers/N3656.txt
 *
 * Changes from the published version:
 *  * Changed namespace from `std` to `our` to avoid potential ambiguity
 *  * Added `std::` qualifier to `unique_ptr`, `forward`, and `remove_extent`
 *
 */

#include <memory>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace our {
    template<class T> struct _Unique_if {
        typedef std::unique_ptr<T> _Single_object;
    };

    template<class T> struct _Unique_if<T[]> {
        typedef std::unique_ptr<T[]> _Unknown_bound;
    };

    template<class T, size_t N> struct _Unique_if<T[N]> {
        typedef void _Known_bound;
    };

    template<class T, class... Args>
    typename _Unique_if<T>::_Single_object
    make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    template<class T>
    typename _Unique_if<T>::_Unknown_bound
    make_unique(size_t n) {
        typedef typename std::remove_extent<T>::type U;
        return std::unique_ptr<T>(new U[n]());
    }

    template<class T, class... Args>
    typename _Unique_if<T>::_Known_bound
    make_unique(Args&&...) = delete;
}

#endif // CHASTEMAKEUNIQUE_HPP_
