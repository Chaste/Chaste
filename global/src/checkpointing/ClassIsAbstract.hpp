/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef CLASSISABSTRACT_HPP_
#define CLASSISABSTRACT_HPP_

/** @file
 *
 * This file defines 4 macros to assist with explicitly declaring
 * to the Boost Serialization library when a class is abstract.
 * The interface for doing this changed in Boost 1.36.0, hence this
 * wrapper.
 *
 * The easy case is for a non-templated class.  For example, if you
 * have a class called AbstractClass, use
 *     CLASS_IS_ABSTRACT(AbstractClass)
 *
 * For classes templated over either 1 or 2 unsigned parameters, there
 * are helper macros TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED and
 * TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED.  For example, with a class
 *     template<unsigned SPACE_DIM, unsigned ELEMENT_DIM>
 *     class AbstractTemplatedClass { ... };
 * use
 *     TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractTemplatedClass)
 *
 * For a general templated class, you have to do a little extra work.
 * For example, with a class
 *     template<class C, unsigned U>
 *     class AbstractTemplatedClass { ... };
 * use
 *     namespace boost {
 *     namespace serialization {
 *     template<class C, unsigned U>
 *     struct is_abstract<AbstractTemplatedClass<C, U> >
 *         TEMPLATED_CLASS_IS_ABSTRACT_DEFN
 *     template<class C, unsigned U>
 *     struct is_abstract<const AbstractTemplatedClass<C, U> >
 *         TEMPLATED_CLASS_IS_ABSTRACT_DEFN
 *     }}
 */

#include <boost/version.hpp>

#if BOOST_VERSION >= 103600

// In Boost since 1.36.0, we need to use assume_abstract...
#include <boost/serialization/assume_abstract.hpp>

/**
 * Explicitly mark a non-templated class as being abstract
 * (Boost 1.36 and later).
 * @param T  the class
 */
#define CLASS_IS_ABSTRACT(T) BOOST_SERIALIZATION_ASSUME_ABSTRACT(T)

/**
 * Content of the is_abstract type to mark a templated class as abstract
 * (Boost 1.36 and later).
 */
#define TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    : boost::true_type {};

#else

// In Boost before 1.36.0, we use is_abstract...
#include <boost/serialization/is_abstract.hpp>

/**
 * Explicitly mark a non-templated class as being abstract
 * (Boost 1.35 and earlier).
 * @param T  the class
 */
#define CLASS_IS_ABSTRACT(T) BOOST_IS_ABSTRACT(T)

/**
 * Content of the is_abstract type to mark a templated class as abstract
 * (Boost 1.35 and earlier).
 */
#define TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    { \
        typedef mpl::bool_<true> type; \
        BOOST_STATIC_CONSTANT(bool, value = true); \
    };

#endif // BOOST_VERSION >= 103600

/**
 * Convenience macro to declare a class templated over a single unsigned
 * as abstract.
 * @param T  the class
 */
#define TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(T) \
    namespace boost { \
    namespace serialization { \
    template<unsigned U> \
    struct is_abstract<T<U> > \
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    template<unsigned U> \
    struct is_abstract<const T<U> > \
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    }}

/**
 * Convenience macro to declare a class templated over 2 unsigneds
 * as abstract.
 * @param T  the class
 */
#define TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(T) \
    namespace boost { \
    namespace serialization { \
    template<unsigned U1, unsigned U2> \
    struct is_abstract<T<U1, U2> > \
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    template<unsigned U1, unsigned U2> \
    struct is_abstract<const T<U1, U2> > \
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    }}

/**
 * Convenience macro to declare a class templated over 3 unsigneds
 * as abstract.
 * @param T  the class
 */
#define TEMPLATED_CLASS_IS_ABSTRACT_3_UNSIGNED(T) \
    namespace boost { \
    namespace serialization { \
    template<unsigned U1, unsigned U2, unsigned U3> \
    struct is_abstract<T<U1, U2, U3> > \
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    template<unsigned U1, unsigned U2, unsigned U3> \
    struct is_abstract<const T<U1, U2, U3> > \
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN \
    }}

#endif /*CLASSISABSTRACT_HPP_*/
