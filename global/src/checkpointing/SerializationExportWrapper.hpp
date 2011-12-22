/*

Copyright (C) University of Oxford, 2005-2011

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

/// gcov doesn't like this file...
#define COVERAGE_IGNORE

/**
 * @file
 *
 * Defines some macros to register versions of templated classes with the
 * serialization library, for all space dimensions.  Also contains wrappers
 * around BOOST_CLASS_EXPORT and related functionality, which take care of
 * the differences introduced in new versions of Boost.
 *
 * In Boost 1.33.1 and 1.34, BOOST_CLASS_EXPORT should be placed in the .hpp
 * file for each class, and archive headers included only in tests, or special
 * 'archiver' class header files (e.g. CardiacSimulationArchiver.hpp).
 *
 * Serialization is broken in Boost 1.35.
 *
 * In Boost 1.36 (and up to 1.40) both the archive header includes and the
 * BOOST_CLASS_EXPORT should go in .cpp files.
 *
 * We don't yet support 1.41 and above, which introduce BOOST_CLASS_EXPORT_KEY
 * and BOOST_CLASS_EXPORT_IMPLEMENT.
 *
 * To handle all situations in Chaste:
 *   1. In .hpp files, include this header after the class definition.
 *   2. In .cpp files, after any other includes, include SerializationExportWrapperForCpp.hpp.
 * In both cases, CHASTE_CLASS_EXPORT should be used instead of BOOST_CLASS_EXPORT.
 *
 * There are also variant macros for common cases of templated classes:
 *  - EXPORT_TEMPLATE_CLASS_SAME_DIMS
 *  - EXPORT_TEMPLATE_CLASS_ALL_DIMS
 *  - EXPORT_TEMPLATE_CLASS1
 *  - EXPORT_TEMPLATE_CLASS2
 *  - EXPORT_TEMPLATE_CLASS3
 *
 * The latter 3 macros are usable in any situation where a class has up to 3 template parameters,
 * and you know what values will be needed.  Unfortunately a fully general template solution seems
 * to be impossible in any Boost version (the library makes use of either explicit instantiation or
 * singleton class instantiation).
 */

#include <boost/version.hpp>

////////////////////////////////////////////////////////////////////////////////
// Make sure includes happen in the correct place.  This has to go before
// the SERIALIZATIONEXPORTWRAPPER_HPP_ guard, since we need it to be seen
// by both .hpp and .cpp files.

#if BOOST_VERSION < 103600
// Boost 1.34 and older - export goes in headers
#ifndef CHASTE_SERIALIZATION_CPP
#include <boost/serialization/export.hpp>
#endif // CHASTE_SERIALIZATION_CPP

#elif BOOST_VERSION < 104100
// Boost 1.36-1.40 - export goes in .cpp, along with archive includes
#ifdef CHASTE_SERIALIZATION_CPP
#include "CheckpointArchiveTypes.hpp"
#include <boost/serialization/export.hpp>
#endif // CHASTE_SERIALIZATION_CPP

#else
// Boost 1.41 and newer - export goes in both, with archive includes in .cpp
#include <boost/serialization/extended_type_info.hpp> // We get compile errors without this...
#include <boost/serialization/export.hpp>
#ifdef CHASTE_SERIALIZATION_CPP
#include "CheckpointArchiveTypes.hpp"
#endif // CHASTE_SERIALIZATION_CPP

#endif
// Done includes

////////////////////////////////////////////////////////////////////////////////
#if BOOST_VERSION >= 104100 && defined(CHASTE_SERIALIZATION_CPP)
// .cpp file needs to use BOOST_CLASS_EXPORT_IMPLEMENT, so we need
// to redefine the macros from the .hpp file.  Hence this can't go
// in the include guard.

#undef CHASTE_CLASS_EXPORT_TEMPLATED
/**
 * General export for templated classes.
 * @param T  a type
 * @param S  a unique string for the class + specific template parameter values
 */
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)    \
   BOOST_CLASS_EXPORT_IMPLEMENT(T)

#undef CHASTE_CLASS_EXPORT_INTERNAL
/**
 * What CHASTE_CLASS_EXPORT expands to when it isn't a no-op.
 * @param T  the class to export
 */
#define CHASTE_CLASS_EXPORT_INTERNAL(T)        \
   BOOST_CLASS_EXPORT_IMPLEMENT(T)

#endif // BOOST_VERSION >= 104100 && defined(CHASTE_SERIALIZATION_CPP)

////////////////////////////////////////////////////////////////////////////////

#ifndef SERIALIZATIONEXPORTWRAPPER_HPP_
#define SERIALIZATIONEXPORTWRAPPER_HPP_
// Code in the next block is only compiled when the .hpp is first seen

////////////////////////////////////////////////////////////////////////////////

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>

// Handle broken BOOST_CLASS_EXPORT in Boost 1.36 & 1.37
#if BOOST_VERSION >= 103600 && BOOST_VERSION < 103800

/**
 * Global unique identifier for Boost 1.36 and 1.37, to avoid using the source line
 * number as a 'unique' id.
 * @param T type this guid is for
 * @param K string version of T
 * @param S uniquifying string
 */
#define CHASTE_CLASS_EXPORT_GUID(T, K, S)                                           \
namespace                                                                           \
{                                                                                   \
    ::boost::archive::detail::guid_initializer< T > const &                         \
        BOOST_PP_CAT(BOOST_PP_CAT(boost_serialization_guid_initializer_, __LINE__), S)               \
        = ::boost::serialization::singleton<                                        \
            ::boost::archive::detail::guid_initializer< T >                         \
          >::get_mutable_instance().export_guid(K);                                 \
}

/**
 * General export for templated classes.
 * @param T  a type
 * @param S  a unique string for the class + specific template parameter values
 */
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)    \
    CHASTE_CLASS_EXPORT_GUID(                  \
        T,                                     \
        BOOST_PP_STRINGIZE(T), S               \
    )                                          \

/**
 * What CHASTE_CLASS_EXPORT expands to when it isn't a no-op.
 * @param T  the class to export
 */
#define CHASTE_CLASS_EXPORT_INTERNAL(T)        \
   CHASTE_CLASS_EXPORT_TEMPLATED(T, T)


// The interface changes yet again in Boost 1.41, and we need something in both .hpp and .cpp...
#elif BOOST_VERSION >= 104100
// .hpp file needs to use BOOST_CLASS_EXPORT_KEY

/**
 * General export for templated classes.
 * @param T  a type
 * @param S  a unique string for the class + specific template parameter values
 */
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)   \
   BOOST_CLASS_EXPORT_KEY(T)

/**
 * What CHASTE_CLASS_EXPORT expands to when it isn't a no-op.
 * @param T  the class to export
 */
#define CHASTE_CLASS_EXPORT_INTERNAL(T)       \
   BOOST_CLASS_EXPORT_KEY(T)

#else // BOOST_VERSION < 103600 || (BOOST_VERSION >= 103800 && BOOST_VERSION < 104100)
//Do exactly as we did before (so that archives created with 1.33 don't have to be re-generated)

/**
 * General export for templated classes.
 * @param T  a type
 * @param S  a unique string for the class + specific template parameter values
 */
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)    \
   BOOST_CLASS_EXPORT(T)

/**
 * What CHASTE_CLASS_EXPORT expands to when it isn't a no-op.
 * @param T  the class to export
 */
#define CHASTE_CLASS_EXPORT_INTERNAL(T)        \
   BOOST_CLASS_EXPORT(T)

#endif // BOOST_VERSION >= 103600 && BOOST_VERSION < 103800

template<class> struct pack;
/**
 * Argument pack for macros.  Used to give a type for templated classes when exporting.
 * See http://lists.boost.org/Archives/boost/2004/08/70149.php for more.
 *
 * \todo Check if we need this on Boost>=1.38.  Even if it's not needed there, we might still need it to load earlier archives.
 */
template<class T> struct pack<void (T)> {
    typedef T type; /**< Type definition. */
};

/**
 * Defines the export key for a class templated over 1 parameter.
 * @param CLASS  the class
 * @param P1  the template parameter
 */
#define CHASTE_EXPORT_KEY_1(CLASS, P1) \
    BOOST_PP_CAT(CLASS, P1)

/**
 * Defines the export key for a class templated over 2 parameters.
 * @param CLASS  the class
 * @param P1  the first template parameter
 * @param P2  the second template parameter
 */
#define CHASTE_EXPORT_KEY_2(CLASS, P1, P2) \
    BOOST_PP_CAT(BOOST_PP_CAT(CLASS, P1), P2)

/**
 * Defines the export key for a class templated over 3 parameters.
 * @param CLASS  the class
 * @param P1  the first template parameter
 * @param P2  the second template parameter
 * @param P3  the third template parameter
 */
#define CHASTE_EXPORT_KEY_3(CLASS, P1, P2, P3) \
    BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(CLASS, P1), P2), P3)

/**
 * Defines the export type for a class templated over 1 parameter.
 * @param CLASS  the class
 * @param P1  the template parameter
 */
#define CHASTE_PACK_1(CLASS, P1) pack<void (CLASS< P1 >)>::type

/**
 * Defines the export type for a class templated over 2 parameters.
 * @param CLASS  the class
 * @param P1  the first template parameter
 * @param P2  the second template parameter
 */
#define CHASTE_PACK_2(CLASS, P1, P2) pack<void (CLASS< P1,P2 >)>::type

/**
 * Defines the export type for a class templated over 3 parameters.
 * @param CLASS  the class
 * @param P1  the first template parameter
 * @param P2  the second template parameter
 * @param P3  the third template parameter
 */
#define CHASTE_PACK_3(CLASS, P1, P2, P3) pack<void (CLASS< P1,P2,P3 >)>::type

// End of include guard - code below here is executed in both .hpp and .cpp
#endif // SERIALIZATIONEXPORTWRAPPER_HPP_
////////////////////////////////////////////////////////////////////////////////

// Since CHASTE_CLASS_EXPORT_TEMPLATED and CHASTE_CLASS_EXPORT_INTERNAL are re-defined for the .cpp file
// in some Boost versions, the following macros will also need re-defining.

#ifdef EXPORT_TEMPLATE_CLASS3_INTERNAL
// Avoid re-definition when called from a .cpp file
#undef EXPORT_TEMPLATE_CLASS3_INTERNAL
#undef EXPORT_TEMPLATE_CLASS2_INTERNAL
#undef EXPORT_TEMPLATE_CLASS1_INTERNAL
#undef EXPORT_TEMPLATE_CLASS_ALL_DIMS_INTERNAL
#undef EXPORT_TEMPLATE_CLASS_SAME_DIMS_INTERNAL
#endif // EXPORT_TEMPLATE_CLASS3_INTERNAL

// Macros for templated classes

/**
 * Export a templated class with 3 parameters.
 * This is the definition of EXPORT_TEMPLATE_CLASS3 when it isn't a no-op.
 * @param CLASS the class (without parameters)
 * @param E  first template parameter
 * @param S  second template parameter
 * @param P  third template parameter
 */
#define EXPORT_TEMPLATE_CLASS3_INTERNAL(CLASS, E, S, P) \
    CHASTE_CLASS_EXPORT_TEMPLATED( CHASTE_PACK_3(CLASS, E, S, P), CHASTE_EXPORT_KEY_3(CLASS, E, S, P) )

/**
 * Export a templated class with 2 parameters.
 * This is the definition of EXPORT_TEMPLATE_CLASS2 when it isn't a no-op.
 * @param CLASS the class (without parameters)
 * @param E  first template parameter
 * @param S  second template parameter
 */
#define EXPORT_TEMPLATE_CLASS2_INTERNAL(CLASS, E, S) \
    CHASTE_CLASS_EXPORT_TEMPLATED( CHASTE_PACK_2(CLASS, E, S), CHASTE_EXPORT_KEY_2(CLASS, E, S) )

/**
 * Export a templated class with 1 parameter.
 * This is the definition of EXPORT_TEMPLATE_CLASS1 when it isn't a no-op.
 * @param CLASS the class (without parameters)
 * @param D  template parameter
 */
#define EXPORT_TEMPLATE_CLASS1_INTERNAL(CLASS, D) \
    CHASTE_CLASS_EXPORT_TEMPLATED( CHASTE_PACK_1(CLASS, D), CHASTE_EXPORT_KEY_1(CLASS, D) )

/**
 * Export a class templated over both element and space dimension, for all valid
 * combinations.
 * This is the definition of EXPORT_TEMPLATE_CLASS_ALL_DIMS when it isn't a no-op.
 * @param CLASS the class (without parameters)
 */
#define EXPORT_TEMPLATE_CLASS_ALL_DIMS_INTERNAL(CLASS) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 1) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 2) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 3) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 2, 2) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 2, 3) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 3, 3)

/**
 * Export a class templated over both element and space dimension, for combinations
 * where ELEMENT_DIM == SPACE_DIM.
 * This is the definition of EXPORT_TEMPLATE_CLASS_SAME_DIMS when it isn't a no-op.
 * @param CLASS the class (without parameters)
 */
#define EXPORT_TEMPLATE_CLASS_SAME_DIMS_INTERNAL(CLASS) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 1) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 2) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 3)

// Now define the macros that users actually call.
// Again this goes outside the include guard, so it is seen by both .hpp and .cpp files.

// However, we don't want to define things twice, so...
#if !defined(CHASTE_CLASS_EXPORT) || defined(CHASTE_SERIALIZATION_CPP)
#ifdef CHASTE_SERIALIZATION_CPP
// Remove the definitions from when we were included via an .hpp file
#undef CHASTE_CLASS_EXPORT
#undef EXPORT_TEMPLATE_CLASS_SAME_DIMS
#undef EXPORT_TEMPLATE_CLASS_ALL_DIMS
#undef EXPORT_TEMPLATE_CLASS1
#undef EXPORT_TEMPLATE_CLASS2
#undef EXPORT_TEMPLATE_CLASS3
#endif // CHASTE_SERIALIZATION_CPP

#if (BOOST_VERSION < 103600  && ! defined(CHASTE_SERIALIZATION_CPP)) || \
    (BOOST_VERSION >= 103600 && defined(CHASTE_SERIALIZATION_CPP)) || \
    (BOOST_VERSION >= 104100)
// Boost 1.34 and older - export goes in headers
// Boost 1.36 and newer - export goes in .cpp
// Boost 1.41 and newer - key goes in .hpp, implement goes in .cpp

/**
 * Define the serialization export key for this class.
 * @param T  the class
 */
#define CHASTE_CLASS_EXPORT(T)                 CHASTE_CLASS_EXPORT_INTERNAL(T)
/**
 * Export a class templated over both element and space dimension, for combinations
 * where ELEMENT_DIM == SPACE_DIM.
 * @param CLASS the class (without parameters)
 */
#define EXPORT_TEMPLATE_CLASS_SAME_DIMS(CLASS) EXPORT_TEMPLATE_CLASS_SAME_DIMS_INTERNAL(CLASS)
/**
 * Export a class templated over both element and space dimension, for all valid
 * combinations.
 * @param CLASS the class (without parameters)
 */
#define EXPORT_TEMPLATE_CLASS_ALL_DIMS(CLASS)  EXPORT_TEMPLATE_CLASS_ALL_DIMS_INTERNAL(CLASS)
/**
 * Export a templated class with 1 parameter.
 * @param CLASS the class (without parameters)
 * @param D  template parameter
 */
#define EXPORT_TEMPLATE_CLASS1(CLASS, D)       EXPORT_TEMPLATE_CLASS1_INTERNAL(CLASS, D)
/**
 * Export a templated class with 2 parameters.
 * @param CLASS the class (without parameters)
 * @param E  first template parameter
 * @param S  second template parameter
 */
#define EXPORT_TEMPLATE_CLASS2(CLASS, E, S)    EXPORT_TEMPLATE_CLASS2_INTERNAL(CLASS, E, S)
/**
 * Export a templated class with 3 parameters.
 * @param CLASS the class (without parameters)
 * @param E  first template parameter
 * @param S  second template parameter
 * @param P  third template parameter
 */
#define EXPORT_TEMPLATE_CLASS3(CLASS, E, S, P) EXPORT_TEMPLATE_CLASS3_INTERNAL(CLASS, E, S, P)

#else

/**
 * No-op for this Boost version in this setting.
 * @param T
 */
#define CHASTE_CLASS_EXPORT(T)
/**
 * No-op for this Boost version in this setting.
 * @param CLASS
 */
#define EXPORT_TEMPLATE_CLASS_SAME_DIMS(CLASS)
/**
 * No-op for this Boost version in this setting.
 * @param CLASS
 */
#define EXPORT_TEMPLATE_CLASS_ALL_DIMS(CLASS)
/**
 * No-op for this Boost version in this setting.
 * @param CLASS
 * @param D
 */
#define EXPORT_TEMPLATE_CLASS1(CLASS, D)
/**
 * No-op for this Boost version in this setting.
 * @param CLASS
 * @param E
 * @param S
 */
#define EXPORT_TEMPLATE_CLASS2(CLASS, E, S)
/**
 * No-op for this Boost version in this setting.
 * @param CLASS
 * @param E
 * @param S
 * @param P
 */
#define EXPORT_TEMPLATE_CLASS3(CLASS, E, S, P)

#endif // Long if!
#endif // !defined(CHASTE_CLASS_EXPORT) || defined(CHASTE_SERIALIZATION_CPP)

#undef COVERAGE_IGNORE
