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

#ifndef CHASTESERIALIZATIONVERSION_HPP_
#define CHASTESERIALIZATIONVERSION_HPP_

/// gcov doesn't like this file...
#define COVERAGE_IGNORE

/**
 * @file
 * Provide a wrapper around Boost's serialization version to cope with changes in library
 * interface.
 * Include this header in place of <boost/serialization/version.hpp>
 *
 * For simple classes T, a version number N can be specified just by using the Boost macro
 * BOOST_CLASS_VERSION(T, N)
 *
 * However, templated classes need to expand the definition of this macro, the contents of
 * which changed in Boost 1.44.  Use the CHASTE_VERSION_CONTENT macro within your template.
 *
 * For example:
\code
namespace boost {
namespace serialization {
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM,  unsigned PROBLEM_DIM>
struct version<AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> >
{
    CHASTE_VERSION_CONTENT(1);
};
} // namespace serialization
} // namespace boost
\endcode
 *
 * @see http://www.boost.org/doc/libs/1_44_0/boost/serialization/version.hpp
 */

#include <boost/version.hpp>
#include <boost/serialization/version.hpp>

#if BOOST_VERSION >= 104400

/**
 * Content for the Boost serialization version template on Boost 1.44 and above.
 * @param N  the version number
 */
#define CHASTE_VERSION_CONTENT(N)                             \
    typedef boost::mpl::int_<N> type;                         \
    typedef boost::mpl::integral_c_tag tag;                   \
    BOOST_STATIC_CONSTANT(int, value = version::type::value)

#else // BOOST_VERSION >= 104400

/**
 * Content for the Boost serialization version template on Boost 1.43 and below.
 * @param N  the version number
 */
#define CHASTE_VERSION_CONTENT(N)                             \
    BOOST_STATIC_CONSTANT(unsigned, value = N)

#endif // BOOST_VERSION >= 104400


// gcov doesn't like this file...
#undef COVERAGE_IGNORE

#endif /*CHASTESERIALIZATIONVERSION_HPP_*/
