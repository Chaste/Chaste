/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CHASTESERIALIZATIONVERSION_HPP_
#define CHASTESERIALIZATIONVERSION_HPP_

/// gcov doesn't like this file...
// LCOV_EXCL_START

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
// LCOV_EXCL_STOP

#endif /*CHASTESERIALIZATIONVERSION_HPP_*/
