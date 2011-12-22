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

#ifndef CHASTESERIALIZATION_HPP_
#define CHASTESERIALIZATION_HPP_

/**
 * @file
 *
 * This header is a wrapper including some of the Boost serialization library
 * headers, along with a couple of standard C++ headers required to fix bugs
 * in Boost.
 *
 * Include this header in place of <boost/serialization/access.hpp>
 */

// Apparently 'new' (for boost's two phase construction) isn't included sometimes...
#include <new>
#include <climits> // See #1024.

#include <boost/serialization/access.hpp>

/**
 * Only Boost 1.37 and above can properly handle serialization of dynamically
 * loaded objects.  We define a convenience macro for code to test if this is
 * possible.
 */
#include <boost/version.hpp>
#ifndef CHASTE_CAN_CHECKPOINT_DLLS
#if BOOST_VERSION >= 103700
#define CHASTE_CAN_CHECKPOINT_DLLS
#endif
#endif

/**
 * Also check for versions of Boost we don't support at all.
 */
#if BOOST_VERSION < 103301
#error "Chaste doesn't support versions of Boost earlier than 1.33.1."
#elif BOOST_VERSION == 103500
// There's a bug in 1.35 which involves a
// #include <boost/serialization/extended_type_info_typeid.hpp>
// missing at the end of <boost/serialization/export.hpp>
// We might have worked around this now, but it's so old it's not worth checking!
#error "Chaste won't work with Boost 1.35 due to a bug in its serialization library."
#elif BOOST_VERSION == 104100
// There's another bug in 1.41 with shared_ptr support; see e.g.
// http://sourceforge.net/apps/trac/easystroke/ticket/21
#error "Chaste won't work with Boost 1.41 due to a bug in its serialization library."
#endif

#endif // CHASTESERIALIZATION_HPP_
