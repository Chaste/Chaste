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
 * Check for versions of Boost we don't support at all.
 */
#include <boost/version.hpp>
#if BOOST_VERSION < 103800
#error "Chaste doesn't support versions of Boost earlier than 1.38."
#elif BOOST_VERSION == 104100
// There's a bug in 1.41 with shared_ptr support; see e.g. http://sourceforge.net/apps/trac/easystroke/ticket/21
#error "Chaste won't work with Boost 1.41 due to a bug in its serialization library."
#endif

/**
 * \todo #2417 Mac OS X doesn't appear to be able to checkpoint dynamically loaded
 * models, so that functionality is switched off here.
 */
#ifndef CHASTE_CAN_CHECKPOINT_DLLS
#ifndef __APPLE__
#define CHASTE_CAN_CHECKPOINT_DLLS
#endif //Not APPLE
#endif

/**
 * Add some missing #includes for boost 1.56 and boost 1.57,
 * these are to work around bugs in these boost versions.
 *
 * See tickets #2585 (boost 1.56) and #2626 (boost 1.57).
 */
#if BOOST_VERSION == 105600
#include <boost/serialization/singleton.hpp>
#include <boost/serialization/extended_type_info.hpp>
#endif

#if BOOST_VERSION == 105700
#include <boost/serialization/type_info_implementation.hpp>
#endif

#endif // CHASTESERIALIZATION_HPP_
