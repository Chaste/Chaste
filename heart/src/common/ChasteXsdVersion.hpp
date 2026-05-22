/*

Copyright (c) 2005-2026, University of Oxford.
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


#ifndef CHASTEXSDVERSION_HPP_
#define CHASTEXSDVERSION_HPP_

#include <xsd/cxx/version.hxx>

#include "Exception.hpp"

/*
 * Compatibility layer for CodeSynthesis XSD version macros.
 *
 * XSD 4.2+ provides:
 *
 *   LIBXSD_VERSION_MAJOR
 *   LIBXSD_VERSION_MINOR
 *   LIBXSD_VERSION_PATCH
 *
 * Older XSD versions provide:
 *
 *   XSD_INT_VERSION
 *
 * where XSD_INT_VERSION has format AABBCCDD:
 *
 *   AA - major version
 *   BB - minor version
 *   CC - patch version
 *   DD - alpha/beta version, or 00 for a final release
 */

#if defined(LIBXSD_VERSION_MAJOR) && \
    defined(LIBXSD_VERSION_MINOR) && \
    defined(LIBXSD_VERSION_PATCH)

    #define CHASTE_XSD_VERSION_MAJOR LIBXSD_VERSION_MAJOR
    #define CHASTE_XSD_VERSION_MINOR LIBXSD_VERSION_MINOR
    #define CHASTE_XSD_VERSION_PATCH LIBXSD_VERSION_PATCH

#elif defined(XSD_INT_VERSION)

    /*
     * Pre-release old-style versions have DD != 00. This is currently unsupported here because:
     *  1. it's very unlikely indeed that someone is using a pre-release XSD < 4.2
     *  2. the old AABBCCDD encoding decrements AABBCC for alpha/beta versions, making
     *     major/minor/patch extraction very fiddly.
     */
    #if (XSD_INT_VERSION % 100L) != 0L
        EXCEPTION("Unsupported pre-release CodeSynthesis XSD version");
    #endif

    #define CHASTE_XSD_VERSION_MAJOR \
        ((XSD_INT_VERSION / 1000000L) % 100L)

    #define CHASTE_XSD_VERSION_MINOR \
        ((XSD_INT_VERSION / 10000L) % 100L)

    #define CHASTE_XSD_VERSION_PATCH \
        ((XSD_INT_VERSION / 100L) % 100L)

#else

    EXCEPTION("Unsupported pre-release CodeSynthesis XSD version");

#endif

#define CHASTE_XSD_VERSION_AT_LEAST(MAJOR, MINOR, PATCH)              \
    ((CHASTE_XSD_VERSION_MAJOR > (MAJOR)) ||                          \
     (CHASTE_XSD_VERSION_MAJOR == (MAJOR) &&                          \
      CHASTE_XSD_VERSION_MINOR > (MINOR)) ||                          \
     (CHASTE_XSD_VERSION_MAJOR == (MAJOR) &&                          \
      CHASTE_XSD_VERSION_MINOR == (MINOR) &&                          \
      CHASTE_XSD_VERSION_PATCH >= (PATCH)))


#endif /*CHASTEXSDVERSION_HPP_*/
