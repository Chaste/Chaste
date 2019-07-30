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

#ifndef BOOSTFILESYSTEM_HPP_
#define BOOSTFILESYSTEM_HPP_

/**
 * @file
 * Include the Boost Filesystem library headers,
 * and set up the 'fs' namespace alias.
 * This header also ensures that we use version 2 of the library when available.
 */

#include <boost/version.hpp>

#if BOOST_VERSION <= 104900
/** Which version of the library to use. */
#define BOOST_FILESYSTEM_VERSION 2
/** How to get a leafname as a string (in version 2). */
#define PATH_LEAF_NAME(path) path.leaf()
#else
/** How to get a leafname as a string (in version 3). */
#define PATH_LEAF_NAME(path) path.leaf().string()
#endif

/*
 * There is potential binary incompatibility when linking the filesystem library (built with C++98) and Chaste (built
 * with C++11).  This is due to the use of scoped enums in the filesystem library, which have to be emulated pre-C++11
 * and thus generate different decorated function names.  See #2811 and https://svn.boost.org/trac10/ticket/6779
 * for details.
 * This was resolved in Boost 1.57, so we special-case the older versions.
 */
#if BOOST_VERSION <= 105600
/** Change boost SCOPED_ENUMS behaviour on older versions */
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#else
#include <boost/filesystem.hpp>
#endif
#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;

#endif // BOOSTFILESYSTEM_HPP_
