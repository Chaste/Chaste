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

#ifndef SMARTPOINTERS_HPP_
#define SMARTPOINTERS_HPP_

/**
 * @file
 * Includes the Boost shared_ptr smart pointer, and defines some useful macros to save
 * typing when using it.
 */

#include <boost/shared_ptr.hpp>

/**
 * Create a new instance of a class and assign it to a smart pointer.
 * @param ABS_TYPE  the type of the base of the class hierarchy
 * @param TYPE  the type of the concrete instance to create
 * @param NAME  the name of the pointer variable
 * @param ARGS  constructor arguments for the instance, in brackets
 */
#define MAKE_PTR_ABS(ABS_TYPE, TYPE, NAME, ARGS) boost::shared_ptr<ABS_TYPE> NAME(new TYPE ARGS)

/**
 * Create a new instance of a class and assign it to a smart pointer.
 * @param TYPE  the type of the concrete instance to create
 * @param NAME  the name of the pointer variable
 * @param ARGS  constructor arguments for the instance, in brackets
 */
#define MAKE_PTR_ARGS(TYPE, NAME, ARGS) MAKE_PTR_ABS(TYPE, TYPE, NAME, ARGS)

/**
 * Create a new instance of a class and assign it to a smart pointer.
 * @param TYPE  the type of the concrete instance to create
 * @param NAME  the name of the pointer variable
 */
#define MAKE_PTR(TYPE, NAME) MAKE_PTR_ABS(TYPE, TYPE, NAME, )

/**
 * Create a new class instance and reset a smart pointer to point at it.
 * @param NAME  the name of the pointer variable
 * @param TYPE  the type of the concrete instance to create
 * @param ARGS  constructor arguments for the instance, in brackets
 */
#define ASSIGN_PTR(NAME, TYPE, ARGS) NAME.reset(new TYPE ARGS)

#endif // SMARTPOINTERS_HPP_
