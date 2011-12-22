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
