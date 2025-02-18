/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef MESHUTILITYFUNCTIONS_HPP_
#define MESHUTILITYFUNCTIONS_HPP_

#include "UblasVectorInclude.hpp"

#include <cfloat>

/**
 * @file
 * This file contains utility functions for mesh operations.
 */


/**
 * Evenly space a number of points along a given path so that the distance (along the path) between each two consecutive
 * output points is identical.
 *
 * @tparam DIM the space dimension of the path vertices
 * @param path a vector of vertices defining the existing path along which to distribute points
 * @param closedPath whether the path forms a closed polygon
 * @param permuteOrder whether to randomly permute the input vertices, if the path is closed
 * @param numPointsToPlace the number of points to distribute evenly along the path
 * @param targetSpacing the target distance (along the path) between output points. If set, this overrides numPointsToPlace.
 * @return
 */
template <std::size_t DIM>
std::vector<c_vector<double, DIM>> EvenlySpaceAlongPath(
        std::vector<c_vector<double, DIM>>& path,
        bool closedPath,
        bool permuteOrder,
        std::size_t numPointsToPlace,
        double targetSpacing = DBL_MAX
);

#endif /*MESHUTILITYFUNCTIONS_HPP_*/
