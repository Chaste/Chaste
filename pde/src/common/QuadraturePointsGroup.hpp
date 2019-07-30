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

#ifndef QUADRATUREPOINTSGROUP_HPP_
#define QUADRATUREPOINTSGROUP_HPP_

#include "UblasCustomFunctions.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "GaussianQuadratureRule.hpp"
#include "LinearBasisFunction.hpp"
#include <vector>

/**
 * A simple class which takes in a mesh and a quadrature rule, and collects
 * all the quadrature points in physical space (rather than in natural element coordinates)
 * together in one data structure, for access.
 *
 * In a distributed mesh, a set of quad points can still be accessed for a given
 * element but then any missing data is marked with DOUBLE_UNSET
 */
template<unsigned DIM>
class QuadraturePointsGroup
{
private:

    /** The quadrature points in physical space */
    std::vector<c_vector<double,DIM> > data;

    /** Number of elements in given mesh */
    unsigned mNumElements;

    /** Number of quad points per element in given rule */
    unsigned mNumQuadPointsPerElement;

public:

    /**
     * Constructor takes in a mesh and a rule and computes and stores all
     * the quad points in physical space.
     *
     * @param rMesh
     * @param rQuadRule
     */
    QuadraturePointsGroup(AbstractTetrahedralMesh<DIM,DIM>& rMesh,
                          GaussianQuadratureRule<DIM>& rQuadRule);

    /**
     * @return a stored quadrature point by element index and quad index in the element.
     *
     * @param elementIndex
     * @param quadIndex
     */
    c_vector<double,DIM>& rGet(unsigned elementIndex, unsigned quadIndex);

    /**
     * @return the i-th stored quadrature point.
     *
     * @param i
     */
    c_vector<double,DIM>& rGet(unsigned i);

    /** @return number of elements in the mesh that was given in the constructor */
    unsigned GetNumElements() const;

    /** @return number of quad points per element in the rule that was given in the constructor */
    unsigned GetNumQuadPointsPerElement() const;

    /** @return total size, i.e. total number of quad points, i.e. num_elem times num_quad_points_per_elem */
    unsigned Size() const;
};

#endif /*QUADRATUREPOINTSGROUP_HPP_*/
