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

#ifndef QUADRATUREPOINTSGROUP_HPP_
#define QUADRATUREPOINTSGROUP_HPP_

#include "UblasCustomFunctions.hpp"
#include "TetrahedralMesh.hpp"
#include "GaussianQuadratureRule.hpp"
#include "LinearBasisFunction.hpp"
#include <vector>

/**
 * A simple class which takes in a mesh and a quadrature rule, and collects
 * are the quadrature points (in physical space ie several for each element)
 * together in one data structure, for access.
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
    QuadraturePointsGroup(TetrahedralMesh<DIM,DIM>& rMesh,
                          GaussianQuadratureRule<DIM>& rQuadRule);

    /**
     * Access the stored quadrature point by element index and quad index in the element.
     *
     * \todo this method should be renamed rGet() as it returns a reference
     *
     * @param elementIndex
     * @param quadIndex
     */
    c_vector<double,DIM>& Get(unsigned elementIndex, unsigned quadIndex);

    /**
     * Get the i-th stored quadrature point.
     *
     * \todo this method should be renamed rGet() as it returns a reference
     *
     * @param i
     */
    c_vector<double,DIM>& Get(unsigned i);

    /** Number of elements in the mesh that was given in the constructor */
    unsigned GetNumElements() const;

    /** Number of quad points per element in the rule that was given in the constructor */
    unsigned GetNumQuadPointsPerElement() const;

    /** Total size, ie total number of quad points, ie num_elem times num_quad_points_per_elem */
    unsigned Size() const;
};

#endif /*QUADRATUREPOINTSGROUP_HPP_*/
