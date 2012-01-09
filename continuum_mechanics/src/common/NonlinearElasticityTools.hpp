/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef NONLINEARELASTICITYTOOLS_HPP_
#define NONLINEARELASTICITYTOOLS_HPP_

#include "TetrahedralMesh.hpp"

/**
 *  A class of helper methods for problems which use NonlinearElasticitySolver.
 */
template<unsigned DIM>
class NonlinearElasticityTools
{
public:

    /**
     * Collect all the nodes which satisfy x[k] = c, for given k and c, in order
     * to be set as fixed (or displacement boundary condition) nodes. Note that
     * this method does not check if the nodes on the required surface are actually
     * boundary nodes. It does however throw an exception if no nodes on the given
     * surface are found.
     *
     * @param rMesh  the mesh
     * @param component  the component k
     * @param value  the value c
     */
    static std::vector<unsigned> GetNodesByComponentValue(TetrahedralMesh<DIM,DIM>& rMesh,
                                                          unsigned component,
                                                          double value);
};

#endif /*NONLINEARELASTICITYTOOLS_HPP_*/
