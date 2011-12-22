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

#include "NonlinearElasticityTools.hpp"

template<unsigned DIM>
std::vector<unsigned> NonlinearElasticityTools<DIM>::GetNodesByComponentValue(TetrahedralMesh<DIM,DIM>& rMesh,
                                                          unsigned component,
                                                          double value)
{
    std::vector<unsigned> fixed_nodes;
    double tol = 1e-8;
    for (unsigned i=0; i<rMesh.GetNumNodes(); i++)
    {
        if ( fabs(rMesh.GetNode(i)->rGetLocation()[component] - value)<1e-8)
        {
            fixed_nodes.push_back(i);
        }
    }

    if (fixed_nodes.size() == 0)
    {
        EXCEPTION("Could not find any nodes on requested surface (note: tolerance = "<<tol<<")");
    }

    return fixed_nodes;
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

//template class NonlinearElasticityTools<1>;
template class NonlinearElasticityTools<2>;
template class NonlinearElasticityTools<3>;
