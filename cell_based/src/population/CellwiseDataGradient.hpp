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

#ifndef CELLWISEDATAGRADIENT_HPP_
#define CELLWISEDATAGRADIENT_HPP_

#include "CellwiseData.hpp"

/**
 *  A class for calculating the gradients of the CellwiseData.
 */
template<unsigned DIM>
class CellwiseDataGradient
{
private:

    /**
     * The final gradients at the nodes
     */
    std::vector<c_vector<double, DIM> > mGradients;

public:

    /**
     * Compute the gradients at the nodes.
     *
     * This is done by averaging the gradients at all the containing (non-ghost)
     * elements for that node. Note that the gradients are piecewise constant-
     * constant in each element
     */
    void SetupGradients();

    /**
     * Get the gradient at a given node. Not set up for ghost nodes.
     *
     * @param nodeIndex
     *
     * @return the gradient at the node.
     */
    c_vector<double, DIM>& rGetGradient(unsigned nodeIndex);

};

#endif /*CELLWISEDATAGRADIENT_HPP_*/
