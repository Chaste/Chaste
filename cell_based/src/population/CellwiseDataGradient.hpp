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

#ifndef CELLWISEDATAGRADIENT_HPP_
#define CELLWISEDATAGRADIENT_HPP_

#include <MeshBasedCellPopulation.hpp>

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
     *
     * @param rCellPopulation population on which to calculate gradients - must be instantiation of MeshBasedCellPopulation
     * @param rItemName is the name of the data from which to form the gradient (e.g. "oxygen").
     */
    void SetupGradients(AbstractCellPopulation<DIM>& rCellPopulation, const std::string& rItemName);

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
