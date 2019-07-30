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

#ifndef ABSTRACTCRYPTSTATISTICS_HPP_
#define ABSTRACTCRYPTSTATISTICS_HPP_

#include "MeshBasedCellPopulation.hpp"

/**
 * Abstract crypt statistics class.
 */
class AbstractCryptStatistics
{
protected:

    /** The crypt. */
    MeshBasedCellPopulation<2>& mrCrypt;

public:

    /**
     * Constructor.
     *
     * @param rCrypt The crypt
     */
    AbstractCryptStatistics(MeshBasedCellPopulation<2>& rCrypt);

    /**
     * Destructor.
     */
    virtual ~AbstractCryptStatistics();

    /**
     * To recreate the virtual labelling experiments performed by Meineke et al
     * (2001) in their off-lattice model of the intestinal crypt
     * (doi:10.1046/j.0960-7722.2001.00216.x).
     *
     * Cells which are in S phase are labelled using the CellLabel cell property.
     *
     * In Owen Sansom's experiments this is called twice; once at the
     * beginning and once at the end of an hour to simulate uptake of the
     * label over an hour, so some cells will already be labelled when this
     * is called the second time.
     *
     * (assumption that S phase lasts longer than one hour is pretty sound)
     */
    void LabelSPhaseCells();

    /**
     * Set all the cells in the crypt to have the mutation
     * state WildTypeCellMutationState.
     */
    void LabelAllCellsAsHealthy();

    /**
     * Get all cells within a cell width of the section defined as the line between points (xBottom,0)
     * and (xTop,yTop). If a patricular cell is labelled then the boolean true is returned.
     *
     * Periodicity can be taken into account (if xTop and xBottom are more than half a crypt
     * width apart then a more realistic section will be across the periodic boundary), using the
     * final parameter. This obviously requires the mesh to be cylindrical.
     *
     * @param rCryptSection  vector of cells in the section (from a call to GetCryptSection in the concrete class)
     *
     * @return  a standard vector of booleans which states whether a labelled cell is present at a corresponding position.
     */
    std::vector<bool> AreCryptSectionCellsLabelled(std::vector<CellPtr>& rCryptSection);
};

#endif /*ABSTRACTCRYPTSTATISTICS_HPP_*/
