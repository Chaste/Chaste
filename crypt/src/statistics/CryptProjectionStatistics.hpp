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

#ifndef CRYPTPROJECTIONSTATISTICS_HPP_
#define CRYPTPROJECTIONSTATISTICS_HPP_

#include "AbstractCryptStatistics.hpp"

/**
 * Concrete statistics class for the crypt projection model.
 */
class CryptProjectionStatistics : public AbstractCryptStatistics
{
private:

    /**
     * CellIsInSection method.
     *
     * @param angle  The angle between the crypt section and the x axis in the projection
     * @param rCellPosition  The vector of a cell's position
     * @param widthOfSection The width of the section
     * @return whether the cell is in the section
     */
    bool CellIsInSection(double angle, const c_vector<double,2>& rCellPosition, double widthOfSection=0.6);

public:

    /**
     * Constructor.
     *
     * @param rCrypt The crypt
     */
    CryptProjectionStatistics(MeshBasedCellPopulation<2>& rCrypt);

    /**
     * GetCryptSection method. Takes in an angle from the
     * interval (-pi, pi].
     *
     * @param angle  The angle between the crypt section and the x axis in the projection
     * @return all cells in the section
     */
    std::vector<CellPtr> GetCryptSection(double angle = DBL_MAX);
};

#endif /*CRYPTPROJECTIONSTATISTICS_HPP_*/
