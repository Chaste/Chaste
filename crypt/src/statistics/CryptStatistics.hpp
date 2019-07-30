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

#ifndef CRYPTSTATISTICS_HPP_
#define CRYPTSTATISTICS_HPP_

#include "AbstractCryptStatistics.hpp"

/**
 * Concrete statistics class for the cylindrical crypt model.
 */
class CryptStatistics : public AbstractCryptStatistics
{
private:

    /**
     * Method computing the perpendicular distance from the cell to the line from (xBottom, 0) to
     * (xTop, yTop), and returning if the distance is within the specified width to the section
     * (defaults to 0.5).
     *
     * @param xBottom the x coordinate of the bottom of the line
     * @param xTop the x coordinate of the top of the line
     * @param yTop the y coordinate of the top of the line
     * @param rCellPosition the location of the cell centre
     * @param widthOfSection the width of the line
     *
     * @return whether the cell is in the section
     */
    bool CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& rCellPosition, double widthOfSection=0.5);

    /**
     * Method computing the perpendicular distance from the cell to the line from (xBottom, 0) to
     * (xTop, yTop), taking into account periodicity, and returning if the distance is within the
     * specified width to the section (defaults to 1.0). Done by considering the two possible lines
     * and checking if cells are within range.
     *
     * @param xBottom the x coordinate of the bottom of the line
     * @param xTop the x coordinate of the top of the line
     * @param yTop the y coordinate of the top of the line
     * @param rCellPosition the location of the cell centre
     * @param widthOfSection the width of the line
     *
     * @return whether the cell is in the section
     */
    bool CellIsInSectionPeriodic(double xBottom, double xTop, double yTop, const c_vector<double,2>& rCellPosition, double widthOfSection=1.0);

public:

    /**
     * Constructor
     *
     * @param rCrypt The crypt
     */
    CryptStatistics(MeshBasedCellPopulation<2>& rCrypt);

    /**
     * Get all cells within a cell width of the section defined as the line between points (xBottom,0)
     * and (xTop,yTop)
     *
     * Periodicity can be taken into account (if xTop and xBottom are more than half a crypt
     * width apart then a more realistic section will be across the periodic boundary), using the
     * final parameter. This obviously requires the mesh to be cylindrical.
     * @param yTop the y coordinate of the top of the line
     * @param xBottom the x coordinate of the bottom of the line (defaults to a random number U[0,crypt_width])
     * @param xTop the x coordinate of the top of the line (defaults to a random number U[0,crypt_width])
     * @param periodic whether periodicity is accounted for (defaults to false)
     *
     * @return  an ordered list of CellPtrs from the bottom to the top of the crypt.
     *
     * Note that placing calls to functions with side-effects (eg. changing the random seed)
     * in the default arguments is DANGEROUS.  There is no guarantee that the compiler will
     * execute these in a sensible order.
     * It appears that Intel goes left-to-right and Gcc goes right-to-left.
     */
     std::vector<CellPtr> GetCryptSection(double yTop,
                                          double xBottom = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*crypt_width,
                                          double xTop = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*crypt_width,
                                          bool periodic = false);

    /**
     * Get all cells with a cell width of the line defined by the points (xBottom,0)
     * and (xTop,yTop), taking into account periodicity
     *
     * If xTop and xBottom are more than half a crypt width apart then a more realistic section
     * will be across the periodic boundary.
     *
     * Note that placing calls to functions with side-effects (eg. changing the random seed)
     * in the default arguments is DANGEROUS.  There is no guarantee that the compiler will
     * execute these in a sensible order. It appears that Intel goes left-to-right and Gcc goes
     * right-to-left.
     *
     * @param yTop the y coordinate of the top of the line
     * @param xBottom the x coordinate of the bottom of the line (defaults to a random number U[0,crypt_width])
     * @param xTop the x coordinate of the top of the line (defaults to a random number U[0,crypt_width])
     *
     * @return an ordered list of CellPtrs from the bottom to the top of the crypt.
     */
    std::vector<CellPtr> GetCryptSectionPeriodic(double yTop,
                                                 double xBottom = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*crypt_width,
                                                 double xTop = DBL_MAX); //RandomNumberGenerator::Instance()->ranf()*crypt_width,
};

#endif /*CRYPTSTATISTICS_HPP_*/
