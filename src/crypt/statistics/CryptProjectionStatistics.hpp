/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef CRYPTPROJECTIONSTATISTICS_HPP_
#define CRYPTPROJECTIONSTATISTICS_HPP_

#include "AbstractCryptStatistics.hpp"

/**
 * Concrete statistics class for the crypt projection model.
 */
class CryptProjectionStatistics : public AbstractCryptStatistics
{
protected:

    /**
     * CellIsInSection method.
     *
     * @param angle  The angle between the crypt section and the x axis in the projection
     * @param rCellPosition  The vector of a cell's position
     * @param widthOfSection The width of the section
     */
    bool CellIsInSection(double angle, const c_vector<double,2>& rCellPosition, double widthOfSection=0.6);

public:

    /**
     *  Constructor
     *
     *  @param rCrypt  The crypt
     */
    CryptProjectionStatistics(MeshBasedTissue<2>& rCrypt)
        : AbstractCryptStatistics(rCrypt)
    {}

    /**
     * GetCryptSection method. Takes in an angle from the
     * interval (-pi, pi].
     *
     * @param angle  The angle between the crypt section and the x axis in the projection
     */
    std::vector<TissueCell*> GetCryptSection(double angle = DBL_MAX);

};

#endif /*CRYPTPROJECTIONSTATISTICS_HPP_*/
