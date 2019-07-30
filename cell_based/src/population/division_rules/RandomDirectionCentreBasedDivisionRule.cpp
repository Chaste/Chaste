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

#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > RandomDirectionCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Get separation parameter
    double separation = rCellPopulation.GetMeinekeDivisionSeparation();

    // Make a random direction vector of the required length
    c_vector<double, SPACE_DIM> random_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */
    switch (SPACE_DIM)
    {
        case 1:
        {
            double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

            random_vector(0) = 0.5*separation*random_direction;
            break;
        }
        case 2:
        {
            double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

            random_vector(0) = 0.5*separation*cos(random_angle);
            random_vector(1) = 0.5*separation*sin(random_angle);
            break;
        }
        case 3:
        {
            /*
             * Note that to pick a random point on the surface of a sphere, it is incorrect
             * to select spherical coordinates from uniform distributions on [0, 2*pi) and
             * [0, pi) respectively, since points picked in this way will be 'bunched' near
             * the poles. See #2230.
             */
            double u = RandomNumberGenerator::Instance()->ranf();
            double v = RandomNumberGenerator::Instance()->ranf();

            double random_azimuth_angle = 2*M_PI*u;
            double random_zenith_angle = std::acos(2*v - 1);

            random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
            random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
            random_vector(2) = 0.5*separation*cos(random_zenith_angle);
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
    }

    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - random_vector;
    c_vector<double, SPACE_DIM> daughter_position = parent_position + random_vector;

    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);

    return positions;
}

// Explicit instantiation
template class RandomDirectionCentreBasedDivisionRule<1,1>;
template class RandomDirectionCentreBasedDivisionRule<1,2>;
template class RandomDirectionCentreBasedDivisionRule<2,2>;
template class RandomDirectionCentreBasedDivisionRule<1,3>;
template class RandomDirectionCentreBasedDivisionRule<2,3>;
template class RandomDirectionCentreBasedDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RandomDirectionCentreBasedDivisionRule)
