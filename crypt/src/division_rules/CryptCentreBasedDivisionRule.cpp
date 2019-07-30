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

#include "CryptCentreBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CryptCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Get separation parameter
    double separation = rCellPopulation.GetMeinekeDivisionSeparation();

    // Make a random direction vector of the required length
    c_vector<double, SPACE_DIM> random_vector;

    c_vector<double, SPACE_DIM> parent_coords = rCellPopulation.GetLocationOfCellCentre(pParentCell);
    c_vector<double, SPACE_DIM> daughter_coords;

    switch (SPACE_DIM)
    {
        case 1:
        {
            double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);
            random_vector(0) = 0.5*separation*random_direction;
            c_vector<double, SPACE_DIM> proposed_new_parent_coords;
            //proposed_new_parent_coords = parent_coords - random_vector;
            proposed_new_parent_coords(0) = parent_coords(0) - random_vector(0);
            c_vector<double, SPACE_DIM> proposed_new_daughter_coords;
            proposed_new_daughter_coords = parent_coords + random_vector;

            if ((proposed_new_parent_coords(0) >= 0.0) && (proposed_new_daughter_coords(0) >= 0.0))
            {
                // We are not too close to the bottom of the cell population, so move parent
                parent_coords = proposed_new_parent_coords;
                daughter_coords = proposed_new_daughter_coords;
            }
            else
            {
                proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
                while (proposed_new_daughter_coords(0) < 0.0)
                {
                    double fresh_random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);
                    random_vector(0) = 0.5*separation*fresh_random_direction;
                    proposed_new_daughter_coords = parent_coords + random_vector;
                }
                daughter_coords = proposed_new_daughter_coords;
            }

            // To make sure that dividing cells stay in the cell population
            assert(daughter_coords(0) >= 0.0);
            assert(parent_coords(0) >= 0.0);

            break;
        }
        case 2:
        {
            /*
             * Pick a random direction and move the parent cell backwards by 0.5*separation
             * in that direction and return the position of the daughter cell 0.5*separation
             * forwards in that direction.
             */
            double random_angle = RandomNumberGenerator::Instance()->ranf();
            random_angle *= 2.0*M_PI;

            random_vector(0) = 0.5*separation*cos(random_angle);
            random_vector(1) = 0.5*separation*sin(random_angle);

            c_vector<double, 2> proposed_new_parent_coords;
            proposed_new_parent_coords = parent_coords - random_vector;
            c_vector<double, 2> proposed_new_daughter_coords;
            proposed_new_daughter_coords = parent_coords + random_vector;

            if ((proposed_new_parent_coords(1) >= 0.0) && (proposed_new_daughter_coords(1) >= 0.0))
            {
                // We are not too close to the bottom of the cell population, so move parent
                parent_coords = proposed_new_parent_coords;
                daughter_coords = proposed_new_daughter_coords;
            }
            else
            {
                proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
                while (proposed_new_daughter_coords(1) < 0.0)
                {
                    random_angle = RandomNumberGenerator::Instance()->ranf();
                    random_angle *= 2.0*M_PI;

                    random_vector(0) = separation*cos(random_angle);
                    random_vector(1) = separation*sin(random_angle);
                    proposed_new_daughter_coords = parent_coords + random_vector;
                }
                daughter_coords = proposed_new_daughter_coords;
            }

            // To make sure that dividing cells stay in the cell population
            assert(daughter_coords(1) >= 0.0);
            assert(parent_coords(1) >= 0.0);
            break;
        }
        case 3:
        {
            EXCEPTION("CryptCentreBasedDivisionRule is not implemented for SPACE_DIM == 3");
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
    }

    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_coords, daughter_coords);
    return positions;
}

// Explicit instantiation
template class CryptCentreBasedDivisionRule<1,1>;
template class CryptCentreBasedDivisionRule<1,2>;
template class CryptCentreBasedDivisionRule<2,2>;
template class CryptCentreBasedDivisionRule<1,3>;
template class CryptCentreBasedDivisionRule<2,3>;
template class CryptCentreBasedDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CryptCentreBasedDivisionRule)


