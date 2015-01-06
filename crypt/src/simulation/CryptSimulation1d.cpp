/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "CryptSimulation1d.hpp"
#include "WntConcentration.hpp"
#include "SmartPointers.hpp"

CryptSimulation1d::CryptSimulation1d(AbstractCellPopulation<1>& rCellPopulation,
                  bool deleteCellPopulationInDestructor,
                  bool initialiseCells)
    : OffLatticeSimulation<1>(rCellPopulation,
                          deleteCellPopulationInDestructor,
                          initialiseCells)
{
    mpStaticCastCellPopulation = static_cast<MeshBasedCellPopulation<1>*>(&mrCellPopulation);

    if (!mDeleteCellPopulationInDestructor)
    {
        // Pass a CryptSimulationBoundaryCondition object into mBoundaryConditions
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<1>, p_bc, (&rCellPopulation));
        AddCellPopulationBoundaryCondition(p_bc);
    }
}

CryptSimulation1d::~CryptSimulation1d()
{
}

c_vector<double, 1> CryptSimulation1d::CalculateCellDivisionVector(CellPtr pParentCell)
{
    // Location of parent and daughter cells
    c_vector<double, 1> parent_coords = mpStaticCastCellPopulation->GetLocationOfCellCentre(pParentCell);
    c_vector<double, 1> daughter_coords;

    // Get separation parameter
    double separation = mpStaticCastCellPopulation->GetMeinekeDivisionSeparation();

    // Make a random direction vector of the required length
    c_vector<double, 1> random_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */

    double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);
    random_vector(0) = 0.5*separation*random_direction;
    c_vector<double, 1> proposed_new_parent_coords = parent_coords - random_vector;
    c_vector<double, 1> proposed_new_daughter_coords = parent_coords + random_vector;

    if (   (proposed_new_parent_coords(0) >= 0.0)
        && (proposed_new_daughter_coords(0) >= 0.0))
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

    assert(daughter_coords(0) >= 0.0); // to make sure dividing cells stay in the cell population
    assert(parent_coords(0) >= 0.0);   // to make sure dividing cells stay in the cell population

    // Set the parent to use this location
    ChastePoint<1> parent_coords_point(parent_coords);

    unsigned node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);
    mrCellPopulation.SetNode(node_index, parent_coords_point);

    return daughter_coords;
}

void CryptSimulation1d::OutputSimulationParameters(out_stream& rParamsFile)
{
    // No parameters to output

    // Call method on direct parent class
    OffLatticeSimulation<1>::OutputSimulationParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation1d)
