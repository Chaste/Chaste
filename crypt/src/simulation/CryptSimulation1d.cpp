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
            double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);
            random_vector(0) = 0.5*separation*random_direction;
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
