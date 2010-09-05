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
#include "CryptProjectionStatistics.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * This global function is to allow the list of cells in to be compared in
 * terms of their y-value and std::list.sort() to be called
 */
bool CellsRadiusComparison(const std::pair<CellPtr, double> lhs, const std::pair<CellPtr, double> rhs)
{
    return lhs.second < rhs.second;
}

CryptProjectionStatistics::CryptProjectionStatistics(MeshBasedCellPopulation<2>& rCrypt)
    : AbstractCryptStatistics(rCrypt)
{
}

bool CryptProjectionStatistics::CellIsInSection(double angle, const c_vector<double,2>& rCellPosition, double widthOfSection)
{
    // Get corresponding 3D position of closest point on line
    c_vector<double,2> line_position;
    line_position[0] = norm_2(rCellPosition)*cos(angle);
    line_position[1] = norm_2(rCellPosition)*sin(angle);

    double distance_between_cell_and_line = norm_2(rCellPosition - line_position);

    return (distance_between_cell_and_line <= widthOfSection);
}

std::vector<CellPtr> CryptProjectionStatistics::GetCryptSection(double angle)
{
    if (angle == DBL_MAX)
    {
        angle = M_PI - 2*M_PI*RandomNumberGenerator::Instance()->ranf();
    }

    assert(angle>=-M_PI && angle<=M_PI);

    std::list<std::pair<CellPtr, double> > cells_list; // the second entry is the radius (needed for sorting)


    // Loop over cells and add to the store if they are within a cell's radius of the
    // specified line
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {

        if ( CellIsInSection(angle, mrCrypt.GetLocationOfCellCentre(*cell_iter)) )
        {
            // Set up a pair, equal to (cell,r) and insert
            std::pair<CellPtr, double> pair(*cell_iter, norm_2(mrCrypt.GetLocationOfCellCentre(*cell_iter)));
            cells_list.push_back(pair);
        }
    }

    // Sort the list
    cells_list.sort(CellsRadiusComparison);

    // Copy to a vector
    std::vector<CellPtr> ordered_cells;
    for (std::list<std::pair<CellPtr, double> >::iterator iter = cells_list.begin();
         iter != cells_list.end();
         ++iter)
    {
        ordered_cells.push_back(iter->first);
    }

    return ordered_cells;
}
