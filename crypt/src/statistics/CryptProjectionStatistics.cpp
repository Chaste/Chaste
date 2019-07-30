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
        if (CellIsInSection(angle, mrCrypt.GetLocationOfCellCentre(*cell_iter)))
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
