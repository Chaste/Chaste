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
#include "CryptStatistics.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * This global function is to allow the list of cells in to be compared in
 * terms of their y-value and std::list.sort() to be called
 */
bool CellsHeightComparison(const std::pair<CellPtr, double> lhs, const std::pair<CellPtr, double> rhs)
{
    return lhs.second < rhs.second;
}

CryptStatistics::CryptStatistics(MeshBasedCellPopulation<2>& rCrypt)
    : AbstractCryptStatistics(rCrypt)
{
}

std::vector<CellPtr> CryptStatistics::GetCryptSection(double yTop, double xBottom, double xTop, bool periodic)
{
    double crypt_width = mrCrypt.rGetMesh().GetWidth(0);

    // Fill in the default values - in a sequential manner
    if (xBottom == DBL_MAX)
    {
        xBottom = RandomNumberGenerator::Instance()->ranf()*crypt_width;
    }

    if (xTop == DBL_MAX)
    {
        xTop = RandomNumberGenerator::Instance()->ranf()*crypt_width;
    }

    assert(yTop>0.0);
    std::list<std::pair<CellPtr, double> > cells_list; // the second entry is the y value (needed for sorting)

    if (fabs(xTop-xBottom)<0.5*crypt_width)
    {
        // The periodic version isn't needed, ignore even if periodic was set to true
        periodic = false;
    }

    // Loop over cells and add to the store if they are within a cell's radius of the specified line
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        if (periodic)
        {
            if (CellIsInSectionPeriodic(xBottom, xTop, yTop, mrCrypt.GetLocationOfCellCentre(*cell_iter)))
            {
                // Set up a pair, equal to (cell,y_val) and insert
                std::pair<CellPtr, double> pair(*cell_iter, mrCrypt.GetLocationOfCellCentre(*cell_iter)[1]);
                cells_list.push_back(pair);
            }
        }
        else
        {
            if (CellIsInSection(xBottom, xTop, yTop, mrCrypt.GetLocationOfCellCentre(*cell_iter)))
            {
                // Set up a pair, equal to (cell,y_val) and insert
                std::pair<CellPtr, double> pair(*cell_iter, mrCrypt.GetLocationOfCellCentre(*cell_iter)[1]);
                cells_list.push_back(pair);
            }
        }
    }

    // Sort the list
    cells_list.sort(CellsHeightComparison);

    // Copy to a vector
    std::vector<CellPtr> ordered_cells;
    for (std::list<std::pair<CellPtr, double> >::iterator iter = cells_list.begin();
         iter!=cells_list.end();
         ++iter)
    {
        ordered_cells.push_back(iter->first);
    }

    return ordered_cells;
}

std::vector<CellPtr> CryptStatistics::GetCryptSectionPeriodic(double yTop, double xBottom, double xTop)
{
   return GetCryptSection(yTop, xBottom, xTop, true);
}
bool CryptStatistics::CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& rCellPosition, double widthOfSection)
{
    c_vector<double,2> intercept;

    if (xBottom == xTop)
    {
        intercept[0] = xTop;
        intercept[1] = rCellPosition[1];
    }
    else
    {
        double m = (yTop)/(xTop-xBottom); // gradient of line

        intercept[0] = (m*m*xBottom + rCellPosition[0] + m*rCellPosition[1])/(1+m*m);
        intercept[1] = m*(intercept[0] - xBottom);
    }

    c_vector<double,2> vec_from_A_to_B = mrCrypt.rGetMesh().GetVectorFromAtoB(intercept, rCellPosition);
    double dist = norm_2(vec_from_A_to_B);

    return (dist <= widthOfSection);
}

bool CryptStatistics::CellIsInSectionPeriodic(double xBottom, double xTop, double yTop, const c_vector<double,2>& rCellPosition, double widthOfSection)
{
    bool is_in_section = false;

    c_vector<double,2> intercept;
    double crypt_width = mrCrypt.rGetMesh().GetWidth(0u);

    double m; // gradient of line
    double offset;

    if (xBottom < xTop)
    {
        offset = -crypt_width;
    }
    else
    {
        offset = crypt_width;
    }

    m = (yTop)/(xTop-xBottom+offset); // gradient of line

    // 1st line
    intercept[0] = (m*m*xBottom + rCellPosition[0] + m*rCellPosition[1])/(1+m*m);
    intercept[1] = m*(intercept[0] - xBottom);

    c_vector<double,2> vec_from_A_to_B = mrCrypt.rGetMesh().GetVectorFromAtoB(intercept, rCellPosition);
    double dist = norm_2(vec_from_A_to_B);

    if (dist < widthOfSection)
    {
        is_in_section = true;
    }

    // 2nd line
    intercept[0] = (m*m*(xBottom-offset) + rCellPosition[0] + m*rCellPosition[1])/(1+m*m);
    intercept[1] = m*(intercept[0] - (xBottom-offset));

    vec_from_A_to_B = mrCrypt.rGetMesh().GetVectorFromAtoB(intercept, rCellPosition);
    dist = norm_2(vec_from_A_to_B);

    if (dist < widthOfSection)
    {
        is_in_section = true;
    }

    return is_in_section;
}
