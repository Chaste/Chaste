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

#ifndef CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_
#define CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "AbstractCellCycleModel.hpp"

/**
 * A helper method that is called in cell-cycle model tests
 * to check that a cell is progressing through the cell cycle
 * correctly. This is a method which should only be called by
 * cxx-test classes (not source code) as it includes TS_ASSERT
 * calls.
 *
 * @param pModel Pointer to the cell-cycle model
 * @param g1Duration Correct duration of the G1 phase, to test against
 * @param g2Duration Correct duration of the G1 phase, to test against
 */
void CheckReadyToDivideAndPhaseIsUpdated(AbstractCellCycleModel* pModel,
                                         double g1Duration,
                                         double g2Duration=DBL_MAX)
{
    if (g2Duration==DBL_MAX)
    {
        g2Duration = pModel->GetG2Duration();
    }


    double age = pModel->GetAge();

    const double G1TOL = 1e-5; // how accurate the expected G1 duration is

    // If the G1 duration is incorrect, print out the mismatch
    if ((pModel->GetCellProliferativeType() != DIFFERENTIATED) &&
        (age >= pModel->GetMDuration()) &&
        (pModel->GetG1Duration() != DOUBLE_UNSET) &&
        (fabs(pModel->GetG1Duration() - g1Duration) > G1TOL))
    {
        std::cout << "G1 duration mismatch: actual = " << pModel->GetG1Duration()
                  << ", expected = " << g1Duration
                  << std::endl;
    }

    if (pModel->GetCellProliferativeType()==DIFFERENTIATED)
    {
        // If the cell is differentiated, then it must be in G0 phase and must never divide
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(), G_ZERO_PHASE);
    }
    else if (age < pModel->GetMDuration())
    {
        // If the cell in M phase, then it must not be ready to divide
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(), M_PHASE);
    }
    else if (age < pModel->GetMDuration() + g1Duration - G1TOL)
    {
        // The next cell cycle phase after M is G1; cells in G1 phase
        // must still not be ready to divide
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(), G_ONE_PHASE);

        // If the cell is not in G1 phase when it should be, print out the mismatch
        if (pModel->GetCurrentCellCyclePhase() != G_ONE_PHASE)
        {
            std::cout << "Expected G1: " << g1Duration
                      << "; actual: " << pModel->GetG1Duration()
                      << "; age = " << age
                      << "; G1-S transition = " << pModel->GetMDuration() + g1Duration
                      << std::endl;
        }
    }
    else if (age < pModel->GetMDuration() + g1Duration + pModel->GetSDuration() - G1TOL)
    {
        // The next cell cycle phase after G1 is S; cells in S phase
        // must still not be ready to divide
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(), S_PHASE);
    }
    else if (age < pModel->GetMDuration() + g1Duration + pModel->GetSDuration() + g2Duration  - G1TOL)
    {
        // The next cell cycle phase after S is G2; cells in G2 phase
        // must still not be ready to divide
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(), G_TWO_PHASE);
    }
    else
    {
        // Cells must be ready to divide as soon as they leave G2 phase
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), true);
    }
}

#endif /*CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_*/
