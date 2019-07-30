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

#ifndef CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_
#define CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "AbstractSimpleCellCycleModel.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

/**
 * A helper method that is called in cell-cycle model tests
 * to check that a cell is progressing through the cell cycle
 * correctly. This is a method which should only be called by
 * cxx-test classes (not source code) as it includes TS_ASSERT
 * calls.
 *
 * @param pModel Pointer to the cell-cycle model (note must be a phase based model)
 * @param g1Duration Correct duration of the G1 phase, to test against
 * @param g2Duration Correct duration of the G1 phase, to test against
 */
void CheckReadyToDivideAndPhaseIsUpdated(AbstractPhaseBasedCellCycleModel* pModel,
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
    if ((pModel->GetCell()->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) &&
        (age >= pModel->GetMDuration()) &&
        (pModel->GetG1Duration() != DOUBLE_UNSET) &&
        (fabs(pModel->GetG1Duration() - g1Duration) > G1TOL))
    {
        std::cout << "G1 duration mismatch: actual = " << pModel->GetG1Duration()
                  << ", expected = " << g1Duration
                  << std::endl;
    }

    if (pModel->GetCell()->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
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

/**
 * A helper method that is called in cell-cycle model tests
 * to check that a cell is progressing through the cell cycle
 * correctly. This is a method which should only be called by
 * cxx-test classes (not source code) as it includes TS_ASSERT
 * calls.
 *
 * @param pModel Pointer to the cell-cycle model (note must be a simple non phase based model)
 * @param cellCycleDuration Correct duration of the cell cycle phase, to test against
 */
void CheckReadyToDivideIsUpdated(AbstractSimpleCellCycleModel* pModel,
                                 double cellCycleDuration)
{
    double age = pModel->GetAge();

    const double CCDTOL = 1e-5; // how accurate the expected CCD duration is

    // If the CCD duration is incorrect, print out the mismatch
    if ((pModel->GetCellCycleDuration() != DOUBLE_UNSET) &&
        (fabs(pModel->GetCellCycleDuration() - cellCycleDuration) > CCDTOL))
    {
        std::cout << "CCD duration mismatch: actual = " << pModel->GetCellCycleDuration()
                  << ", expected = " << cellCycleDuration
                  << std::endl;
    }

    // Check ReadyToDivide is updated correctly.
    if (pModel->GetCell()->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        // If the cell is differentiated, then it must be in G0 phase and must never divide
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
    }
    else if (age < cellCycleDuration)
    {
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), false);
    }
    else
    {
        TS_ASSERT_EQUALS(pModel->ReadyToDivide(), true);
    }
}


#endif /*CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_*/
