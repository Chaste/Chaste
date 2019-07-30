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


#ifndef _TESTMONODOMAINNOSTIMULUS_HPP_
#define _TESTMONODOMAINNOSTIMULUS_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>

#include "PropagationPropertiesCalculator.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"

class ZeroStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
public:
    ZeroStimulusCellFactory() : AbstractCardiacCellFactory<1>()
    {}

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);

    }
};


/* TestMonodomainNoStimulus - based on TestMonodomainConductionVelocity
 *
 * No initial stimulus applied.
 * Check that the voltage of all cells is constant thoroughout the mesh
 * at any point in time and never lower than the resting potential
 * of the LR cell = -85.0 mV
 *
 * Best run with optimisation on.
 */
class TestMonodomainNoStimulus : public CxxTest::TestSuite
{
public:

    void TestZeroStimulus()
    {
        HeartConfig::Instance()->SetSimulationDuration(30); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_20_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoNoStim");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainNoStimLR91_1d");

        ZeroStimulusCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double constant_voltage = voltage_replicated[0];
        TS_ASSERT_LESS_THAN(-85.0, constant_voltage);

        for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index] , constant_voltage, 1E-5);
        }
    }
};
#endif //_TESTMONODOMAINNOSTIMULUS_HPP_
