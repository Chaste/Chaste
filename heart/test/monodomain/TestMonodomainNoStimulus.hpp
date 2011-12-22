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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
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
