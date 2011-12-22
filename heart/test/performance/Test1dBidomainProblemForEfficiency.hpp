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


#ifndef TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_
#define TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_


#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudy1991.hpp"


class Test1dBidomainProblemForEfficiency : public CxxTest::TestSuite
{
public:
    void TestBidomainDg01WithNoOutput()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.00005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.00005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_1000_elements");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        bidomain_problem.PrintOutput(false);

        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            TS_FAIL(e.GetMessage());
        }

        DistributedVector striped_voltage = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe voltage(striped_voltage, 0);

        for (DistributedVector::Iterator index = striped_voltage.Begin();
             index != striped_voltage.End();
             ++index)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV

            TS_ASSERT_LESS_THAN_EQUALS( voltage[index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-voltage[index] + (Ek-30), 0);

            std::vector<double>& r_ode_vars = bidomain_problem.GetBidomainTissue()->GetCardiacCell(index.Global)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=0) && (j!=7))
                {
                    TS_ASSERT_LESS_THAN_EQUALS( r_ode_vars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS(-r_ode_vars[j], 0.0);
                }
            }


            // final voltages for six nodes at the beginning of the mesh with a stride of 10
            double test_values[6]={11.5550, -78.3303, -83.7585, -83.8568,  -83.8570, -83.8568};

            for (unsigned i=0; i<=5; i++)
            {
                unsigned node=10*i; //Step through every 10th node
                if (index.Global == node)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], test_values[i], 1e-1);
                }
            }
        }
    }
};

#endif /*TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_*/
