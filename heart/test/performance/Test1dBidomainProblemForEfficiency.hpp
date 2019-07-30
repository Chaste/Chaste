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

            std::vector<double> ode_vars = bidomain_problem.GetBidomainTissue()->GetCardiacCell(index.Global)->GetStdVecStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=0) && (j!=7))
                {
                    TS_ASSERT_LESS_THAN_EQUALS( ode_vars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS(-ode_vars[j], 0.0);
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
