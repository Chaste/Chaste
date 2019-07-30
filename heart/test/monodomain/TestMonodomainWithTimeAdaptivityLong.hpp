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

#ifndef TESTMONODOMAINWITHTIMEADAPTIVITYLONG_HPP_
#define TESTMONODOMAINWITHTIMEADAPTIVITYLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "ReplicatableVector.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscVecTools.hpp"


// Toy controller which just goes alters the timestep from 0.01ms to 1ms after a given
// threshold time.
class FixedTimeAdaptivityController : public AbstractTimeAdaptivityController
{
private:
    double mThresholdTime;

    double ComputeTimeStep(double currentTime, Vec currentSolution)
    {
        if (currentTime < mThresholdTime)
        {
            return 0.01; // ms
        }
        else
        {
            return 1;
        }
    }


public:
    FixedTimeAdaptivityController(double thresholdTime)
      : AbstractTimeAdaptivityController(0.01, 1.0),
        mThresholdTime(thresholdTime)
    {
        assert(thresholdTime > 0);
    }
};


class TestMonodomainWithTimeAdaptivity : public CxxTest::TestSuite
{
public:
    void Test1dApd()
    {
        HeartConfig::Instance()->SetPrintingTimeStep(1.0);
        HeartConfig::Instance()->SetSimulationDuration(400); //ms

        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.01, 1.0); // h=0.01cm, width=1cm

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory(-600.0*1000);

        //////////////////////////////////////////////////////////////////////////
        // run original simulation - no adaptivity, dt=0.01 all the way through
        //////////////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("MonoWithTimeAdaptivity1dLong/OrigNoAdapt");
        MonodomainProblem<1> problem(&cell_factory);
        problem.SetMesh(&mesh);

        problem.Initialise();
        problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        //////////////////////////////////////////////////////////////////////////
        // run adaptive simulation - dt=0.01 for first 2ms, then dt=1
        //////////////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("MonoWithTimeAdaptivity1dLong/SimpleAdapt");
        MonodomainProblem<1> adaptive_problem(&cell_factory);
        adaptive_problem.SetMesh(&mesh);

        FixedTimeAdaptivityController controller(25);
        adaptive_problem.SetUseTimeAdaptivityController(true, &controller);
        adaptive_problem.Initialise();
        adaptive_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        Hdf5DataReader reader_no_adapt("MonoWithTimeAdaptivity1dLong/OrigNoAdapt","SimulationResults");
        Hdf5DataReader reader_adapt("MonoWithTimeAdaptivity1dLong/SimpleAdapt","SimulationResults");

        unsigned num_timesteps = reader_no_adapt.GetUnlimitedDimensionValues().size();
        assert(num_timesteps == reader_adapt.GetUnlimitedDimensionValues().size());

        DistributedVectorFactory factory(mesh.GetNumNodes());
        Vec voltage_no_adapt = factory.CreateVec();
        Vec voltage_adapt = factory.CreateVec();

        Vec difference;
        VecDuplicate(voltage_adapt, &difference);

        for (unsigned timestep=0; timestep<num_timesteps; timestep++)
        {
            reader_no_adapt.GetVariableOverNodes(voltage_no_adapt, "V", timestep);
            reader_adapt.GetVariableOverNodes(voltage_adapt, "V", timestep);

            PetscVecTools::WAXPY(difference, -1.0, voltage_adapt, voltage_no_adapt);
            double l_inf_norm;
            VecNorm(difference, NORM_INFINITY, &l_inf_norm);

            //std::cout << l_inf_norm << "\n";
            if (timestep < 25)
            {
                TS_ASSERT_DELTA(l_inf_norm, 0.0, 1e-10); // first 25 ms, there should be no difference
            }
            else
            {
                TS_ASSERT_DELTA(l_inf_norm, 0.0, 2.25); // the difference is at most ~2mv, which occurs during the downstroke
            }
        }

        PetscTools::Destroy(voltage_no_adapt);
        PetscTools::Destroy(voltage_adapt);
    }
};

#endif /*TESTMONODOMAINWITHTIMEADAPTIVITYLONG_HPP_*/

//
//
//Entering Test1dApd
//         InMesh             Init           AssSys              Ode            Comms           AssRhs           NeuBCs           DirBCs              Ksp           Output         PostProc            User1            Total
//   0.000 (  0%)     0.000 (  0%)     0.003 (  0%)    13.342 ( 61%)     0.904 (  4%)     0.830 (  4%)     0.048 (  0%)     0.000 (  0%)     4.915 ( 23%)     0.171 (  1%)     0.124 (  1%)     0.000 (  0%)    21.720 (100%)  (seconds)
//         InMesh             Init           AssSys              Ode            Comms           AssRhs           NeuBCs           DirBCs              Ksp           Output         PostProc            User1            Total
//   0.000 (  0%)     0.001 (  0%)     0.006 (  0%)     7.018 ( 84%)     0.077 (  1%)     0.066 (  1%)     0.004 (  0%)     0.000 (  0%)     0.388 (  5%)     0.506 (  6%)     0.121 (  1%)     0.000 (  0%)     8.323 (100%)  (seconds)
