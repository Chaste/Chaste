/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTPARABOLICPDEANDBOUNDARYCONDITIONS_HPP_
#define TESTPARABOLICPDEANDBOUNDARYCONDITIONS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "Timer.hpp"
#include "ParabolicPdeAndBoundaryConditions.hpp"
#include "ConstBoundaryCondition.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "ReplicatableVector.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

template <int SPACE_DIM>
class HeatEquation : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& , double, Element<SPACE_DIM,SPACE_DIM>*)
    {
        return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return identity_matrix<double>(SPACE_DIM);
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& )
    {
        return 1;
    }

};

/**
 * For use in TestParabolicPdeAndBoundaryConditions::TestWithBoundaryConditionVaryingInSpace.
 */
double bc_func1(const ChastePoint<2>& p)
{
    return p[1]*p[1];
}

/**
 * For use in TestParabolicPdeAndBoundaryConditions::TestWithBoundaryConditionVaryingInTime.
 */
double bc_func2(const ChastePoint<2>& p)
{
    SimulationTime* p_time = SimulationTime::Instance();
    double value = 1.0 + 0.5*p_time->GetTime();
    return value;
}

class TestParabolicPdeAndBoundaryConditions : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestMethods() throw(Exception)
    {
        // Create a ParabolicPdeAndBoundaryConditions object
        HeatEquation<2> pde;
        ConstBoundaryCondition<2> bc(15.0);
        bool is_neumann_bc = false;

        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, is_neumann_bc);

        // Make sure it has no name
        TS_ASSERT_EQUALS(pde_and_bc.rGetDependentVariableName(), "");
        pde_and_bc.SetDependentVariableName("something");
        TS_ASSERT_EQUALS(pde_and_bc.rGetDependentVariableName(), "something");

        // Test Get methods
        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 15.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), false);

        bool solution_exists = pde_and_bc.GetSolution();
        TS_ASSERT_EQUALS(solution_exists, false);

        AbstractLinearParabolicPde<2,2>* p_pde = pde_and_bc.GetPde();
        TS_ASSERT_EQUALS(p_pde, &pde);

        // Set mCurrentSolution
        std::vector<double> data(10);
        for (unsigned i=0; i<10; i++)
        {
            data[i] = i + 0.45;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        // Test mCurrentSolution has been correctly set
        Vec solution =  pde_and_bc.GetSolution();
        ReplicatableVector solution_repl(solution);

        TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-12);
        }

        PetscInt size_of_solution = 0;
        VecGetSize(pde_and_bc.GetSolution(), &size_of_solution);
        TS_ASSERT_EQUALS(size_of_solution, 10);
    }

    void TestMethodsNeumann() throw(Exception)
    {
        // Create a ParabolicPdeAndBoundaryConditions object
        HeatEquation<2> pde;
        ConstBoundaryCondition<2> bc(0.0);

        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc); // third argument defaults to Neumann

        // Test Get methods
        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 0.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), true);

        bool solution_exists = pde_and_bc.GetSolution();
        TS_ASSERT_EQUALS(solution_exists, false);

        AbstractLinearParabolicPde<2,2>* p_pde = pde_and_bc.GetPde();
        TS_ASSERT_EQUALS(p_pde, &pde);

        // Set mCurrentSolution
        std::vector<double> data(10);
        for (unsigned i=0; i<10; i++)
        {
            data[i] = i + 0.45;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        // Test mCurrentSolution has been correctly set
        Vec solution =  pde_and_bc.GetSolution();
        ReplicatableVector solution_repl(solution);

        TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-12);
        }

        PetscInt size_of_solution = 0;
        VecGetSize(pde_and_bc.GetSolution(), &size_of_solution);
        TS_ASSERT_EQUALS(size_of_solution, 10);
    }

    void TestWithBoundaryConditionVaryingInSpace() throw(Exception)
    {
        // Create a ParabolicPdeAndBoundaryConditions object with spatially varying boundary condition
        HeatEquation<2> pde;
        FunctionalBoundaryCondition<2> functional_bc(&bc_func1);
        bool is_neumann_bc = false;

        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &functional_bc, is_neumann_bc);

        ChastePoint<2> point1;
        point1.rGetLocation()[0] = 0.0;
        point1.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point1), 0.0, 1e-6);

        ChastePoint<2> point2;
        point2.rGetLocation()[0] = 1.0;
        point2.rGetLocation()[1] = 5.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point2), 25.0, 1e-6);

        ChastePoint<2> point3;
        point3.rGetLocation()[0] = 3.0;
        point3.rGetLocation()[1] = -3.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point3), 9.0, 1e-6);
    }

    void TestWithBoundaryConditionVaryingInTime() throw(Exception)
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 2);

        // Create a ParabolicPdeAndBoundaryConditions object with time-dependent boundary condition
        HeatEquation<2> pde;
        FunctionalBoundaryCondition<2> functional_bc(&bc_func2);
        bool is_neumann_bc = false;

        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &functional_bc, is_neumann_bc);

        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        // At t=0, the boundary condition should take the value 1.0
        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 1.0, 1e-6);

        // At t=5, the boundary condition should take the value 3.5
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 3.5, 1e-6);

        // At t=10, the boundary condition should take the value 6.0
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 6.0, 1e-6);
    }

    void TestIn3d() throw(Exception)
    {
        // Create a 3D ParabolicPdeAndBoundaryConditions object
        HeatEquation<3> pde;
        ConstBoundaryCondition<3> bc(0.0);
        ParabolicPdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc);

        ChastePoint<3> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;
        point.rGetLocation()[2] = 1.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 0.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), true);
    }

    void TestArchivingWithoutSolution() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ParabolicPdeAndBoundaryConditions.arch";

        {
            // Create a ParabolicPdeAndBoundaryConditions object
            UniformSourceParabolicPde<2> pde(0.75);
            ConstBoundaryCondition<2> bc(2.45);
            bool is_neumann_bc = false;

            ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, is_neumann_bc);
            ParabolicPdeAndBoundaryConditions<2>* const p_const_pde_and_bc = &pde_and_bc;

            // Archive the object
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_const_pde_and_bc;
        }

        {
            ParabolicPdeAndBoundaryConditions<2>* p_pde_and_bc;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore object from the archive
            input_arch >> p_pde_and_bc;

            // Test that the object was archived correctly
            TS_ASSERT_EQUALS(p_pde_and_bc->IsNeumannBoundaryCondition(), false);

            ChastePoint<2> point;
            TS_ASSERT_DELTA(p_pde_and_bc->GetBoundaryCondition()->GetValue(point), 2.45, 1e-6);

            AbstractLinearParabolicPde<2,2>* p_pde = p_pde_and_bc->GetPde();
            TS_ASSERT(dynamic_cast<UniformSourceParabolicPde<2>*>(p_pde) != NULL);
            TS_ASSERT_DELTA(static_cast<UniformSourceParabolicPde<2>*>(p_pde)->GetCoefficient(), 0.75, 1e-6);

            // Avoid memory leaks
            delete p_pde_and_bc;
        }
    }

    void TestArchivingWithSolution() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ParabolicPdeAndBoundaryConditions.arch";

        {
            // Create a ParabolicPdeAndBoundaryConditions object
            UniformSourceParabolicPde<2> pde(0.75);
            ConstBoundaryCondition<2> bc(2.45);
            bool is_neumann_bc = false;

            ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, is_neumann_bc);
            pde_and_bc.SetDependentVariableName("quantity");

            std::vector<double> data(10);
            for (unsigned i=0; i<10; i++)
            {
                data[i] = i + 0.45;
            }

            Vec vector = PetscTools::CreateVec(data);
            pde_and_bc.SetSolution(vector);

            ParabolicPdeAndBoundaryConditions<2>* const p_const_pde_and_bc = &pde_and_bc;

            // Archive the object
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_const_pde_and_bc;
        }

        {
            ParabolicPdeAndBoundaryConditions<2>* p_pde_and_bc;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore object from the archive
            input_arch >> p_pde_and_bc;

            // Test that the object was archived correctly
            TS_ASSERT_EQUALS(p_pde_and_bc->IsNeumannBoundaryCondition(), false);
            TS_ASSERT_EQUALS(p_pde_and_bc->rGetDependentVariableName(), "quantity");

            ChastePoint<2> point;
            TS_ASSERT_DELTA(p_pde_and_bc->GetBoundaryCondition()->GetValue(point), 2.45, 1e-6);

            AbstractLinearParabolicPde<2,2>* p_pde = p_pde_and_bc->GetPde();
            TS_ASSERT(dynamic_cast<UniformSourceParabolicPde<2>*>(p_pde) != NULL);
            TS_ASSERT_DELTA(static_cast<UniformSourceParabolicPde<2>*>(p_pde)->GetCoefficient(), 0.75, 1e-6);

            Vec solution = p_pde_and_bc->GetSolution();
            ReplicatableVector solution_repl(solution);

            TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
            for (unsigned i=0; i<10; i++)
            {
                TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-6);
            }

            // Avoid memory leaks
            delete p_pde_and_bc;
        }
    }
};

#endif /* TESTPARABOLICPDEANDBOUNDARYCONDITIONS_HPP_ */
