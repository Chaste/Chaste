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

#ifndef TESTBENCHMARKMESHINDEPENDENCE_HPP_
#define TESTBENCHMARKMESHINDEPENDENCE_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "BidomainProblem.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "CuboidMeshConstructor.hpp"

#include "PetscSetupAndFinalize.hpp"


static const unsigned num_meshes = 2; // Use 4 if needed to assess performance;
static const double mesh_size = 0.01; //0.07;

/**
 * MultiMeshSolver
 * Runs simulations over multiple meshes for a particular cell type, mono/bidomain and dimension
 */
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class MultiMeshSolver
{
private:
    double mMeshWidth;
    unsigned mNumMeshes;

    void DisplayRun(unsigned numMesh)
    {
        unsigned num_ele_across = (unsigned) pow(2, numMesh+1);// number of elements in each dimension
        double scaling = mMeshWidth/(double) num_ele_across;

        std::cout<<"================================================================================"<<std::endl;
        std::cout<<"Solving in "<<DIM<<"D\t";
        std::cout<<"Space step "<< scaling << " cm (mesh " << numMesh << ")" << "\n";
        std::cout<<"PDE step "<<HeartConfig::Instance()->GetPdeTimeStep()<<" ms"<<"\t";
        std::cout<<"ODE step "<<HeartConfig::Instance()->GetOdeTimeStep()<<" ms"<<"\t";
        if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
        {
            std::cout<<"KSP absolute "<<HeartConfig::Instance()->GetAbsoluteTolerance()<<"\t";
        }
        else
        {
            std::cout<<"KSP relative "<<HeartConfig::Instance()->GetRelativeTolerance()<<"\t";
        }
    }

public:

    MultiMeshSolver(double meshWidth, unsigned numMeshes):
        mMeshWidth(meshWidth),
        mNumMeshes(numMeshes)
    {
    }

    void Solve()
    {
        // Loop over all the values of h requested
        for (unsigned mesh_index=0; mesh_index<mNumMeshes; mesh_index++)
        {
            // DisplayRun(mesh_index);

            CuboidMeshConstructor<DIM> constructor;
            DistributedTetrahedralMesh<DIM, DIM> mesh;
            constructor.Construct(mesh, mesh_index, mMeshWidth);

            // Do I need a unique name for the output file???
            //HeartConfig::Instance()->SetOutputFilenamePrefix ("Results");

            unsigned num_ele_across = (unsigned) pow(2, mesh_index+1); // number of elements in each dimension
            GeneralPlaneStimulusCellFactory<CELL, DIM> cell_factory(num_ele_across, constructor.GetWidth());

            CARDIAC_PROBLEM cardiac_problem(&cell_factory);
            cardiac_problem.SetMesh(&mesh);

            cardiac_problem.Initialise();

            try
            {
                cardiac_problem.Solve();
            }
            catch(Exception& e)
            {
                std::cout << "Simulation threw an exception!" << std::endl;
                std::cout << e.GetMessage() << std::endl;
            }

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
        }
    }

};

class TestBenchmarkPreconditionersForMeshIndependence : public CxxTest::TestSuite
{
private:

    void SetParametersMeshIndependent()
    {
        HeartEventHandler::Reset();

        // Timesteps and simulation duration
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(0.1);  //ms

        // Output directory
        HeartConfig::Instance()->SetOutputDirectory("BenchmarkMeshIndependence");

        // Solver and preconditioner selection through Chaste parameter system.
        // Not all the possible methods can be selected via HeartConfig. If an error like the following
        // happens, see next comment.
        //    heart/test/TestRealisticLinearAlgebra.hpp:105: Error: Test failed:
        //    Chaste error: heart/src/problem/HeartConfig.cpp:866: Unknown solver type provided
        HeartConfig::Instance()->SetKSPSolver("cg");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-10);
        //HeartConfig::Instance()->SetUseRelativeTolerance(1e-12);

        // In the case you want to select a solver or preconditioner not supported in HeartConfig,
        // you should talk to the KSP object directly. Uncomment and modify accordingly
//        PetscTools::SetOption("-ksp_type", "bicg");
//        PetscTools::SetOption("-pc_type", "asm");

        // If extra parameters need to be passed to the solver/preconditioner (e.g. number of levels of
        // fill to use in ILU preconditioner), they can be added to parameters database in the following
        // way
//        PetscTools::SetOption("-pc_factor_levels", "3");

        // Traces KSP solution (# of iterations, residual, etc)
        //PetscTools::SetOption("-ksp_monitor", "");

        // Traces true (non-preconditioned) residual
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
        PetscTools::SetOption("-ksp_monitor_true_residual", "");
#else
        PetscTools::SetOption("-ksp_truemonitor", "");
#endif

        // Enables extra logging (# of flops, messages, reductions, etc)
//        PetscTools::SetOption("-log_summary", "");

        //PetscTools::SetOption("-options_table", "");
        PetscTools::SetOption("-ksp_norm_type", "unpreconditioned");
        PetscTools::SetOption("-ksp_max_it", "200");
    }

public:

    void TestMeshIndependentPreconditionersBJ()
    {
        SetParametersMeshIndependent();
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainMeshIndependencePEBJ");

        MultiMeshSolver<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester(mesh_size, num_meshes);

        tester.Solve();
    }

    void TestMeshIndependentPreconditionersBD()
    {
        SetParametersMeshIndependent();
        HeartConfig::Instance()->SetKSPPreconditioner("blockdiagonal");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainMeshIndependencePEBD");

        MultiMeshSolver<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester(mesh_size, num_meshes);

        tester.Solve();
    }

    void TestMeshIndependentPreconditionersLDU()
    {
        SetParametersMeshIndependent();
        HeartConfig::Instance()->SetKSPPreconditioner("ldufactorisation");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainMeshIndependencePELDU");

        MultiMeshSolver<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester(mesh_size, num_meshes);

        tester.Solve();
    }
};


#endif /*TESTBENCHMARKMESHINDEPENDENCE_HPP_*/
