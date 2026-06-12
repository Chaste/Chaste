/*

Copyright (c) 2005-2026, University of Oxford.
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

#include <string>
#include "AbstractCardiacCellFactory.hpp"
#include "BidomainProblem.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "LuoRudy1991BackwardEulerOpt.hpp"
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

    // Problem settings (formerly HeartConfig)
    double mOdeTimeStep;
    double mPdeTimeStep;
    double mPrintingTimeStep;
    double mSimulationDuration;
    std::string mOutputDirectory;
    std::string mOutputFilenamePrefix;
    std::string mKspSolver;
    std::string mKspPreconditioner;
    double mKspAbsoluteTolerance;

public:

    MultiMeshSolver(double meshWidth, unsigned numMeshes):
        mMeshWidth(meshWidth),
        mNumMeshes(numMeshes),
        mOdeTimeStep(0.01),
        mPdeTimeStep(0.01),
        mPrintingTimeStep(0.1),
        mSimulationDuration(0.1),
        mOutputDirectory("BenchmarkMeshIndependence"),
        mOutputFilenamePrefix("Results"),
        mKspSolver("cg"),
        mKspPreconditioner("bjacobi"),
        mKspAbsoluteTolerance(1e-10)
    {
    }

    void SetOutputFilenamePrefix(const std::string& prefix) { mOutputFilenamePrefix = prefix; }
    void SetKspPreconditioner(const std::string& pc)        { mKspPreconditioner = pc; }

    void Solve()
    {
        // Loop over all the values of h requested
        for (unsigned mesh_index=0; mesh_index<mNumMeshes; mesh_index++)
        {
            CuboidMeshConstructor<DIM> constructor;
            DistributedTetrahedralMesh<DIM, DIM> mesh;
            constructor.Construct(mesh, mesh_index, mMeshWidth);

            unsigned num_ele_across = (unsigned) pow(2, mesh_index+1); // number of elements in each dimension
            GeneralPlaneStimulusCellFactory<CELL, DIM> cell_factory(num_ele_across, constructor.GetWidth());

            CARDIAC_PROBLEM cardiac_problem(&cell_factory);
            cardiac_problem.SetMesh(&mesh);
            cardiac_problem.SetOdePdeAndPrintingTimeSteps(mOdeTimeStep, mPdeTimeStep, mPrintingTimeStep);
            cardiac_problem.SetSimulationDuration(mSimulationDuration);
            cardiac_problem.SetOutputDirectory(mOutputDirectory);
            cardiac_problem.SetOutputFilenamePrefix(mOutputFilenamePrefix);
            cardiac_problem.SetKspSolverType(mKspSolver);
            cardiac_problem.SetKspPreconditionerType(mKspPreconditioner);
            cardiac_problem.SetKspAbsoluteTolerance(mKspAbsoluteTolerance);

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
        // Settings are now applied to the problem in MultiMeshSolver::Solve()
        // KSP solver type, tolerance, timesteps, etc. are set in the tester constructor defaults.
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
        MultiMeshSolver<CellLuoRudy1991FromCellMLBackwardEulerOpt, BidomainProblem<3>, 3, 2> tester(mesh_size, num_meshes);
        tester.SetKspPreconditioner("bjacobi");
        tester.SetOutputFilenamePrefix("BidomainMeshIndependencePEBJ");
        tester.Solve();
    }

    void TestMeshIndependentPreconditionersBD()
    {
        SetParametersMeshIndependent();
        MultiMeshSolver<CellLuoRudy1991FromCellMLBackwardEulerOpt, BidomainProblem<3>, 3, 2> tester(mesh_size, num_meshes);
        tester.SetKspPreconditioner("blockdiagonal");
        tester.SetOutputFilenamePrefix("BidomainMeshIndependencePEBD");
        tester.Solve();
    }

    void TestMeshIndependentPreconditionersLDU()
    {
        SetParametersMeshIndependent();
        MultiMeshSolver<CellLuoRudy1991FromCellMLBackwardEulerOpt, BidomainProblem<3>, 3, 2> tester(mesh_size, num_meshes);
        tester.SetKspPreconditioner("ldufactorisation");
        tester.SetOutputFilenamePrefix("BidomainMeshIndependencePELDU");
        tester.Solve();
    }
};


#endif /*TESTBENCHMARKMESHINDEPENDENCE_HPP_*/
