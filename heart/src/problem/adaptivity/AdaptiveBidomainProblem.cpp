/*

Copyright (C) Fujitsu Laboratories of Europe, 2009-2010

*/

/*

Copyright (c) 2005-2013, University of Oxford.
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


#define COVERAGE_IGNORE /// \todo #2367 We no longer test adaptivity library

#ifdef CHASTE_ADAPTIVITY

#include "AdaptiveBidomainProblem.hpp"
#include "BidomainSolver.hpp"

#include "VtkMeshWriter.hpp"
#include "TetrahedralMesh.hpp"
#include "BidomainTissue.hpp"
#include "HeartRegionCodes.hpp"
#include "HeartConfig.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"
#include "ProgressReporter.hpp"
#include "PetscTools.hpp"
#include "RegularStimulus.hpp"

AdaptiveBidomainProblem::AdaptiveBidomainProblem(
            AbstractCardiacCellFactory<3>* pCellFactory, bool hasBath)
    : BidomainProblem<3>(pCellFactory, hasBath),
      mIsMeshAdapting(true),
      mInitializeFromVtu(false),
      mUseNeumannBoundaryConditions(false),
      mNeumannStimulusIndex(0),
      mNeumannStimulusLowerValue(DBL_MAX),
      mNeumannStimulusUpperValue(-DBL_MAX),
      mNeumannStimulusMagnitude(0.0),
      mNeumannStimulusDuration(0.0),
      mNeumannStimulusPeriod(DBL_MAX)
//      mGoodEdgeRange(0.0),
//      mBadEdgeCriterion(0.0)
{
    mFixedExtracellularPotentialNodes.resize(0);
    mpAdaptiveMesh = new AdaptiveTetrahedralMesh;
}

AdaptiveBidomainProblem::~AdaptiveBidomainProblem()
{
    delete mpAdaptiveMesh;
}

void AdaptiveBidomainProblem::DoNotAdaptMesh()
{
    mIsMeshAdapting = false;
}

//void AdaptiveBidomainProblem::SetAdaptCriterion(double range, double criterion)
//{
//    mGoodEdgeRange = range;
//    mBadEdgeCriterion = criterion;
//}

void AdaptiveBidomainProblem::AddCurrentSolutionToAdaptiveMesh( Vec solution )
{
    HeartEventHandler::BeginEvent(HeartEventHandler::USER1);

    ReplicatableVector replicatable_solution( solution );
    std::vector<double> vm_for_vtk, phi_for_vtk;
    vm_for_vtk.resize(mpMesh->GetNumNodes());
    phi_for_vtk.resize(mpMesh->GetNumNodes());

    for (AbstractTetrahedralMesh<3,3>::NodeIterator it=mpMesh->GetNodeIteratorBegin();
         it != mpMesh->GetNodeIteratorEnd();
         ++it)
    {
        vm_for_vtk[it->GetIndex()]  = replicatable_solution[2*it->GetIndex()];
        phi_for_vtk[it->GetIndex()] = replicatable_solution[2*it->GetIndex()+1];
    }
    mpAdaptiveMesh->AddPointData("Vm", vm_for_vtk);
    mpAdaptiveMesh->AddPointData("Phi", phi_for_vtk);

    std::vector<std::string> state_variable_names = mpBidomainTissue->GetCardiacCell(0)->rGetStateVariableNames();
    unsigned number_of_state_variables = state_variable_names.size();
    std::vector< std::vector<double> > state_variable_values;
    state_variable_values.resize( mpMesh->GetNumNodes() );

    // add state variable to vtu file as point data
    for (unsigned variable=0; variable<number_of_state_variables; variable++)
    {
        if (variable != mpBidomainTissue->GetCardiacCell(0)->GetVoltageIndex())
        {
            std::vector<double> variable_for_vtk;
            variable_for_vtk.resize(mpMesh->GetNumNodes());
            for (AbstractTetrahedralMesh<3,3>::NodeIterator it=mpMesh->GetNodeIteratorBegin();
                 it != mpMesh->GetNodeIteratorEnd();
                 ++it)
            {
                variable_for_vtk[it->GetIndex()] = mpBidomainTissue->GetCardiacCell(it->GetIndex())->GetStdVecStateVariables()[variable];
            }
            mpAdaptiveMesh->AddPointData(state_variable_names[variable], variable_for_vtk);
        }
    }

    HeartEventHandler::EndEvent(HeartEventHandler::USER1);
}

void AdaptiveBidomainProblem::InitializeSolutionOnAdaptedMesh( VtkMeshReader<3,3>* reader )
{
    HeartEventHandler::BeginEvent(HeartEventHandler::USER1);

    std::vector<double> adapted_vm, adapted_phi;

    reader->GetPointData("Vm", adapted_vm);
    reader->GetPointData("Phi", adapted_phi);

    Vec solution = mpMesh->GetDistributedVectorFactory()->CreateVec(2);
    DistributedVector nic = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(solution);
    std::vector<DistributedVector::Stripe> stripe;
    stripe.reserve(2);

    for (unsigned i=0; i<2; i++)
    {
        stripe.push_back(DistributedVector::Stripe(nic, i));
    }

    for (DistributedVector::Iterator it = nic.Begin();
             it != nic.End();
             ++it)
    {
        stripe[0][it] = adapted_vm[it.Global];
        stripe[1][it] = adapted_phi[it.Global];
    }

    nic.Restore();

    std::vector<std::string> state_variable_names = mpBidomainTissue->GetCardiacCell(0)->rGetStateVariableNames();
    unsigned number_of_state_variables = state_variable_names.size();

    for (unsigned variable=0; variable<number_of_state_variables; variable++)
    {
        if (variable != mpBidomainTissue->GetCardiacCell(0)->GetVoltageIndex())
        {
            std::vector<double> adapted_state_variable;
            reader->GetPointData( state_variable_names[variable], adapted_state_variable);

            for (DistributedVector::Iterator it = nic.Begin();
                     it != nic.End();
                     ++it)
            {
                mpBidomainTissue->GetCardiacCell(it.Global)->SetStateVariable(variable, adapted_state_variable[it.Global]);
            }
        }
        else
        {
            for (DistributedVector::Iterator it = nic.Begin();
                     it != nic.End();
                     ++it)
            {
                mpBidomainTissue->GetCardiacCell(it.Global)->SetStateVariable(variable, adapted_vm[it.Global]);
            }

        }
    }

    if (mSolution)
    {
        PetscTools::Destroy(mSolution);
    }

    mSolution = solution;

    HeartEventHandler::EndEvent(HeartEventHandler::USER1);
}

void AdaptiveBidomainProblem::AdaptMesh()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::USER1);
    mpAdaptiveMesh->AdaptMesh();
    HeartEventHandler::EndEvent(HeartEventHandler::USER1);

    if ( mpAdaptiveMesh->GetAdaptSuccess() )
    {
        if (mWriteInfo)
        {
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            std::cout << "Adapt completed. New mesh has " << mpAdaptiveMesh->GetNumNodes() << " nodes" << std::endl;
            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }
        // Adapt succeeded, need new mesh, boundary conditions, solver
        delete mpMesh;
        DistributedTetrahedralMesh<3,3>* p_new_mesh = new DistributedTetrahedralMesh<3,3>;
        VtkMeshReader<3,3> mesh_reader( mpAdaptiveMesh->GetVtkUnstructuredGrid() );
        HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
        p_new_mesh->ConstructFromMeshReader(mesh_reader);
        HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
        mpMesh = p_new_mesh;

        mpCellFactory->SetMesh( mpMesh );

        if (mUseNeumannBoundaryConditions)
        {
            SetupNeumannBoundaryConditionOnMesh();
        }
        else
        {
            boost::shared_ptr<BoundaryConditionsContainer<3, 3, 2> > p_new_bcc(new BoundaryConditionsContainer<3, 3, 2>);
            for (unsigned problem_index=0; problem_index<2; problem_index++)
            {
                p_new_bcc->DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
            }
            mpBoundaryConditionsContainer = p_new_bcc;
        }

        delete mpBidomainTissue;
        BidomainTissue<3>* new_bidomain_tissue;
        new_bidomain_tissue = new BidomainTissue<3>(mpCellFactory);
        mpBidomainTissue = new_bidomain_tissue;
        mpCardiacTissue = mpBidomainTissue;

        delete mpSolver;
        BidomainSolver<3,3>* p_new_solver;
        p_new_solver = new BidomainSolver<3,3>(false, mpMesh, mpBidomainTissue,
                                               mpBoundaryConditionsContainer.get());
        mpSolver = p_new_solver;
        mpSolver->SetTimeStep(HeartConfig::Instance()->GetPdeTimeStep());

        InitializeSolutionOnAdaptedMesh( &mesh_reader );
    }
    else
    {
        NEVER_REACHED;
    }
}

void AdaptiveBidomainProblem::SetNeumannStimulusMagnitudeAndDuration(double magnitude, double duration, double period)
{
    mNeumannStimulusMagnitude = magnitude;
    mNeumannStimulusDuration  = duration;
    mNeumannStimulusPeriod = std::min( period, HeartConfig::Instance()->GetSimulationDuration() );
}

void AdaptiveBidomainProblem::UseNeumannBoundaryCondition(unsigned index)
{
    mUseNeumannBoundaryConditions = true;
    mNeumannStimulusIndex = index;
}

double AdaptiveBidomainProblem::GetTargetError()
{
    return HeartConfig::Instance()->GetTargetErrorForAdaptivity();
}

double AdaptiveBidomainProblem::GetSigma()
{
    return HeartConfig::Instance()->GetSigmaForAdaptivity();
}

double AdaptiveBidomainProblem::GetMaxEdgeLength()
{
    return HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity();
}

double AdaptiveBidomainProblem::GetMinEdgeLength()
{
    return HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity();
}

double AdaptiveBidomainProblem::GetGradation()
{
    return HeartConfig::Instance()->GetGradationForAdaptivity();
}

unsigned AdaptiveBidomainProblem::GetMaxMeshNodes()
{
    return HeartConfig::Instance()->GetMaxNodesForAdaptivity();
}

unsigned AdaptiveBidomainProblem::GetNumAdaptSweeps()
{
    return HeartConfig::Instance()->GetNumberOfAdaptiveSweeps();
}

void AdaptiveBidomainProblem::SetupNeumannBoundaryConditionOnMesh()
{
    boost::shared_ptr<BoundaryConditionsContainer<3, 3, 2> > p_new_bcc(new BoundaryConditionsContainer<3, 3, 2>);

    mpNeumannStimulusBoundaryCondition = new StimulusBoundaryCondition<3>(mpNeumannStimulus);
    ConstBoundaryCondition<3>* p_zero_bc = new ConstBoundaryCondition<3>(0.0);

    // loop over boundary elements
    AbstractTetrahedralMesh<3,3>::BoundaryElementIterator iter;
    iter = mpMesh->GetBoundaryElementIteratorBegin();

    while (iter != mpMesh->GetBoundaryElementIteratorEnd())
    {
        double x = ((*iter)->CalculateCentroid())[mNeumannStimulusIndex];
        ///\todo remove magic number? (#1884)
        if ( (x-mNeumannStimulusLowerValue)*(x-mNeumannStimulusLowerValue) <= 1e-10 )
        {
            p_new_bcc->AddNeumannBoundaryCondition(*iter, mpNeumannStimulusBoundaryCondition);
        }
        iter++;
    }

    // Ground other end of domain
    for (AbstractTetrahedralMesh<3,3>::NodeIterator node_iter=mpMesh->GetNodeIteratorBegin();
         node_iter != mpMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        if (fabs((*node_iter).rGetLocation()[mNeumannStimulusIndex] - mNeumannStimulusUpperValue) < 1e-6)
        {
            p_new_bcc->AddDirichletBoundaryCondition(&(*node_iter), p_zero_bc, 1);
        }
    }

    mpBoundaryConditionsContainer = p_new_bcc;
}

void AdaptiveBidomainProblem::LoadSimulationFromVtuFile()
{
    mInitializeFromVtu = true;
    AbstractCardiacProblem<3,3,2>::Initialise();
}

void AdaptiveBidomainProblem::Solve()
{
    OutputFileHandler file_handler(HeartConfig::Instance()->GetOutputDirectory(), false);

    mOutputDirectory = file_handler.GetOutputDirectoryFullPath();
    mOutputFilenamePrefix = HeartConfig::Instance()->GetOutputFilenamePrefix();

    PreSolveChecks();

    mpNeumannStimulus = new RegularStimulus(mNeumannStimulusMagnitude, mNeumannStimulusDuration, mNeumannStimulusPeriod, 0.0);

    if (mUseNeumannBoundaryConditions)
    {
        // Determine the values a and b to apply Neumann bcs at x_i=a (stimulus), x_i=b (ground)
        double local_min = DBL_MAX;
        double local_max = -DBL_MAX;

        for (AbstractTetrahedralMesh<3,3>::NodeIterator iter=mpMesh->GetNodeIteratorBegin();
             iter != mpMesh->GetNodeIteratorEnd();
             ++iter)
        {
            double value = (*iter).rGetLocation()[mNeumannStimulusIndex];
            if(value < local_min)
            {
                local_min = value;
            }
            if(value > local_max)
            {
               local_max = value;
            }
        }

        int mpi_ret = MPI_Allreduce(&local_min, &mNeumannStimulusLowerValue, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
        assert(mpi_ret == MPI_SUCCESS);
        mpi_ret = MPI_Allreduce(&local_max, &mNeumannStimulusUpperValue, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
        assert(mpi_ret == MPI_SUCCESS);
        UNUSED_OPT(mpi_ret);
        SetupNeumannBoundaryConditionOnMesh();
    }
    else
    {
        if(mpBoundaryConditionsContainer == NULL) // the user didnt supply a bcc
        {
            // set up the default bcc
            boost::shared_ptr<BoundaryConditionsContainer<3, 3, 2> > p_allocated_memory(new BoundaryConditionsContainer<3, 3, 2>);
            mpDefaultBoundaryConditionsContainer = p_allocated_memory;
            for (unsigned problem_index=0; problem_index<2; problem_index++)
            {
                mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
            }
            mpBoundaryConditionsContainer = mpDefaultBoundaryConditionsContainer;
        }
    }

    mpSolver = new BidomainSolver<3,3>(false, mpMesh, mpBidomainTissue,
                                       mpBoundaryConditionsContainer.get());
    mSolution = CreateInitialCondition();

    TimeStepper stepper(0.0, HeartConfig::Instance()->GetSimulationDuration(),
                        HeartConfig::Instance()->GetPrintingTimeStep());

    std::string progress_reporter_dir;

    assert( !mPrintOutput );

    progress_reporter_dir = ""; // progress printed to CHASTE_TEST_OUTPUT

    // create a progress reporter so users can track how much has gone and
    // estimate how much time is left. (Note this has to be done after the
    // InitialiseWriter above (if mPrintOutput==true)
    ProgressReporter progress_reporter(progress_reporter_dir, 0.0, HeartConfig::Instance()->GetSimulationDuration());
    progress_reporter.Update(0);

    // Initialize adaptive mesh using Chaste mesh
//    mpAdaptiveMesh->ConstructFromMesh( mpMesh );
//    mpAdaptiveMesh->CalculateSENListAndSids();

    // Get initial condition from file, or add Chaste initial condition to adaptive mesh
    if ( mInitializeFromVtu )
    {
        mpAdaptiveMesh->ConstructFromVtuFile( HeartConfig::Instance()->GetMeshName() );
        mpAdaptiveMesh->CalculateSENListAndSids();

        VtkMeshReader<3,3> mesh_reader( HeartConfig::Instance()->GetMeshName() );

        InitializeSolutionOnAdaptedMesh( &mesh_reader );
    }
    else
    {
        mpAdaptiveMesh->ConstructFromMesh( mpMesh );
        mpAdaptiveMesh->CalculateSENListAndSids();
        AddCurrentSolutionToAdaptiveMesh( mSolution );
    }

    // Use printing time step as adaptive time step...

    unsigned count = 0;

    // Set up filenames for convenient ParaView visualization
    std::ostringstream adapt_file_name, solution_file_name;

//    {
//        std::cout << "Values of Vm at node 2,000,000" << std::endl;
//        ReplicatableVector replicatable_solution( mSolution );
//        std::cout << replicatable_solution[2*2e6] << std::endl;        // Vm at node x is mSolution[2*x] (phi is mSolution[2*x+1])
//    }

    mpSolver->SetTimeStep(HeartConfig::Instance()->GetPdeTimeStep());

    while ( !stepper.IsTimeAtEnd() )
    {
        {
            solution_file_name.str("");
            solution_file_name << mOutputFilenamePrefix << std::setw(4) << std::setfill('0') << count << ".vtu";
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            mpAdaptiveMesh->WriteMeshToFile( mOutputDirectory, solution_file_name.str() );
            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }


        // Adapt the mesh
        if ( mIsMeshAdapting && ( count > 0 || mInitializeFromVtu ) )
        {
            AdaptMesh();
        }

//        // For non-adapting heart simulation, we may want to refine the mesh *ONCE* after the end of the stimulus
        // (in order to change the resolution).
//        if ( (! mIsMeshAdapting) && (count == 1) )
//        {
//            AdaptMesh();
//        }

        // Solve from now up to the next printing time step
        mpSolver->SetTimes(stepper.GetTime(), stepper.GetNextTime());
        mpSolver->SetInitialCondition( mSolution );

        try
        {
            Vec new_solution = mpSolver->Solve();
            PetscTools::Destroy(mSolution);
            mSolution = new_solution;
//            ReplicatableVector replicatable_solution( mSolution );
//            std::cout << replicatable_solution[2*2e6] << std::endl;        // Vm at node x is mSolution[2*x]
        }
        catch (Exception &e)
        {
            // Free memory.

            delete mpSolver;
            mpSolver=NULL;
            if (!mUseNeumannBoundaryConditions)
            {
                mpDefaultBoundaryConditionsContainer = mpBoundaryConditionsContainer;
            }
            delete mpNeumannStimulus;

            PetscTools::ReplicateException(true);
            // Re-throw
            HeartEventHandler::Reset();//EndEvent(HeartEventHandler::EVERYTHING);

            CloseFilesAndPostProcess();
            throw e;
        }
        PetscTools::ReplicateException(false);

        // Add current solution into the node "point" data
        AddCurrentSolutionToAdaptiveMesh( mSolution );

        // update the current time
        stepper.AdvanceOneTimeStep();

        if (mWriteInfo)
        {
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            WriteInfo(stepper.GetTime());
            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }

        progress_reporter.Update(stepper.GetTime());

        OnEndOfTimestep(stepper.GetTime());

        count++;
    }

    {
        solution_file_name.str("");
        solution_file_name << mOutputFilenamePrefix << std::setw(4) << std::setfill('0') << count << ".vtu";
        HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
        mpAdaptiveMesh->WriteMeshToFile( mOutputDirectory, solution_file_name.str() );
        HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
    }

    // Free solver
    delete mpSolver;

    if (!mUseNeumannBoundaryConditions)
    {
        mpDefaultBoundaryConditionsContainer = mpBoundaryConditionsContainer;
    }
    delete mpNeumannStimulus;

    // close the file that stores voltage values
    progress_reporter.PrintFinalising();
    CloseFilesAndPostProcess();
    HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
}

#endif //CHASTE_ADAPTIVITY
#undef COVERAGE_IGNORE /// \todo #2367 We no longer test adaptivity library
