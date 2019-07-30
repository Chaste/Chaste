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
#ifndef LINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_
#define LINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_

#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "CvodeAdaptor.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "Warnings.hpp"
#include "VtkMeshWriter.hpp"

#include <boost/shared_ptr.hpp>

/**
 * A class for solving systems of parabolic PDEs and ODEs, which may be coupled
 * via their source terms:
 *
 * d/dt (u_i) = div (D(x) grad (u_i)) + f_i (x, u_1, ..., u_p, v_1, ..., v_q),  i=1,...,p,
 * d/dt (v_j) = g_j(x, u_1, ..., u_p, v_1, ..., v_q),  j=1,...,q.
 *
 * The solver class is templated over spatial dimension and PDE problem dimension (p).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM, unsigned PROBLEM_DIM=1>
class LinearParabolicPdeSystemWithCoupledOdeSystemSolver
    : public AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>,
      public AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    /** Pointer to the mesh. */
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** The PDE system to be solved. */
    AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpPdeSystem;

    /** Vector of pointers to ODE systems, defined at nodes. */
    std::vector<AbstractOdeSystemForCoupledPdeSystem*> mOdeSystemsAtNodes;

    /** The values of the ODE system state variables, interpolated at a quadrature point. */
    std::vector<double> mInterpolatedOdeStateVariables;

    /** The ODE solver. */
    boost::shared_ptr<AbstractIvpOdeSolver> mpOdeSolver;

    /**
     * A sampling timestep for writing results to file. Set to
     * PdeSimulationTime::GetPdeTimeStep() in the constructor;
     * may be overwritten using the SetSamplingTimeStep() method.
     */
    double mSamplingTimeStep;

    /** Whether ODE systems are present (if not, then the system comprises coupled PDEs only). */
    bool mOdeSystemsPresent;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

    /**
     * Whether the output directory should be cleared before solve or not. False by default.
     * Can be changed when setting the output directory
     */
    bool mClearOutputDirectory;

    /**
     * Write the current results to mpVtkMetaFile.
     */
    void WriteVtkResultsToFile();

    /**
     * @return the term to be added to the element stiffness matrix.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * @return the term to be added to the element stiffness vector.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Reset the member variable mInterpolatedOdeStateVariables.
     */
    void ResetInterpolatedQuantities();

    /**
     * Update the member variable mInterpolatedOdeStateVariables by computing the
     * interpolated value of each ODE state variable at each Gauss point.
     *
     * @param phiI
     * @param pNode pointer to a Node
     */
    void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode);

    /**
     * Initialise method: sets up the linear system (using the mesh to
     * determine the number of unknowns per row to preallocate) if it is not
     * already set up. Can use an initial solution as PETSc template,
     * or base it on the mesh size.
     *
     * @param initialSolution Initial solution (defaults to NULL) for PETSc to use as a template.
     */
    void InitialiseForSolve(Vec initialSolution=NULL);

    /**
     * Completely set up the linear system that has to be solved each timestep.
     *
     * @param currentSolution The current solution which can be used in setting up
     *  the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *   (mainly for dynamic solves).
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPdeSystem pointer to the PDE system
     * @param pBoundaryConditions pointer to the boundary conditions.
     * @param odeSystemsAtNodes optional vector of pointers to ODE systems, defined at nodes
     * @param pOdeSolver optional pointer to an ODE solver (defaults to NULL)
     */
    LinearParabolicPdeSystemWithCoupledOdeSystemSolver(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                                       AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
                                                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
                                                       std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystemsAtNodes=std::vector<AbstractOdeSystemForCoupledPdeSystem*>(),
                                                       boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver=boost::shared_ptr<AbstractIvpOdeSolver>());

    /**
     * Destructor.
     * If an ODE system is present, the pointers to the ODE system objects are deleted here.
     */
    ~LinearParabolicPdeSystemWithCoupledOdeSystemSolver();

    /**
     * Overridden PrepareForSetupLinearSystem() method.
     * Pass the current solution to the PDE system to the ODE system and solve it over the next timestep.
     *
     * @param currentPdeSolution the solution to the PDE system at the current time
     */
    void PrepareForSetupLinearSystem(Vec currentPdeSolution);

    /**
     * Set mOutputDirectory.
     *
     * @param outputDirectory the output directory to use
     * @param clearDirectory whether to clear outputDirectory or not. Note that the actual clearing happens when you call SolveAndWriteResultsToFile().
     *                       False by default.
     */
    void SetOutputDirectory(std::string outputDirectory, bool clearDirectory=false);

    /**
     * Set mSamplingTimeStep.
     *
     * @param samplingTimeStep the sampling timestep to use
     */
    void SetSamplingTimeStep(double samplingTimeStep);

    /**
     * Solve the coupled PDE/ODE system over the pre-specified time interval,
     * and record results using mSamplingTimeStep.
     */
    void SolveAndWriteResultsToFile();

    /**
     * Write the solution to VTK. Called by SolveAndWriteResultsToFile().
     *
     * @param solution the solution of the coupled PDE/ODE system
     * @param numTimeStepsElapsed the number of timesteps that have elapsed
     */
    void WriteVtkResultsToFile(Vec solution, unsigned numTimeStepsElapsed);

    /**
     * Get a pointer to the ODE system defined at a given node.
     *
     * @param index the global index of a node in the mpMesh
     * @return mOdeSystemsAtNodes[index]
     */
    AbstractOdeSystemForCoupledPdeSystem* GetOdeSystemAtNode(unsigned index);
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeMatrixTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double,PROBLEM_DIM>& rU,
    c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    double timestep_inverse = PdeSimulationTime::GetPdeTimeStepInverse();
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> matrix_term = zero_matrix<double>(PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate matrix_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double this_dudt_coefficient = mpPdeSystem->ComputeDuDtCoefficientFunction(rX, pde_index);
        c_matrix<double, SPACE_DIM, SPACE_DIM> this_pde_diffusion_term = mpPdeSystem->ComputeDiffusionTerm(rX, pde_index, pElement);
        c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> this_stiffness_matrix =
            prod(trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(this_pde_diffusion_term, rGradPhi)) )
                + timestep_inverse * this_dudt_coefficient * outer_prod(rPhi, rPhi);

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                matrix_term(i*PROBLEM_DIM + pde_index, j*PROBLEM_DIM + pde_index) = this_stiffness_matrix(i,j);
            }
        }
    }
    return matrix_term;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeVectorTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double,PROBLEM_DIM>& rU,
    c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    double timestep_inverse = PdeSimulationTime::GetPdeTimeStepInverse();
    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> vector_term;
    vector_term = zero_vector<double>(PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate vector_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double this_dudt_coefficient = mpPdeSystem->ComputeDuDtCoefficientFunction(rX, pde_index);
        double this_source_term = mpPdeSystem->ComputeSourceTerm(rX, rU, mInterpolatedOdeStateVariables, pde_index);
        c_vector<double, ELEMENT_DIM+1> this_vector_term;
        this_vector_term = (this_source_term + timestep_inverse*this_dudt_coefficient*rU(pde_index))* rPhi;

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            vector_term(i*PROBLEM_DIM + pde_index) = this_vector_term(i);
        }
    }

    return vector_term;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ResetInterpolatedQuantities()
{
    mInterpolatedOdeStateVariables.clear();

    if (mOdeSystemsPresent)
    {
        unsigned num_state_variables = mOdeSystemsAtNodes[0]->GetNumberOfStateVariables();
        mInterpolatedOdeStateVariables.resize(num_state_variables, 0.0);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
{
    if (mOdeSystemsPresent)
    {
        unsigned num_state_variables = mOdeSystemsAtNodes[0]->GetNumberOfStateVariables();

        for (unsigned i=0; i<num_state_variables; i++)
        {
            mInterpolatedOdeStateVariables[i] += phiI * mOdeSystemsAtNodes[pNode->GetIndex()]->rGetStateVariables()[i];
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem == NULL)
    {
        unsigned preallocation = mpMesh->CalculateMaximumContainingElementsPerProcess() + ELEMENT_DIM;
        if (ELEMENT_DIM > 1)
        {
            // Highest connectivity is closed
            preallocation--;
        }
        preallocation *= PROBLEM_DIM;

        /*
         * Use the current solution (ie the initial solution) as the
         * template in the alternative constructor of LinearSystem.
         * This is to avoid problems with VecScatter.
         */
        this->mpLinearSystem = new LinearSystem(initialSolution, preallocation);
    }

    assert(this->mpLinearSystem);
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetKspType("cg");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    this->SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::LinearParabolicPdeSystemWithCoupledOdeSystemSolver(
        TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
        AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
        BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystemsAtNodes,
        boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver)
    : AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>(pMesh, pBoundaryConditions),
      AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mpMesh(pMesh),
      mpPdeSystem(pPdeSystem),
      mOdeSystemsAtNodes(odeSystemsAtNodes),
      mpOdeSolver(pOdeSolver),
      mSamplingTimeStep(DOUBLE_UNSET),
      mOdeSystemsPresent(false),
      mClearOutputDirectory(false)
{
    this->mpBoundaryConditions = pBoundaryConditions;

    /*
     * If any ODE systems are passed in to the constructor, then we aren't just
     * solving a coupled PDE system, in which case the number of ODE system objects
     * must match the number of nodes in the finite element mesh.
     */
    if (!mOdeSystemsAtNodes.empty())
    {
        mOdeSystemsPresent = true;
        assert(mOdeSystemsAtNodes.size() == mpMesh->GetNumNodes());

        /*
         * In this case, if an ODE solver is not explicitly passed into the
         * constructor, then we create a default solver.
         */
        if (!mpOdeSolver)
        {
#ifdef CHASTE_CVODE
            mpOdeSolver.reset(new CvodeAdaptor);
#else
            mpOdeSolver.reset(new BackwardEulerIvpOdeSolver(mOdeSystemsAtNodes[0]->GetNumberOfStateVariables()));
#endif //CHASTE_CVODE
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~LinearParabolicPdeSystemWithCoupledOdeSystemSolver()
{
    if (mOdeSystemsPresent)
    {
        for (unsigned i=0; i<mOdeSystemsAtNodes.size(); i++)
        {
            delete mOdeSystemsAtNodes[i];
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PrepareForSetupLinearSystem(Vec currentPdeSolution)
{
    if (mOdeSystemsPresent)
    {
        double time = PdeSimulationTime::GetTime();
        double next_time = PdeSimulationTime::GetNextTime();
        double dt = PdeSimulationTime::GetPdeTimeStep();

        ReplicatableVector soln_repl(currentPdeSolution);
        std::vector<double> current_soln_this_node(PROBLEM_DIM);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mpMesh->GetNumNodes(); node_index++)
        {
            // Store the current solution to the PDE system at this node
            for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                double current_soln_this_pde_this_node = soln_repl[PROBLEM_DIM*node_index + pde_index];
                current_soln_this_node[pde_index] = current_soln_this_pde_this_node;
            }

            // Pass it into the ODE system at this node
            mOdeSystemsAtNodes[node_index]->SetPdeSolution(current_soln_this_node);

            // Solve ODE system at this node
            mpOdeSolver->SolveAndUpdateStateVariable(mOdeSystemsAtNodes[node_index], time, next_time, dt);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputDirectory(std::string outputDirectory, bool clearDirectory)
{
    mClearOutputDirectory = clearDirectory;
    this->mOutputDirectory = outputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetSamplingTimeStep(double samplingTimeStep)
{
    assert(samplingTimeStep >= this->mIdealTimeStep);
    mSamplingTimeStep = samplingTimeStep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SolveAndWriteResultsToFile()
{
    // A number of methods must have been called prior to this method
    if (this->mOutputDirectory == "")
    {
        EXCEPTION("SetOutputDirectory() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (this->mTimesSet == false)
    {
        EXCEPTION("SetTimes() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (this->mIdealTimeStep <= 0.0)
    {
        EXCEPTION("SetTimeStep() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (mSamplingTimeStep == DOUBLE_UNSET)
    {
        EXCEPTION("SetSamplingTimeStep() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (!this->mInitialCondition)
    {
        EXCEPTION("SetInitialCondition() must be called prior to SolveAndWriteResultsToFile()");
    }

#ifdef CHASTE_VTK
    // Create a .pvd output file
    OutputFileHandler output_file_handler(this->mOutputDirectory, mClearOutputDirectory);
    mpVtkMetaFile = output_file_handler.OpenOutputFile("results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";

    // Write initial condition to VTK
    Vec initial_condition = this->mInitialCondition;
    WriteVtkResultsToFile(initial_condition, 0);

    // The helper class TimeStepper deals with issues such as small final timesteps so we don't have to
    TimeStepper stepper(this->mTstart, this->mTend, mSamplingTimeStep);

    // Main time loop
    while (!stepper.IsTimeAtEnd())
    {
        // Reset start and end times
        this->SetTimes(stepper.GetTime(), stepper.GetNextTime());

        // Solve the system up to the new end time
        Vec soln = this->Solve();

        // Reset the initial condition for the next timestep
        if (this->mInitialCondition != initial_condition)
        {
            PetscTools::Destroy(this->mInitialCondition);
        }
        this->mInitialCondition = soln;

        // Move forward in time
        stepper.AdvanceOneTimeStep();

        // Write solution to VTK
        WriteVtkResultsToFile(soln, stepper.GetTotalTimeStepsTaken());
    }

    // Restore saved initial condition to avoid user confusion!
    if (this->mInitialCondition != initial_condition)
    {
        PetscTools::Destroy(this->mInitialCondition);
    }
    this->mInitialCondition = initial_condition;

    // Close .pvd output file
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#else //CHASTE_VTK
// LCOV_EXCL_START // We only test this in weekly builds
    WARNING("VTK is not installed and is required for this functionality");
// LCOV_EXCL_STOP
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WriteVtkResultsToFile(Vec solution, unsigned numTimeStepsElapsed)
{
#ifdef CHASTE_VTK

    // Create a new VTK file for this time step
    std::stringstream time;
    time << numTimeStepsElapsed;
    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(this->mOutputDirectory, "results_"+time.str(), false);

    /*
     * We first loop over PDEs. For each PDE we store the solution
     * at each node in a vector, then pass this vector to the mesh
     * writer.
     */
    ReplicatableVector solution_repl(solution);
    unsigned num_nodes = mpMesh->GetNumNodes();
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        // Store the solution of this PDE at each node
        std::vector<double> pde_index_data;
        pde_index_data.resize(num_nodes, 0.0);
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            pde_index_data[node_index] = solution_repl[PROBLEM_DIM*node_index + pde_index];
        }

        // Add this data to the mesh writer
        std::stringstream data_name;
        data_name << "PDE variable " << pde_index;
        mesh_writer.AddPointData(data_name.str(), pde_index_data);
    }

    if (mOdeSystemsPresent)
    {
        /*
         * We cannot loop over ODEs like PDEs, since the solutions are not
         * stored in one place. Therefore we build up a large 'vector of
         * vectors', then pass each component of this vector to the mesh
         * writer.
         */
        std::vector<std::vector<double> > ode_data;
        unsigned num_odes = mOdeSystemsAtNodes[0]->rGetStateVariables().size();
        ode_data.resize(num_odes);
        for (unsigned ode_index=0; ode_index<num_odes; ode_index++)
        {
            ode_data[ode_index].resize(num_nodes, 0.0);
        }

        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            std::vector<double> all_odes_this_node = mOdeSystemsAtNodes[node_index]->rGetStateVariables();
            for (unsigned i=0; i<num_odes; i++)
            {
                ode_data[i][node_index] = all_odes_this_node[i];
            }
        }

        for (unsigned ode_index=0; ode_index<num_odes; ode_index++)
        {
            std::vector<double> ode_index_data = ode_data[ode_index];

            // Add this data to the mesh writer
            std::stringstream data_name;
            data_name << "ODE variable " << ode_index;
            mesh_writer.AddPointData(data_name.str(), ode_index_data);
        }
    }

    mesh_writer.WriteFilesUsingMesh(*mpMesh);
    *mpVtkMetaFile << "        <DataSet timestep=\"";
    *mpVtkMetaFile << numTimeStepsElapsed;
    *mpVtkMetaFile << "\" group=\"\" part=\"0\" file=\"results_";
    *mpVtkMetaFile << numTimeStepsElapsed;
    *mpVtkMetaFile << ".vtu\"/>\n";
#endif // CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractOdeSystemForCoupledPdeSystem* LinearParabolicPdeSystemWithCoupledOdeSystemSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOdeSystemAtNode(unsigned index)
{
    return mOdeSystemsAtNodes[index];
}

#endif /*LINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_*/
