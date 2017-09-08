/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef ABSTRACTPDESYSTEMMODIFIER_HPP_
#define ABSTRACTPDESYSTEMMODIFIER_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearPdeSystem.hpp"
#include "AbstractBoundaryCondition.hpp"
#include "VtkMeshWriter.hpp"
#include "ReplicatableVector.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde.hpp"

/**
 * An abstract modifier class containing functionality common to AbstractBoxDomainPdeSystemModifier,
 * AbstractGrowingDomainPdeSystemModifier and their subclasses, which solve a linear elliptic or
 * parabolic PDE coupled to a cell-based simulation.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class AbstractPdeSystemModifier : public AbstractCellBasedSimulationModifier<DIM>
{
protected:

    /**
     * Shared pointer to a linear PDE object.
     */
    boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > mpPdeSystem;

    /**
     * Shared pointer to a boundary condition object.
     */
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > mpBoundaryConditions;

    /**
     * Whether the boundary condition is Neumann (false corresponds to a Dirichlet boundary condition).
     *
     * \todo Generalize to allow mixed boundary conditions
     */
    bool mIsNeumannBoundaryCondition;

    /**
     * For use in PDEs where we wish to record a name for each dependent variable, e.g. oxygen concentration.
     */
    std::vector<std::string> mDependentVariableNames;

    /** The solution to the PDE problem at the current time step. */
    Vec mSolution;

    /** Pointer to the finite element mesh on which to solve the PDE. */
    TetrahedralMesh<DIM,DIM>* mpFeMesh;

    /** Store the output directory name. */
    std::string mOutputDirectory;

    /** Whether or not to calculate and output the gradient of the solution. */
    bool mOutputGradient;

    /**
     * Whether to output the PDE solution at each node of the FE mesh at output time steps.
     * Defaults to false.
     */
    bool mOutputSolutionAtPdeNodes;

    /** File that the values of the PDE solution are written out to. */
    out_stream mpVizPdeSolutionResultsFile;

    /**
     * Whether to delete the finite element mesh when we are destroyed.
     */
    bool mDeleteFeMesh;

public:

    /**
     * Constructor.
     *
     * @param pPdeSystem A shared pointer to a linear PDE system object (defaults to NULL)
     * @param pBoundaryConditions A vector of shared pointers to abstract boundary conditions
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param solution solution vector (defaults to NULL)
     */
    AbstractPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > pPdeSystem=nullptr,
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~AbstractPdeSystemModifier();

    /**
     * @return mpPdeSystem
     */
    boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > GetPdeSystem();

    /**
     * @return boundary conditions for the dependent variable with given index.
     *
     * @param pdeIndex the index (defaults to 0)
     */
    boost::shared_ptr<AbstractBoundaryCondition<DIM> > GetBoundaryCondition(unsigned pdeIndex=0);

    /**
     * @return mIsNeumannBoundaryCondition
     */
    bool IsNeumannBoundaryCondition();

    /**
     * Set the name of the dependent variable with given index.
     *
     * @param rName the name
     * @param pdeIndex the index (defaults to 0)
     */
    void SetDependentVariableName(const std::string& rName, unsigned pdeIndex=0);

    /**
     * Get the name of the dependent variable with given index.
     *
     * @param pdeIndex the index (defaults to 0)
     * @return the name
     */
    std::string& rGetDependentVariableName(unsigned pdeIndex=0);

    /**
     * @return whether the PDE has an averaged source
     */
    bool HasAveragedSourcePde();

    /**
     * In the case where the PDE has an averaged source, set the source terms
     * using the information in the given mesh.
     *
     * @param pMesh Pointer to a tetrahedral mesh
     * @param pCellPdeElementMap map between cells and elements
     */
    void SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM, DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);

    /**
     * @return mSolution.
     */
    Vec GetSolution();

    /**
     * @return mSolution (used in archiving)
     */
    Vec GetSolution() const;

    /**
     * @return mpFeMesh.
     */
    TetrahedralMesh<DIM,DIM>* GetFeMesh() const;

    /**
     * Overridden SetupSolve() method.
     *
     * Set mOutputDirectory and, if mOutputSolutionAtPdeNodes is set to true, open mpVizPdeSolutionResultsFile.
     * This method is overridden in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * As this method is pure virtual, it must be overridden in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)=0;

    /**
     * Overridden UpdateAtEndOfOutputTimeStep() method,
     * after UpdateAtEndOfTimeStep() has been called.
     *
     * Output the solution to the PDE at each cell to VTK and, if mOutputSolutionAtPdeNodes is set to true,
     * output the solution to the PDE at each node of mpFeMesh to mpVizPdeSolutionResultsFile.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden UpdateAtEndOfSolve() method.
     *
     * If mOutputSolutionAtPdeNodes is set to true, close mpVizPdeSolutionResultsFile.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Set whether to calculate and save the gradient of the solution to CellData.
     *
     * @return mOutputGradient
     */
    bool GetOutputGradient();

    /**
     * Set whether to calculate and save the gradient of the solution to CellData.
     *
     * @param outputGradient whether to output the gradient
     */
    void SetOutputGradient(bool outputGradient);

    /**
     * Set mOutputSolutionAtPdeNodes.
     *
     * @param outputSolutionAtPdeNodes whether to output the PDE solution at each node of the FE mesh at output time steps
     */
    void SetOutputSolutionAtPdeNodes(bool outputSolutionAtPdeNodes);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     *
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

/*
 * As this class is templated over PROBLEM_DIM, we put the implementation
 * in the header file to avoid explicit instantiation.
 */

template<unsigned DIM, unsigned PROBLEM_DIM>
AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::AbstractPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    Vec solution)
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpPdeSystem(pPdeSystem),
      mpBoundaryConditions(pBoundaryConditions),
      mIsNeumannBoundaryCondition(isNeumannBoundaryCondition),
      mSolution(nullptr),
      mOutputDirectory(""),
      mOutputGradient(false),
      mOutputSolutionAtPdeNodes(false),
      mDeleteFeMesh(false)
{
    if (solution)
    {
        mSolution = solution;
    }
    mDependentVariableNames.clear();
    mDependentVariableNames.resize(PROBLEM_DIM);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::~AbstractPdeSystemModifier()
{
    if (mDeleteFeMesh and mpFeMesh!=nullptr)
    {
        delete mpFeMesh;
    }
    if (mSolution)
    {
        PetscTools::Destroy(mSolution);
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::GetPdeSystem()
{
    return mpPdeSystem;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
boost::shared_ptr<AbstractBoundaryCondition<DIM> > AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::GetBoundaryCondition(unsigned pdeIndex)
{
    return mpBoundaryConditions[pdeIndex];
}

template<unsigned DIM, unsigned PROBLEM_DIM>
bool AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::IsNeumannBoundaryCondition()
{
    return mIsNeumannBoundaryCondition;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::SetDependentVariableName(const std::string& rName, unsigned pdeIndex)
{
    assert(pdeIndex < mDependentVariableNames.size());
    mDependentVariableNames[pdeIndex] = rName;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
std::string& AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::rGetDependentVariableName(unsigned pdeIndex)
{
    assert(pdeIndex < mDependentVariableNames.size());
    return mDependentVariableNames[pdeIndex];
}

template<unsigned DIM, unsigned PROBLEM_DIM>
bool AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::HasAveragedSourcePde()
{
    return ((boost::dynamic_pointer_cast<AveragedSourceEllipticPde<DIM> >(mpPdeSystem) != nullptr) ||
            (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<DIM> >(mpPdeSystem) != nullptr));
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap)
{
    assert(HasAveragedSourcePde());
    if (boost::dynamic_pointer_cast<AveragedSourceEllipticPde<DIM> >(mpPdeSystem) != nullptr)
    {
        boost::static_pointer_cast<AveragedSourceEllipticPde<DIM> >(mpPdeSystem)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
    }
    else if (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<DIM> >(mpPdeSystem) != nullptr)
    {
        boost::static_pointer_cast<AveragedSourceParabolicPde<DIM> >(mpPdeSystem)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
Vec AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
Vec AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::GetSolution() const
{
    return mSolution;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
TetrahedralMesh<DIM,DIM>* AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::GetFeMesh() const
{
    return mpFeMesh;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    // Cache the output directory
    this->mOutputDirectory = outputDirectory;

    if (mOutputSolutionAtPdeNodes)
    {
        if (PetscTools::AmMaster())
        {
            OutputFileHandler output_file_handler(outputDirectory+"/", false);
            mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizpdesolution");
        }
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (mOutputSolutionAtPdeNodes)
    {
        if (PetscTools::AmMaster())
        {
            (*mpVizPdeSolutionResultsFile) << SimulationTime::Instance()->GetTime() << "\t";

            assert(mpFeMesh != nullptr);
            assert(mDependentVariableNames[0] != "");

            for (unsigned i=0; i<mpFeMesh->GetNumNodes(); i++)
            {
                (*mpVizPdeSolutionResultsFile) << i << " ";
                const c_vector<double,DIM>& r_location = mpFeMesh->GetNode(i)->rGetLocation();
                for (unsigned k=0; k<DIM; k++)
                {
                    (*mpVizPdeSolutionResultsFile) << r_location[k] << " ";
                }

                assert(mSolution != nullptr);
                ReplicatableVector solution_repl(mSolution);
                (*mpVizPdeSolutionResultsFile) << solution_repl[i] << " ";
            }

            (*mpVizPdeSolutionResultsFile) << "\n";
        }
    }
#ifdef CHASTE_VTK
    if (DIM > 1)
    {
        std::ostringstream time_string;
        time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
        std::string results_file = "pde_results_" + mDependentVariableNames[0] + "_" + time_string.str();
        VtkMeshWriter<DIM,DIM>* p_vtk_mesh_writer = new VtkMeshWriter<DIM,DIM>(mOutputDirectory, results_file, false);

        ReplicatableVector solution_repl(mSolution);
        std::vector<double> pde_solution;
        for (unsigned i=0; i<mpFeMesh->GetNumNodes(); i++)
        {
           pde_solution.push_back(solution_repl[i]);
        }

        p_vtk_mesh_writer->AddPointData(mDependentVariableNames[0], pde_solution);

        p_vtk_mesh_writer->WriteFilesUsingMesh(*mpFeMesh);
        delete p_vtk_mesh_writer;
    }
#endif //CHASTE_VTK
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (mOutputSolutionAtPdeNodes)
    {
        if (PetscTools::AmMaster())
        {
            mpVizPdeSolutionResultsFile->close();
        }
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
bool AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::GetOutputGradient()
{
    return mOutputGradient;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::SetOutputGradient(bool outputGradient)
{
    mOutputGradient = outputGradient;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::SetOutputSolutionAtPdeNodes(bool outputSolutionAtPdeNodes)
{
    mOutputSolutionAtPdeNodes = outputSolutionAtPdeNodes;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*ABSTRACTPDESYSTEMMODIFIER_HPP_*/
