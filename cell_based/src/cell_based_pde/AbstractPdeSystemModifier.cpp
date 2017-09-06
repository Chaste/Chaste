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

#include "AbstractPdeSystemModifier.hpp"
#include "VtkMeshWriter.hpp"
#include "ReplicatableVector.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde.hpp"

template<unsigned DIM, unsigned PROBLEM_DIM>
AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::AbstractPdeSystemModifier(boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem,
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
}

template<unsigned DIM, unsigned PROBLEM_DIM>
AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::~AbstractPdeSystemModifier()
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
boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::GetPdeSystem()
{
    return mpPdeSystem;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
boost::shared_ptr<AbstractBoundaryCondition<DIM> > AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::GetBoundaryCondition()
{
    return mpBoundaryConditions[0];
}

template<unsigned DIM, unsigned PROBLEM_DIM>
bool AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::IsNeumannBoundaryCondition()
{
    return mIsNeumannBoundaryCondition;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::SetDependentVariableName(const std::string& rName, unsigned pdeIndex)
{
    assert(pdeIndex < mDependentVariableNames.size());
    mDependentVariableNames[pdeIndex] = rName;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
std::string& AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::rGetDependentVariableName(unsigned pdeIndex)
{
    return mDependentVariableNames[pdeIndex];
}

template<unsigned DIM, unsigned PROBLEM_DIM>
bool AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::HasAveragedSourcePde()
{
    return ((boost::dynamic_pointer_cast<AveragedSourceEllipticPde<DIM> >(mpPdeSystem) != nullptr) ||
            (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<DIM> >(mpPdeSystem) != nullptr));
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap)
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
Vec AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::GetSolution()
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
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
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
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
bool AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::GetOutputGradient()
{
    return mOutputGradient;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::SetOutputGradient(bool outputGradient)
{
    mOutputGradient = outputGradient;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::SetOutputSolutionAtPdeNodes(bool outputSolutionAtPdeNodes)
{
    mOutputSolutionAtPdeNodes = outputSolutionAtPdeNodes;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}
