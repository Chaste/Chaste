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

#include "AbstractPdeModifier.hpp"
#include "VtkMeshWriter.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
AbstractPdeModifier<DIM>::AbstractPdeModifier(boost::shared_ptr<PdeAndBoundaryConditions<DIM> > pPdeAndBcs)
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpPdeAndBcs(pPdeAndBcs),
      mSolution(NULL),
      mOutputDirectory(""),
      mOutputGradient(false)
{
    assert(DIM == 2);
}

template<unsigned DIM>
AbstractPdeModifier<DIM>::~AbstractPdeModifier()
{
}

template<unsigned DIM>
const boost::shared_ptr<PdeAndBoundaryConditions<DIM> > AbstractPdeModifier<DIM>::GetPdeAndBcs() const
{
    return mpPdeAndBcs;
}

template<unsigned DIM>
void AbstractPdeModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
#ifdef CHASTE_VTK
    if (DIM > 1)
    {
        std::ostringstream time_string;
        time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
        std::string results_file = "pde_results_" + mpPdeAndBcs->rGetDependentVariableName() + "_" + time_string.str();
        VtkMeshWriter<DIM,DIM>* p_vtk_mesh_writer = new VtkMeshWriter<DIM,DIM>(mOutputDirectory, results_file, false);

        ReplicatableVector solution_repl(mSolution);
        std::vector<double> pde_solution;
        for (unsigned i=0; i<mpFeMesh->GetNumNodes(); i++)
        {
           pde_solution.push_back(solution_repl[i]);
        }

        p_vtk_mesh_writer->AddPointData(mpPdeAndBcs->rGetDependentVariableName(), pde_solution);

        p_vtk_mesh_writer->WriteFilesUsingMesh(*mpFeMesh);
        delete p_vtk_mesh_writer;
    }
#endif //CHASTE_VTK
}

template<unsigned DIM>
void AbstractPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Cache the output directory
    this->mOutputDirectory = outputDirectory;
}

template<unsigned DIM>
bool AbstractPdeModifier<DIM>::GetOutputGradient()
{
    return mOutputGradient;
}

template<unsigned DIM>
void AbstractPdeModifier<DIM>::SetOutputGradient(bool outputGradient)
{
    mOutputGradient = outputGradient;
}

template<unsigned DIM>
void AbstractPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractPdeModifier<1>;
template class AbstractPdeModifier<2>;
template class AbstractPdeModifier<3>;
