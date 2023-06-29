/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "DivisionBiasTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

template<unsigned DIM>
DivisionBiasTrackingModifier<DIM>::DivisionBiasTrackingModifier(c_vector<double, DIM> divisionBiasVector)
    : AbstractCellBasedSimulationModifier<DIM>(),
      mDivisionBiasVector(divisionBiasVector)
{
    // mDivisionBiasVector must be a unit vector
    assert(fabs(norm_2(mDivisionBiasVector) - 1.0) < 1e-6);
}

template<unsigned DIM>
DivisionBiasTrackingModifier<DIM>::~DivisionBiasTrackingModifier()
{
}

template<unsigned DIM>
const c_vector<double, DIM>& DivisionBiasTrackingModifier<DIM>::rGetDivisionBiasVector() const
{
    return mDivisionBiasVector;
}

template<unsigned DIM>
void DivisionBiasTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DivisionBiasTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DivisionBiasTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    /**
     * This hack is needed because in the case of a MeshBasedCellPopulation in which
     * multiple cell divisions have occurred over one time step, the Voronoi tessellation
     * (while existing) is out-of-date. Thus, if we did not regenerate the Voronoi
     * tessellation here, an assertion may trip as we try to access a Voronoi element
     * whose index exceeds the number of elements in the out-of-date tessellation.
     *
     * \todo work out how to properly fix this (#1986)
     */
    if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    }

    // Find centroid of cell population
    c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

    /**
     * Iterate over cell population and store each cell's signed distance along mDivisionBiasVector
     * through the centroid of the cell population, where zero corresponds to a cell located at the 
     * centroid of the cell population.
     */
    std::vector<double> biases;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        double distance = inner_prod(rCellPopulation.GetLocationOfCellCentre(*cell_iter) - centroid, mDivisionBiasVector);
        biases.push_back(distance);
    }

    // Map signed distances into [0,1]
    double min_distance = *std::min_element(biases.begin(), biases.end());
    double max_distance = *std::max_element(biases.begin(), biases.end());
    for (unsigned i = 0; i < biases.size(); i++)
    {
        biases[i] = (biases[i] - min_distance) / (max_distance - min_distance);
    }

    // Iterate over cell population
    unsigned i = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Store the cell's volume in CellData
        cell_iter->GetCellData()->SetItem("bias", biases[i]);
        i++;
    }
}

template<unsigned DIM>
void DivisionBiasTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DivisionBiasVector>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mDivisionBiasVector[index] << ",";
    }
    *rParamsFile << mDivisionBiasVector[DIM-1] << "</DivisionBiasVector>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DivisionBiasTrackingModifier<1>;
template class DivisionBiasTrackingModifier<2>;
template class DivisionBiasTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DivisionBiasTrackingModifier)

