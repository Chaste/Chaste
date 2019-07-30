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

#include "RadialCellDataDistributionWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::RadialCellDataDistributionWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("radial_dist.dat"),
      mVariableName(""),
      mNumRadialBins(UNSIGNED_UNSET)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{
    // Calculate the centre of the cell population
    c_vector<double, SPACE_DIM> centre = pCellPopulation->GetCentroidOfCellPopulation();

    // Calculate the distance between each node and the centre of the cell population, as well as the maximum of these
    std::map<double, CellPtr> radius_cell_map;
    double max_distance_from_centre = 0.0;
    for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        double distance = norm_2(pCellPopulation->GetLocationOfCellCentre(*cell_iter) - centre);
        radius_cell_map[distance] = *cell_iter;

        if (distance > max_distance_from_centre)
        {
            max_distance_from_centre = distance;
        }
    }

    // Create vector of radius intervals
    std::vector<double> radius_intervals;
    for (unsigned i=0; i<mNumRadialBins; i++)
    {
        double upper_radius = max_distance_from_centre*((double) i+1)/((double) mNumRadialBins);
        radius_intervals.push_back(upper_radius);
    }

    // Calculate PDE solution in each radial interval
    double lower_radius = 0.0;
    for (unsigned i=0; i<mNumRadialBins; i++)
    {
        unsigned counter = 0;
        double average_solution = 0.0;

        for (std::map<double, CellPtr>::iterator iter = radius_cell_map.begin(); iter != radius_cell_map.end(); ++iter)
        {
            if (iter->first > lower_radius && iter->first <= radius_intervals[i])
            {
                average_solution += (iter->second)->GetCellData()->GetItem(mVariableName);
                counter++;
            }
        }
        if (counter != 0)
        {
            average_solution /= (double) counter;
        }

        // Write results to file
        *this->mpOutStream << radius_intervals[i] << " " << average_solution << " ";
        lower_radius = radius_intervals[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Calculate the centre of the cell population
    c_vector<double, SPACE_DIM> centre = pCellPopulation->GetCentroidOfCellPopulation();

    // Calculate the distance between each node and the centre of the cell population, as well as the maximum of these
    std::map<double, CellPtr> radius_cell_map;
    double max_distance_from_centre = 0.0;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        double distance = norm_2(pCellPopulation->GetLocationOfCellCentre(*cell_iter) - centre);
        radius_cell_map[distance] = *cell_iter;

        if (distance > max_distance_from_centre)
        {
            max_distance_from_centre = distance;
        }
    }

    // Create vector of radius intervals
    std::vector<double> radius_intervals;
    for (unsigned i=0; i<mNumRadialBins; i++)
    {
        double upper_radius = max_distance_from_centre*((double) i+1)/((double) mNumRadialBins);
        radius_intervals.push_back(upper_radius);
    }

    // Calculate PDE solution in each radial interval
    double lower_radius = 0.0;
    for (unsigned i=0; i<mNumRadialBins; i++)
    {
        unsigned counter = 0;
        double average_solution = 0.0;

        for (std::map<double, CellPtr>::iterator iter = radius_cell_map.begin(); iter != radius_cell_map.end(); ++iter)
        {
            if (iter->first > lower_radius && iter->first <= radius_intervals[i])
            {
                average_solution += (iter->second)->GetCellData()->GetItem(mVariableName);
                counter++;
            }
        }
        if (counter != 0)
        {
            average_solution /= (double) counter;
        }

        // Write results to file
        *this->mpOutStream << radius_intervals[i] << " " << average_solution << " ";
        lower_radius = radius_intervals[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::SetVariableName(std::string variableName)
{
    mVariableName = variableName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::GetVariableName() const
{
    return mVariableName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::SetNumRadialBins(unsigned numRadialBins)
{
    mNumRadialBins = numRadialBins;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RadialCellDataDistributionWriter<ELEMENT_DIM, SPACE_DIM>::GetNumRadialBins() const
{
    return mNumRadialBins;
}

// Explicit instantiation
template class RadialCellDataDistributionWriter<1,1>;
template class RadialCellDataDistributionWriter<1,2>;
template class RadialCellDataDistributionWriter<2,2>;
template class RadialCellDataDistributionWriter<1,3>;
template class RadialCellDataDistributionWriter<2,3>;
template class RadialCellDataDistributionWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RadialCellDataDistributionWriter)

