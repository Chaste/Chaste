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

#include "BuskeCompressionForce.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
BuskeCompressionForce<DIM>::BuskeCompressionForce()
    : AbstractForce<DIM>(),
      mCompressionEnergyParameter(5.0)
{
}

template<unsigned DIM>
double BuskeCompressionForce<DIM>::GetCompressionEnergyParameter()
{
    return mCompressionEnergyParameter;
}

template<unsigned DIM>
void BuskeCompressionForce<DIM>::SetCompressionEnergyParameter(double compressionEnergyParameter)
{
    mCompressionEnergyParameter = compressionEnergyParameter;
}

template<unsigned DIM>
void BuskeCompressionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This force class is defined for NodeBasedCellPopulations only
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr);

    NodeBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

    c_vector<double, DIM> unit_vector;

    // Loop over cells in the population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        Node<DIM>* p_node_i = rCellPopulation.GetNode(node_index);

        // Get the location of this node
        const c_vector<double, DIM>& r_node_i_location = p_node_i->rGetLocation();

        // Get the radius of this cell
        double radius_of_cell_i = p_node_i->GetRadius();

        double delta_V_c = 0.0;
        c_vector<double, DIM> dVAdd_vector = zero_vector<double>(DIM);

        // Get the set of node indices corresponding to this cell's neighbours
        std::set<unsigned> neighbouring_node_indices = p_static_cast_cell_population->GetNeighbouringNodeIndices(node_index);

        // Loop over this set
        for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
             iter != neighbouring_node_indices.end();
             ++iter)
        {
            Node<DIM>* p_node_j = rCellPopulation.GetNode(*iter);

            // Get the location of this node
            const c_vector<double, DIM>& r_node_j_location = p_node_j->rGetLocation();

            // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
            unit_vector = r_node_j_location - r_node_i_location;

            // Calculate the distance between the two nodes
            double dij = norm_2(unit_vector);

            unit_vector /= dij;

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = p_node_j->GetRadius();

            // If the cells are close enough to exert a force on each other...
            if (dij < radius_of_cell_i + radius_of_cell_j)
            {
                // ...then compute the adhesion force and add it to the vector of forces...
                double xij = 0.5*(radius_of_cell_i*radius_of_cell_i - radius_of_cell_j*radius_of_cell_j + dij*dij)/dij;
                double dxijdd = 1.0 - xij/dij;
                double dVAdd = M_PI*dxijdd*(5.0*pow(radius_of_cell_i,2.0) + 3.0*pow(xij,2.0) - 8.0*radius_of_cell_i*xij)/3.0;

                dVAdd_vector += dVAdd*unit_vector;

                // ...and add the contribution to the compression force acting on cell i
                delta_V_c += M_PI*pow(radius_of_cell_i - xij,2.0)*(2*radius_of_cell_i - xij)/3.0;
            }
        }

        double V_A = 4.0/3.0*M_PI*pow(radius_of_cell_i,3.0) - delta_V_c;

        /**
         * Target volume of the cell
         * \todo Doesn't say in the Buske paper how they calculate this, so
         * we need to look at this to be sure it's what we want (#1764)
         */
        double V_T = 5.0;

        // Note: the sign in force_magnitude is different from the one in equation (A3) in the Buske paper
        c_vector<double, DIM> applied_force = -mCompressionEnergyParameter/V_T*(V_T - V_A)*dVAdd_vector;
        p_node_i->AddAppliedForceContribution(applied_force);
    }
}

template<unsigned DIM>
void BuskeCompressionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CompressionEnergyParameter>" << mCompressionEnergyParameter << "</CompressionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class BuskeCompressionForce<1>;
template class BuskeCompressionForce<2>;
template class BuskeCompressionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeCompressionForce)
