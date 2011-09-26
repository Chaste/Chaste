/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
void BuskeCompressionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                      AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This force class is defined for NodeBasedCellPopulations only
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL);

    NodeBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);


    c_vector<double, DIM> unit_vector;

    // Loop over cells in the population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Get the location of this node
        c_vector<double, DIM> node_i_location = rCellPopulation.GetNode(node_index)->rGetLocation();

        // Get the radius of this cell
        double radius_of_cell_i = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh().GetCellRadius(node_index);

        double delta_V_c = 0.0;
        c_vector<double, DIM> dVAdd_vector = zero_vector<double>(DIM);

        // Get the set of node indices corresponding to this cell's neighbours
        std::set<unsigned> neighbouring_node_indices = p_static_cast_cell_population->GetNeighbouringNodeIndices(node_index);

        // Loop over this set
        for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
             iter != neighbouring_node_indices.end();
             ++iter)
        {
            // Get the location of this node
            c_vector<double, DIM> node_j_location = rCellPopulation.GetNode(*iter)->rGetLocation();

            // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
            unit_vector = node_j_location - node_i_location;

            // Calculate the distance between the two nodes
            double dij = norm_2(unit_vector);

            unit_vector /= dij;

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh().GetCellRadius(*iter);

            // If the cells are close enough to exert a force on each other...
            if (dij < radius_of_cell_i + radius_of_cell_j)
            {
                // ...then compute the adhesion force and add it to the vector of forces...
				double xij = 0.5*(radius_of_cell_i*radius_of_cell_i - radius_of_cell_j*radius_of_cell_j + dij*dij)/dij;
				double dxijdd = 1.0 - xij/dij;
				double dVAdd = M_PI*dxijdd*(5.0*pow(radius_of_cell_i,2) + 3.0*pow(xij,2) - 8.0*radius_of_cell_i*xij)/3.0;

                dVAdd_vector += dVAdd*unit_vector;

                // ...and add the contribution to the compression force acting on cell i
                delta_V_c += M_PI*pow(radius_of_cell_i - xij,2)*(2*radius_of_cell_i - xij)/3.0;
            }
        }

        double V_A = 4.0/3.0*M_PI*pow(radius_of_cell_i,3) - delta_V_c;

        /**
         * Target volume of the cell
         * \todo Doesn't say in the Buske paper how they calculate this, so
         * we need to look at this to be sure it's what we want (#1764)
         */
        double V_T = 5.0;

        // Note: the sign in force_magnitude is different from the one in equation (A3) in the Buske paper
        rForces[node_index] += -mCompressionEnergyParameter/V_T*(V_T - V_A)*dVAdd_vector;
    }
}

template<unsigned DIM>
void BuskeCompressionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CompressionEnergyParameter>" << mCompressionEnergyParameter << "</CompressionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BuskeCompressionForce<1>;
template class BuskeCompressionForce<2>;
template class BuskeCompressionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeCompressionForce)
