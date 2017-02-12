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

#include "FarhadifarForce.hpp"

template<unsigned DIM>
FarhadifarForce<DIM>::FarhadifarForce()
   : AbstractForce<DIM>(),
     mAreaElasticityParameter(1.0), // These parameters are Case I in Farhadifar's paper
     mPerimeterContractilityParameter(0.04),
     mLineTensionParameter(0.12),
     mBoundaryLineTensionParameter(0.12) // this parameter as such does not exist in Farhadifar's model.
{
}

template<unsigned DIM>
FarhadifarForce<DIM>::~FarhadifarForce()
{
}

template<unsigned DIM>
void FarhadifarForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL);

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    VertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = r_mesh.GetNumNodes();
    unsigned num_elements = r_mesh.GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = r_mesh.GetElementIteratorBegin();
         elem_iter != r_mesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = r_mesh.GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = r_mesh.GetSurfaceAreaOfElement(elem_index);

        // Get the target area of the cell associated with this element, throwing an exception if this is not stored as a CellData item
        try
        {
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarForce");
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        c_vector<double, DIM> area_elasticity_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> perimeter_contractility_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> line_tension_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = r_mesh.GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = r_mesh.GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's area elasticity (note the minus sign)
            c_vector<double, DIM> element_area_gradient = r_mesh.GetAreaGradientOfElementAtNode(p_element, local_index);
            area_elasticity_contribution -= GetAreaElasticityParameter()*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;
        }

        r_mesh.GetNode(node_index)->AddAppliedForceContribution(area_elasticity_contribution);
    }

    // Iterate over all edges and add line tension and perimeter contractility force contributions
    for (typename VertexMesh<DIM, DIM>::EdgeIterator edge_iter = r_mesh.EdgesBegin();
         edge_iter != r_mesh.EdgesEnd();
         ++edge_iter)
    {
    	// Compute the line tension parameter for this edge (this is reset to the value for a non-boundary edge below, if required)
        double line_tension_parameter = GetBoundaryLineTensionParameter();

        unsigned node_A_index = edge_iter.GetNodeA()->GetIndex();
        unsigned node_B_index = edge_iter.GetNodeB()->GetIndex();

        const c_vector<double, DIM>& r_node_A_location = edge_iter.GetNodeA()->rGetLocation();
        const c_vector<double, DIM>& r_node_B_location = edge_iter.GetNodeB()->rGetLocation();

        c_vector<double, DIM> edge_vector = r_mesh.GetVectorFromAtoB(r_node_B_location, r_node_A_location);
        double edge_length = norm_2(edge_vector);

        assert(edge_length > DBL_EPSILON);

        c_vector<double, DIM> edge_gradient = edge_vector/edge_length;

        // Perimeter contractility term
        unsigned elem_index = edge_iter.GetElemIndex();
        c_vector<double, DIM> perimeter_contractility_force_on_node_A = -GetPerimeterContractilityParameter()* element_perimeters[elem_index]*edge_gradient;
        c_vector<double, DIM> perimeter_contractility_force_on_node_B = -perimeter_contractility_force_on_node_A;
        r_mesh.GetNode(node_A_index)->AddAppliedForceContribution(perimeter_contractility_force_on_node_A);
        r_mesh.GetNode(node_B_index)->AddAppliedForceContribution(perimeter_contractility_force_on_node_B);

        unsigned other_elem_index = edge_iter.GetOtherElemIndex();
        if (other_elem_index != UINT_MAX)
        {
        	line_tension_parameter = GetLineTensionParameter();

            perimeter_contractility_force_on_node_A = -GetPerimeterContractilityParameter()* element_perimeters[other_elem_index]*edge_gradient;
            perimeter_contractility_force_on_node_B = -perimeter_contractility_force_on_node_A;
            r_mesh.GetNode(node_A_index)->AddAppliedForceContribution(perimeter_contractility_force_on_node_A);
            r_mesh.GetNode(node_B_index)->AddAppliedForceContribution(perimeter_contractility_force_on_node_B);
        }

        // Line tension term
        c_vector<double, DIM> line_tension_force_on_node_A = -line_tension_parameter*edge_gradient;
        c_vector<double, DIM> line_tension_force_on_node_B = -line_tension_force_on_node_A;
        r_mesh.GetNode(node_A_index)->AddAppliedForceContribution(line_tension_force_on_node_A);
        r_mesh.GetNode(node_B_index)->AddAppliedForceContribution(line_tension_force_on_node_B);
    }
}

template<unsigned DIM>
double FarhadifarForce<DIM>::GetAreaElasticityParameter()
{
    return mAreaElasticityParameter;
}

template<unsigned DIM>
double FarhadifarForce<DIM>::GetPerimeterContractilityParameter()
{
    return mPerimeterContractilityParameter;
}

template<unsigned DIM>
double FarhadifarForce<DIM>::GetLineTensionParameter()
{
    return mLineTensionParameter;
}

template<unsigned DIM>
double FarhadifarForce<DIM>::GetBoundaryLineTensionParameter()
{
    return mBoundaryLineTensionParameter;
}

template<unsigned DIM>
void FarhadifarForce<DIM>::SetAreaElasticityParameter(double areaElasticityParameter)
{
    mAreaElasticityParameter = areaElasticityParameter;
}

template<unsigned DIM>
void FarhadifarForce<DIM>::SetPerimeterContractilityParameter(double perimeterContractilityParameter)
{
    mPerimeterContractilityParameter = perimeterContractilityParameter;
}

template<unsigned DIM>
void FarhadifarForce<DIM>::SetLineTensionParameter(double lineTensionParameter)
{
    mLineTensionParameter = lineTensionParameter;
}

template<unsigned DIM>
void FarhadifarForce<DIM>::SetBoundaryLineTensionParameter(double boundaryLineTensionParameter)
{
    mBoundaryLineTensionParameter = boundaryLineTensionParameter;
}

template<unsigned DIM>
void FarhadifarForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AreaElasticityParameter>" << mAreaElasticityParameter << "</AreaElasticityParameter>\n";
    *rParamsFile << "\t\t\t<PerimeterContractilityParameter>" << mPerimeterContractilityParameter << "</PerimeterContractilityParameter>\n";
    *rParamsFile << "\t\t\t<LineTensionParameter>" << mLineTensionParameter << "</LineTensionParameter>\n";
    *rParamsFile << "\t\t\t<BoundaryLineTensionParameter>" << mBoundaryLineTensionParameter << "</BoundaryLineTensionParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class FarhadifarForce<1>;
template class FarhadifarForce<2>;
template class FarhadifarForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FarhadifarForce)
