/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "ImmersedBoundaryMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                   std::vector<ImmersedBoundaryElement<ELEMENT_DIM,SPACE_DIM>*> elements,
                                                                   unsigned numGridPtsX,
                                                                   unsigned numGridPtsY,
                                                                   unsigned membraneIndex)
        : mNumGridPtsX(numGridPtsX),
          mNumGridPtsY(numGridPtsY),
          mMembraneIndex(membraneIndex)
{
    // Clear mNodes and mElements
    Clear();

    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);

    switch (SPACE_DIM)
    {
        case 2:
            m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
            break;

        case 3:
            EXCEPTION("Not implemented yet in 3D");
            break;

        default:
            NEVER_REACHED;
    }

    // If the membrane index is UINT_MAX, there is no membrane; if not, there is
    mMeshHasMembrane = mMembraneIndex != UINT_MAX;

    // Populate mNodes and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index=0; elem_index<elements.size(); elem_index++)
    {
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_temp_element = elements[elem_index];
        mElements.push_back(p_temp_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[index];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_index=0; node_index<num_nodes_in_element; node_index++)
        {
            p_element->GetNode(node_index)->AddElement(element_index);
        }
    }

    // Set characteristic node spacing to the average distance between nodes
    double total_perimeter = 0.0;
    unsigned total_nodes = 0;
    for (unsigned elem_index = 0 ; elem_index < elements.size() ; elem_index++)
    {
        if (elem_index != mMembraneIndex)
        {
            total_perimeter += this->GetSurfaceAreaOfElement(elem_index);
            total_nodes += mElements[elem_index]->GetNumNodes();
        }
    }
    mCharacteristicNodeSpacing = total_perimeter / double(total_nodes);

    // Position source nodes at centroid of each cell, and set 'radius' (strength) to zero
    for (unsigned elem_it = 0 ; elem_it < elements.size() ; elem_it++)
    {
        unsigned this_elem_idx = mElements[elem_it]->GetIndex();

        // Each element other than the membrane element will have a source associated with it
        if (this_elem_idx != mMembraneIndex)
        {
            unsigned source_idx = mElementFluidSources.size();
            c_vector<double, SPACE_DIM> source_location = this->GetCentroidOfElement(this_elem_idx);

            mElementFluidSources.push_back(new FluidSource<SPACE_DIM>(source_idx, source_location));

            // Set source parameters
            mElementFluidSources.back()->SetAssociatedElementIndex(this_elem_idx);

            mElementFluidSources.back()->SetStrength(0.0);
        }
    }

    mElementFluidSources[2]->SetStrength(5.0 * 1e6);

    /*
     * Set up a number of sources to balance any active sources associated with elements
     */
    double balancing_source_spacing = 4.0 / (double)numGridPtsX;

    // We start 1/2 a grid-spacing in from the left-hand end, and place a source every 4-grid-spacings
    double current_location = balancing_source_spacing / 8.0;

    while (current_location < 1.0)
    {
        // Create a new fluid source at the current x-location and zero y-location
        unsigned source_idx = mBalancingFluidSources.size();
        mBalancingFluidSources.push_back(new FluidSource<SPACE_DIM>(source_idx, current_location));

        // Increment the current location
        current_location += balancing_source_spacing;
    }

    this->mMeshChangesDuringSimulation = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElongationShapeFactorOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    double discriminant = sqrt((moments(0) - moments(1))*(moments(0) - moments(1)) + 4.0*moments(2)*moments(2));

    // Note that as the matrix of second moments of area is symmetric, both its eigenvalues are real
    double largest_eigenvalue = (moments(0) + moments(1) + discriminant)*0.5;
    double smallest_eigenvalue = (moments(0) + moments(1) - discriminant)*0.5;

    double elongation_shape_factor = sqrt(largest_eigenvalue/smallest_eigenvalue);
    return elongation_shape_factor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryMesh()
{
    this->mMeshChangesDuringSimulation = false;
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::~ImmersedBoundaryMesh()
{
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Delete elements
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete nodes
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetCharacteristicNodeSpacing() const
{
    return mCharacteristicNodeSpacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetSpacingRatio() const
{
    ///todo If we ever permit mNumGridPtsX != mNumGridPtsY, need to decide how SpacingRatio is defined
    return mCharacteristicNodeSpacing / (1.0 / double(mNumGridPtsX));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumGridPtsX() const
{
    return mNumGridPtsX;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumGridPtsY() const
{
    return mNumGridPtsY;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNumGridPtsX(unsigned mesh_points_x)
{
    mNumGridPtsX = mesh_points_x;
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNumGridPtsY(unsigned mesh_points_y)
{
    mNumGridPtsY = mesh_points_y;
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNumGridPtsXAndY(unsigned numGridPts)
{
    mNumGridPtsX = numGridPts;
    mNumGridPtsY = numGridPts;
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetCharacteristicNodeSpacing(double node_spacing)
{
    mCharacteristicNodeSpacing = node_spacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetMembraneIndex(unsigned membrane_index)
{
    mMembraneIndex = membrane_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetMembraneElement()
{
    if (mMembraneIndex < UINT_MAX)
    {
        return this->GetElement(mMembraneIndex);
    }
    else
    {
        return NULL;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetMembraneIndex()
{
    return mMembraneIndex;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<FluidSource<SPACE_DIM>*>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetElementFluidSources()
{
    return mElementFluidSources;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<FluidSource<SPACE_DIM>*>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetBalancingFluidSources()
{
    return mBalancingFluidSources;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const multi_array<double, 3>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGet2dVelocityGrids() const
{
    return m2dVelocityGrids;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const multi_array<double, 4>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGet3dVelocityGrids() const
{
    return m3dVelocityGrids;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
multi_array<double, 3>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetModifiable2dVelocityGrids()
{
    return m2dVelocityGrids;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
multi_array<double, 4>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetModifiable3dVelocityGrids()
{
    return m3dVelocityGrids;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<Node<SPACE_DIM>*>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodes()
{
    return AbstractMesh<ELEMENT_DIM, SPACE_DIM>::mNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2)
{
    assert(SPACE_DIM == 2);

    // This code currently assumes the grid is precisely [0,1)x[0,1)
    c_vector<double, SPACE_DIM> vector = rLocation2 - rLocation1;

    /*
     * Handle the periodic condition here: if the points are more
     * than 0.5 apart in either direction, choose -(1.0-dist).
     */
    if (fabs(vector[0]) > 0.5)
    {
        vector[0] = copysign(fabs(vector[0]) - 1.0, -vector[0]);
    }
    if (fabs(vector[1]) > 0.5)
    {
        vector[1] = copysign(fabs(vector[1]) - 1.0, -vector[1]);
    }
    return vector;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement(unsigned index)
{
    // Only implemented in 2D
    assert(SPACE_DIM == 2);

    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    // The membrane must be treated differently
    if (index == mMembraneIndex)
    {
        return zero_vector<double>(2);
    }

    else
    {
        unsigned num_nodes = p_element->GetNumNodes();
        c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

        double centroid_x = 0;
        double centroid_y = 0;

        // Note that we cannot use GetVolumeOfElement() below as it returns the absolute, rather than signed, area
        double element_signed_area = 0.0;

        // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
        c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

        // Loop over vertices
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
            c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            double signed_area_term = this_x * next_y - this_y * next_x;

            centroid_x += (this_x + next_x) * signed_area_term;
            centroid_y += (this_y + next_y) * signed_area_term;
            element_signed_area += 0.5 * signed_area_term;

            pos_1 = pos_2;
        }

        assert(element_signed_area != 0.0);

        // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
        centroid = first_node_location;
        centroid[0] += centroid_x / (6.0 * element_signed_area);
        centroid[1] += centroid_y / (6.0 * element_signed_area);

        centroid[0] = centroid[0] < 0 ? centroid[0] + 1.0 : fmod(centroid[0], 1.0);
        centroid[1] = centroid[1] < 0 ? centroid[1] + 1.0 : fmod(centroid[1], 1.0);

        return centroid;
    }
}


/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void ImmersedBoundaryMesh<1,1>::ConstructFromMeshReader(AbstractMeshReader<1,1>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void ImmersedBoundaryMesh<1,2>::ConstructFromMeshReader(AbstractMeshReader<1,2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void ImmersedBoundaryMesh<1,3>::ConstructFromMeshReader(AbstractMeshReader<1,3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void ImmersedBoundaryMesh<2,3>::ConstructFromMeshReader(AbstractMeshReader<2,3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void ImmersedBoundaryMesh<2,2>::ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    ImmersedBoundaryMeshReader<2,2>& rIBMeshReader = dynamic_cast<ImmersedBoundaryMeshReader<2,2>&>(rMeshReader);

    assert(rIBMeshReader.HasNodePermutation() == false);
    // Store numbers of nodes and elements
    unsigned num_nodes = rIBMeshReader.GetNumNodes();
    unsigned num_elements = rIBMeshReader.GetNumElements();
    this->mCharacteristicNodeSpacing = rIBMeshReader.GetCharacteristicNodeSpacing();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rIBMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i=0; i<num_nodes; i++)
    {
        node_data = rIBMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool) node_data[2];
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(i, node_data, is_boundary_node));
    }

    rIBMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rIBMeshReader.GetNumElements());

    // Initially ensure there is no boundary element - this will be updated in the next loop if there is
    this->mMembraneIndex = UINT_MAX;
    this->mMeshHasMembrane = false;

    // Add elements
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        // Get the data for this element
        ImmersedBoundaryElementData element_data = rIBMeshReader.GetNextImmersedBoundaryElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        ImmersedBoundaryElement<2,2>* p_element = new ImmersedBoundaryElement<2,2>(elem_index, nodes);
        mElements.push_back(p_element);

        if (element_data.MembraneElement)
        {
            this->mMeshHasMembrane = true;
            this->mMembraneIndex = elem_index;
        }

        if (rIBMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rIBMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned) element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    // Get grid dimensions from grid file and set up grids accordingly
    this->mNumGridPtsX = rIBMeshReader.GetNumGridPtsX();
    this->mNumGridPtsY = rIBMeshReader.GetNumGridPtsY();
    m2dVelocityGrids.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);

    // Construct the velocity grids from mesh reader
    for (unsigned dim = 0 ; dim < 2 ; dim++)
    {
        for (unsigned grid_row = 0; grid_row < mNumGridPtsY; grid_row++)
        {
            std::vector<double> next_row = rIBMeshReader.GetNextGridRow();
            assert(next_row.size() == mNumGridPtsX);

            for (unsigned i = 0; i < mNumGridPtsX; i++)
            {
                m2dVelocityGrids[dim][i][grid_row] = next_row[i];
            }
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void ImmersedBoundaryMesh<3,3>::ConstructFromMeshReader(AbstractMeshReader<3,3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    // Get pointer to this element
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double element_volume = 0.0;

    // We map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
    c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
    c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

    unsigned num_nodes = p_element->GetNumNodes();
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index+1)%num_nodes);
        c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

        double this_x = pos_1[0];
        double this_y = pos_1[1];
        double next_x = pos_2[0];
        double next_y = pos_2[1];

        element_volume += 0.5*(this_x*next_y - next_x*this_y);

        pos_1 = pos_2;
    }

    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_volume);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    // Get pointer to this element
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double surface_area = 0.0;
    unsigned num_nodes = p_element->GetNumNodes();
    unsigned this_node_index = p_element->GetNodeGlobalIndex(0);
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        unsigned next_node_index = p_element->GetNodeGlobalIndex((local_index+1)%num_nodes);

        surface_area += this->GetDistanceBetweenNodes(this_node_index, next_node_index);
        this_node_index = next_node_index;
    }

    return surface_area;
}


//////////////////////////////////////////////////////////////////////
//                        2D-specific methods                       //
//////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMomentsOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    // Define helper variables
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes = p_element->GetNumNodes();
    c_vector<double, 3> moments = zero_vector<double>(3);

    // Since we compute I_xx, I_yy and I_xy about the centroid, we must shift each vertex accordingly
    c_vector<double, SPACE_DIM> centroid = GetCentroidOfElement(index);

    c_vector<double, SPACE_DIM> this_node_location = p_element->GetNodeLocation(0);
    c_vector<double, SPACE_DIM> pos_1 = this->GetVectorFromAtoB(centroid, this_node_location);

    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        unsigned next_index = (local_index+1)%num_nodes;
        c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation(next_index);
        c_vector<double, SPACE_DIM> pos_2 = this->GetVectorFromAtoB(centroid, next_node_location);

        double signed_area_term = pos_1(0)*pos_2(1) - pos_2(0)*pos_1(1);
        // Ixx
        moments(0) += (pos_1(1)*pos_1(1) + pos_1(1)*pos_2(1) + pos_2(1)*pos_2(1) ) * signed_area_term;

        // Iyy
        moments(1) += (pos_1(0)*pos_1(0) + pos_1(0)*pos_2(0) + pos_2(0)*pos_2(0)) * signed_area_term;

        // Ixy
        moments(2) += (pos_1(0)*pos_2(1) + 2*pos_1(0)*pos_1(1) + 2*pos_2(0)*pos_2(1) + pos_2(0)*pos_1(1)) * signed_area_term;

        pos_1 = pos_2;
    }

    moments(0) /= 12;
    moments(1) /= 12;
    moments(2) /= 24;

    /*
     * If the nodes owned by the element were supplied in a clockwise rather
     * than anticlockwise manner, or if this arose as a result of enforcing
     * periodicity, then our computed quantities will be the wrong sign, so
     * we need to fix this.
     */
    if (moments(0) < 0.0)
    {
        moments(0) = -moments(0);
        moments(1) = -moments(1);
        moments(2) = -moments(2);
    }
    return moments;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetShortAxisOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);

    // Calculate the moments of the element about its centroid (recall that I_xx and I_yy must be non-negative)
    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    // If the principal moments are equal...
    double discriminant = (moments(0) - moments(1))*(moments(0) - moments(1)) + 4.0*moments(2)*moments(2);
    if (fabs(discriminant) < 1e-10) ///\todo remove magic number? (see #1884 and #2401)
    {
        // ...then every axis through the centroid is a principal axis, so return a random unit vector
        short_axis(0) = RandomNumberGenerator::Instance()->ranf();
        short_axis(1) = sqrt(1.0 - short_axis(0)*short_axis(0));
    }
    else
    {
        // If the product of inertia is zero, then the coordinate axes are the principal axes
        if (moments(2) == 0.0)
        {
            if (moments(0) < moments(1))
            {
                short_axis(0) = 0.0;
                short_axis(1) = 1.0;
            }
            else
            {
                short_axis(0) = 1.0;
                short_axis(1) = 0.0;
            }
        }
        else
        {
            // Otherwise we find the eigenvector of the inertia matrix corresponding to the largest eigenvalue
            double lambda = 0.5*(moments(0) + moments(1) + sqrt(discriminant));

            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - lambda)/moments(2);

            double magnitude = norm_2(short_axis);
            short_axis = short_axis / magnitude;
        }
    }

    return short_axis;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class ImmersedBoundaryMesh<1,1>;
template class ImmersedBoundaryMesh<1,2>;
template class ImmersedBoundaryMesh<1,3>;
template class ImmersedBoundaryMesh<2,2>;
template class ImmersedBoundaryMesh<2,3>;
template class ImmersedBoundaryMesh<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ImmersedBoundaryMesh)
