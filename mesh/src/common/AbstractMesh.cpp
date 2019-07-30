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

#include "AbstractMesh.hpp"
#include "Exception.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::AbstractMesh()
        : mpDistributedVectorFactory(nullptr),
          mMeshFileBaseName(""),
          mMeshChangesDuringSimulation(false)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractMesh()
{
    // Iterate over nodes and free the memory
    for (unsigned i = 0; i < mNodes.size(); i++)
    {
        delete mNodes[i];
    }
    if (mpDistributedVectorFactory)
    {
        delete mpDistributedVectorFactory;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes() const
{
    return mBoundaryNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodeAttributes() const
{
    /* Note, this implementation assumes that all nodes have the same number of attributes
     * so that the first node in the container is indicative of the others.*/
    assert(mNodes.size() != 0u);
    return mNodes[0]->GetNumNodeAttributes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    unsigned local_index = SolveNodeMapping(index);
    return mNodes[local_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeOrHaloNode(unsigned index) const
{
    return GetNode(index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeFromPrePermutationIndex(unsigned index) const
{
    if (mNodePermutation.empty())
    {
        return GetNode(index);
    }
    else
    {
        return GetNode(mNodePermutation[index]);
    }
}

// LCOV_EXCL_START
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile)
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships()
{
    // Does nothing, since an AbstractMesh has no elements
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistributedVectorFactory* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetDistributedVectorFactory()
{
    if (mpDistributedVectorFactory == nullptr)
    {
        mpDistributedVectorFactory = new DistributedVectorFactory(GetNumNodes());
        if (PetscTools::IsParallel())
        {
            SetElementOwnerships(); // So any parallel implementation with shared mesh has owners set
        }
    }
    return mpDistributedVectorFactory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::SetDistributedVectorFactory(DistributedVectorFactory* pFactory)
{
    if (mpDistributedVectorFactory)
    {
        EXCEPTION("Cannot change the mesh's distributed vector factory once it has been set.");
    }
    if (pFactory->GetNumProcs() != PetscTools::GetNumProcs())
    {
        EXCEPTION("The distributed vector factory provided to the mesh is for the wrong number of processes.");
    }
    mpDistributedVectorFactory = pFactory;
}

// LCOV_EXCL_START
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorBegin() const
{
    return mBoundaryNodes.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorEnd() const
{
    return mBoundaryNodes.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName() const
{
    if (!IsMeshOnDisk())
    {
        EXCEPTION("This mesh was not constructed from a file.");
    }

    return mMeshFileBaseName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::IsMeshOnDisk() const
{
    return (mMeshFileBaseName != "");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodePermutation() const
{
    return mNodePermutation;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;
    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceBetweenNodes(unsigned indexA, unsigned indexB)
{
    c_vector<double, SPACE_DIM> vector = GetVectorFromAtoB(mNodes[indexA]->rGetLocation(),
                                                           mNodes[indexB]->rGetLocation());
    return norm_2(vector);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    return CalculateBoundingBox(mNodes).GetWidth(rDimension);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ChasteCuboid<SPACE_DIM> AbstractMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundingBox(const std::vector<Node<SPACE_DIM>*>& rNodes) const
{
    // Work out the max and min location in each co-ordinate direction.
    c_vector<double, SPACE_DIM> minimum_point;
    c_vector<double, SPACE_DIM> maximum_point;

    // Deal with the special case of no nodes by returning a cuboid with zero volume.
    if (rNodes.empty())
    {
        minimum_point = scalar_vector<double>(SPACE_DIM, 0.0);
        maximum_point = scalar_vector<double>(SPACE_DIM, 0.0);
    }
    else
    {
        minimum_point = scalar_vector<double>(SPACE_DIM, DBL_MAX);
        maximum_point = scalar_vector<double>(SPACE_DIM, -DBL_MAX);

        // Iterate through nodes
        /// \todo #1322 use a const version of NodeIterator here
        for (unsigned index = 0; index < rNodes.size(); index++)
        {
            if (!rNodes[index]->IsDeleted())
            {
                // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                c_vector<double, SPACE_DIM> position;
                position = rNodes[index]->rGetLocation();

                // Update max/min
                for (unsigned i = 0; i < SPACE_DIM; i++)
                {
                    if (position[i] < minimum_point[i])
                    {
                        minimum_point[i] = position[i];
                    }
                    if (position[i] > maximum_point[i])
                    {
                        maximum_point[i] = position[i];
                    }
                }
            }
        }
    }

    ChastePoint<SPACE_DIM> min(minimum_point);
    ChastePoint<SPACE_DIM> max(maximum_point);

    return ChasteCuboid<SPACE_DIM>(min, max);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ChasteCuboid<SPACE_DIM> AbstractMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundingBox() const
{
    ChasteCuboid<SPACE_DIM> bounding_cuboid = CalculateBoundingBox(mNodes);

    return bounding_cuboid;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestNodeIndex(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    if (mNodes.empty())
    {
        /*
         * This happens in parallel if a process isn't assigned any nodes.
         * This case is covered in TestDistributedTetrahedralMesh::TestConstructLinearMeshSmallest
         * but only when there are 3 or more processes
         */
        return UINT_MAX; // LCOV_EXCL_LINE
    }
    // Hold the best distance from node to point found so far
    // and the (local) node at which this was recorded
    unsigned best_node_index = 0u;
    double best_node_point_distance = DBL_MAX;

    const c_vector<double, SPACE_DIM>& test_location = rTestPoint.rGetLocation();
    // Now loop through the nodes, calculating the distance and updating best_node_point_distance
    for (unsigned node_index = 0; node_index < mNodes.size(); node_index++)
    {
        // Calculate the distance from the chosen point to the current node
        double node_point_distance = norm_2(GetVectorFromAtoB(mNodes[node_index]->rGetLocation(), test_location));
        // Update the "best" distance and node index if necessary
        if (node_point_distance < best_node_point_distance)
        {
            best_node_index = node_index;
            best_node_point_distance = node_point_distance;
        }
    }

    // Return the global index of the closest node to the point
    // (which differs from the local index "best_node_index" in parallel)
    // In the distributed case, we'll have to do an AllReduce next
    return mNodes[best_node_index]->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Scale(const double xScale, const double yScale, const double zScale)
{
    unsigned num_nodes = mNodes.size();

    for (unsigned i = 0; i < num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        if (SPACE_DIM >= 3)
        {
            r_location[2] *= zScale;
        }
        if (SPACE_DIM >= 2)
        {
            r_location[1] *= yScale;
        }
        r_location[0] *= xScale;
    }

    RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Translate(
    const double xMovement,
    const double yMovement,
    const double zMovement)
{
    c_vector<double, SPACE_DIM> displacement;

    switch (SPACE_DIM)
    {
        case 3:
            displacement[2] = zMovement;
        // purposeful fallthrough - Note the presence of this comment stops compiler warning for gcc >= 7!
        // See https://developers.redhat.com/blog/2017/03/10/wimplicit-fallthrough-in-gcc-7/
        case 2:
            displacement[1] = yMovement;
        // fallthrough on purpose
        case 1:
            displacement[0] = xMovement;
            // fallthrough on purpose
    }

    Translate(displacement);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Translate(const c_vector<double, SPACE_DIM>& rDisplacement)
{
    unsigned num_nodes = this->mNodes.size();

    for (unsigned i = 0; i < num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = this->mNodes[i]->rGetModifiableLocation();
        r_location += rDisplacement;
    }

    RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(c_matrix<double, SPACE_DIM, SPACE_DIM> rotationMatrix)
{
    unsigned num_nodes = this->mNodes.size();

    for (unsigned i = 0; i < num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = this->mNodes[i]->rGetModifiableLocation();
        r_location = prod(rotationMatrix, r_location);
    }

    RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(c_vector<double, 3> axis, double angle)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE
    double norm = norm_2(axis);
    c_vector<double, 3> unit_axis = axis / norm;

    c_matrix<double, SPACE_DIM, SPACE_DIM> rotation_matrix;

    double c = cos(angle);
    double s = sin(angle);

    rotation_matrix(0, 0) = unit_axis(0) * unit_axis(0) + c * (1 - unit_axis(0) * unit_axis(0));
    rotation_matrix(0, 1) = unit_axis(0) * unit_axis(1) * (1 - c) - unit_axis(2) * s;
    rotation_matrix(1, 0) = unit_axis(0) * unit_axis(1) * (1 - c) + unit_axis(2) * s;
    rotation_matrix(1, 1) = unit_axis(1) * unit_axis(1) + c * (1 - unit_axis(1) * unit_axis(1));
    rotation_matrix(0, 2) = unit_axis(0) * unit_axis(2) * (1 - c) + unit_axis(1) * s;
    rotation_matrix(1, 2) = unit_axis(1) * unit_axis(2) * (1 - c) - unit_axis(0) * s;
    rotation_matrix(2, 0) = unit_axis(0) * unit_axis(2) * (1 - c) - unit_axis(1) * s;
    rotation_matrix(2, 1) = unit_axis(1) * unit_axis(2) * (1 - c) + unit_axis(0) * s;
    rotation_matrix(2, 2) = unit_axis(2) * unit_axis(2) + c * (1 - unit_axis(2) * unit_axis(2));

    Rotate(rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::RotateX(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double, SPACE_DIM, SPACE_DIM> x_rotation_matrix = identity_matrix<double>(SPACE_DIM);

    x_rotation_matrix(1, 1) = cos(theta);
    x_rotation_matrix(1, 2) = sin(theta);
    x_rotation_matrix(2, 1) = -sin(theta);
    x_rotation_matrix(2, 2) = cos(theta);
    Rotate(x_rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::RotateY(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double, SPACE_DIM, SPACE_DIM> y_rotation_matrix = identity_matrix<double>(SPACE_DIM);

    y_rotation_matrix(0, 0) = cos(theta);
    y_rotation_matrix(0, 2) = -sin(theta);
    y_rotation_matrix(2, 0) = sin(theta);
    y_rotation_matrix(2, 2) = cos(theta);

    Rotate(y_rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::RotateZ(const double theta)
{
    if (SPACE_DIM < 2)
    {
        EXCEPTION("This rotation is not valid in less than 2D");
    }
    c_matrix<double, SPACE_DIM, SPACE_DIM> z_rotation_matrix = identity_matrix<double>(SPACE_DIM);

    z_rotation_matrix(0, 0) = cos(theta);
    z_rotation_matrix(0, 1) = sin(theta);
    z_rotation_matrix(1, 0) = -sin(theta);
    z_rotation_matrix(1, 1) = cos(theta);

    Rotate(z_rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(double theta)
{
    RotateZ(theta);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::IsMeshChanging() const
{
    return mMeshChangesDuringSimulation;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMaximumContainingElementsPerProcess() const
{
    unsigned max_num = 0u;
    for (unsigned local_node_index = 0; local_node_index < mNodes.size(); local_node_index++)
    {
        unsigned num = mNodes[local_node_index]->GetNumContainingElements();
        if (num > max_num)
        {
            max_num = num;
        }
    }
    return max_num;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::SetMeshHasChangedSinceLoading()
{
    // We just forget what the original file was, which has the desired effect
    mMeshFileBaseName = "";
}

// Explicit instantiation
template class AbstractMesh<1, 1>;
template class AbstractMesh<1, 2>;
template class AbstractMesh<1, 3>;
template class AbstractMesh<2, 2>;
template class AbstractMesh<2, 3>;
template class AbstractMesh<3, 3>;
