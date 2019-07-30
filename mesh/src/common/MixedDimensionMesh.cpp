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

#include "MixedDimensionMesh.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::MixedDimensionMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod)
    : DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DistributedTetrahedralMesh(partitioningMethod)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::~MixedDimensionMesh()
{
    for (unsigned i=0; i<mCableElements.size(); i++)
    {
        delete mCableElements[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader)
{
    DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ConstructFromMeshReader(rMeshReader);
    //Note that the above method may permute the node (for a parMETIS partition) after construction
    //We have to convert to the permuted form first

    // Add cable elements
    mNumCableElements = rMeshReader.GetNumCableElements();
    //this->mCableElements.reserve(mNumCableElements);

    for (unsigned element_index=0; element_index < mNumCableElements; element_index++)
    {
        ElementData element_data = rMeshReader.GetNextCableElementData();
        //Convert the node indices from the original to the permuted
        if (!this->mNodePermutation.empty())
        {
            for (unsigned j=0; j<2; j++) // cables are always 1d
            {
                element_data.NodeIndices[j] = this->mNodePermutation[ element_data.NodeIndices[j] ];
            }
        }

        //Determine if we own any nodes on this cable element
        bool node_owned = false;
        for (unsigned j=0; j<2; j++) // cables are always 1d
        {
            try
            {
                this->SolveNodeMapping(element_data.NodeIndices[j]);
                node_owned = true;
                break;
            }
            catch (Exception &)
            {
                //We deal with non-owned nodes in the next part
            }
        }

        //If we don't locally own either node, then we don't construct the cable
        if (node_owned)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.reserve(2u);

            for (unsigned j=0; j<2; j++) // cables are always 1d
            {
                //Note that if we own one node on a cable element then we are likely to own the other.
                //If not, we are likely to have a halo.
                //If not, (free-running Purkinje with monodomain mesh?), then this will terminate.
                try
                {
                    nodes.push_back(this->GetNodeOrHaloNode(element_data.NodeIndices[j]) );
                }
                // LCOV_EXCL_START
                catch (Exception&)
                {
                    NEVER_REACHED;
                }
                // LCOV_EXCL_STOP
            }

            Element<1u, SPACE_DIM>* p_element = new Element<1u,SPACE_DIM>(element_index, nodes, false);
            RegisterCableElement(element_index);
            this->mCableElements.push_back(p_element);
            for (unsigned node_index=0; node_index<p_element->GetNumNodes(); ++node_index)
            {
                mNodeToCablesMapping.insert(std::pair<Node<SPACE_DIM>*, Element<1u, SPACE_DIM>*>(
                        p_element->GetNode(node_index), p_element));
            }

            if (rMeshReader.GetNumCableElementAttributes() > 0)
            {
                assert(rMeshReader.GetNumCableElementAttributes() == 1);
                p_element->SetAttribute(element_data.AttributeValue);
            }
        }
    }

    rMeshReader.Reset();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::RegisterCableElement(unsigned index)
{
    mCableElementsMapping[index] = this->mCableElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
   return mNumCableElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalCableElements() const
{
   return mCableElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<1u, SPACE_DIM>* MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElement(unsigned globalElementIndex) const
{
    std::map<unsigned, unsigned>::const_iterator element_position = mCableElementsMapping.find(globalElementIndex);

    if (element_position == mCableElementsMapping.end())
    {
        EXCEPTION("Requested cable element " << globalElementIndex << " does not belong to processor " << PetscTools::GetMyRank());
    }

    unsigned index = element_position->second;

    return mCableElements[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfCableElement( unsigned elementIndex )
{
    //This should not throw in the distributed parallel case
    try
    {
        unsigned tie_break_index = this->GetCableElement(elementIndex)->GetNodeGlobalIndex(0);

        //if it is in my range
        if (this->GetDistributedVectorFactory()->IsGlobalIndexLocal(tie_break_index))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch (Exception &)
    {
        //We don't own this cable element
        return false;
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CableRangeAtNode MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCablesAtNode(const Node<SPACE_DIM>* pNode)
{
    return mNodeToCablesMapping.equal_range(pNode);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CableElementIterator MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElementIteratorBegin() const
{
    return mCableElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CableElementIterator MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElementIteratorEnd() const
{
    return mCableElements.end();
}

// Explicit instantiation
template class MixedDimensionMesh<1,1>;
template class MixedDimensionMesh<1,2>;
template class MixedDimensionMesh<1,3>;
template class MixedDimensionMesh<2,2>;
template class MixedDimensionMesh<2,3>;
template class MixedDimensionMesh<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MixedDimensionMesh)
