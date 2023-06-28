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

#include "SemMesh.hpp"

template<unsigned DIM>
SemMesh<DIM>::SemMesh(std::vector<Node<DIM>*> nodes,
                      std::vector<SemElement<DIM>*> semElements)
{
    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        Node<DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index = 0; elem_index < semElements.size(); elem_index++)
    {
        SemElement<DIM>* p_temp_sem_element = semElements[elem_index];
        mElements.push_back(p_temp_sem_element);
    }
}

template<unsigned DIM>
SemMesh<DIM>::SemMesh()
{
    // Reset member variables and clear mNodes and mElements
    Clear();
}

template<unsigned DIM>
SemMesh<DIM>::~SemMesh()
{
    // Reset member variables and clear mNodes and mElements
    Clear();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned DIM>
unsigned SemMesh<DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned DIM>
SemElement<DIM>* SemMesh<DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

template<unsigned DIM>
c_vector<double, DIM> SemMesh<DIM>::GetCentroidOfElement(unsigned index)
{
    ///\todo
    c_vector<double, DIM> centroid = zero_vector<double>(DIM);
    return centroid;
}

template<unsigned DIM>
SemMesh<DIM>* SemMesh<DIM>::GetMeshForVtk()
{
    return this;
}

template<unsigned DIM>
void SemMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader)
{
    ///\todo
}

template<unsigned DIM>
void SemMesh<DIM>::Clear()
{
    // Delete elements
    for (unsigned i = 0; i < mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete nodes
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();
}

template<unsigned DIM>
double SemMesh<DIM>::GetVolumeOfElement(unsigned index)
{
    ///\todo
    return 0.0;
}

template<unsigned DIM>
unsigned SemMesh<DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template<unsigned DIM>
unsigned SemMesh<DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template<unsigned DIM>
unsigned SemMesh<DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    return index;
}
