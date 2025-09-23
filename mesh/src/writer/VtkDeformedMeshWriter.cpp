/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "VtkDeformedMeshWriter.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"

template<unsigned DIM>
VtkDeformedMeshWriter<DIM>::VtkDeformedMeshWriter(AbstractTetrahedralMesh<DIM,DIM>* pDeformedOutputMesh, 
    const std::string& rDirectory, 
    const std::string& rBaseName, 
    const bool& rCleanDirectory)
    : VtkMeshWriter<DIM,DIM>(rDirectory,rBaseName,rCleanDirectory)
{
    // if (PetscTools::IsSequential())
    // {
        mpDeformedMesh = new TetrahedralMesh<DIM,DIM>();
   // }
    // else
    // {
    //     mpDeformedMesh = new DistributedTetrahedralMesh<DIM,DIM>();
    //     this->SetParallelFiles(*mpDeformedMesh);
    // }
    
    mpDeformedMesh->ConstructFromMesh(*pDeformedOutputMesh);
}

template<unsigned DIM>
void VtkDeformedMeshWriter<DIM>::SetOutputBaseFileName(const std::string& rOutputBaseFileName)
{
    this->mBaseName = rOutputBaseFileName;
}

template<unsigned DIM>
void VtkDeformedMeshWriter<DIM>::WriteDeformedFiles()
{
    this->WriteFilesUsingMesh(*mpDeformedMesh);
}

template<unsigned DIM>
void VtkDeformedMeshWriter<DIM>::ApplyDeformation(std::vector<c_vector<double,DIM> >& rPositions)
{
    // for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = mpDeformedMesh->GetNodeIteratorBegin();
    //      node_iter != mpDeformedMesh->GetNodeIteratorEnd();
    //      ++node_iter)
    for (unsigned node_index = 0u; node_index<rPositions.size();++node_index)
    {
        ChastePoint<DIM> new_position(rPositions[node_index]);
        try{
            //The next line will change node coordinates. 
            //Note that we do not update element Jacobians (and related quantities)
            //as this mesh is used only in this class for visualization purposes.
            mpDeformedMesh->GetNode(node_index)->SetPoint(new_position);
        }
        // LCOV_EXCL_START
        catch (Exception& e)
        {
            //PetscTools::ReplicateException(true);
            //throw e;

        }
        // LCOV_EXCL_STOP
    }
}
// Explicit instantiation
template class VtkDeformedMeshWriter<2>;
template class VtkDeformedMeshWriter<3>;
