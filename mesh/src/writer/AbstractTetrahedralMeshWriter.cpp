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

// Disable PETSc logging of MPI calls (we don't use this, anyway) to fix
// "right-hand operand of comma has no effect" warnings when building with
// PETSc 2.2.1.
#define PETSC_HAVE_BROKEN_RECURSIVE_MACRO

#include "AbstractTetrahedralMeshWriter.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "MixedDimensionMesh.hpp"
#include "Version.hpp"
#include "Exception.hpp"

#include <mpi.h> // For MPI_Send, MPI_Recv

const char* MeshEventHandler::EventName[] = { "Tri write","BinTri write","VTK write","PVTK write", "node write", "ele write", "face write", "ncl write", "comm1","comm2","Total"};

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over elements */
    typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator* pElemIter;
    /** Iterator over boundary elements */
    typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator* pBoundaryElemIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMeshWriter(const std::string& rDirectory,
                   const std::string& rBaseName,
                   const bool clearOutputDir)
    : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
      mpNodeMap(nullptr),
      mNodesPerElement(ELEMENT_DIM+1),
      mNodesPerBoundaryElement(ELEMENT_DIM),
      mpMesh(nullptr),
      mpDistributedMesh(nullptr),
      mpMixedMesh(nullptr),
      mpIters(new MeshWriterIterators<ELEMENT_DIM,SPACE_DIM>),
      mNodeCounterForParallelMesh(0),
      mElementCounterForParallelMesh(0),
      mBoundaryElementCounterForParallelMesh(0),
      mCableElementCounterForParallelMesh(0),
      mFilesAreBinary(false)
{
    mpIters->pNodeIter = nullptr;
    mpIters->pElemIter = nullptr;
    mpIters->pBoundaryElemIter = nullptr;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::~AbstractTetrahedralMeshWriter()
{
    if (mpIters->pNodeIter)
    {
        delete mpIters->pNodeIter;
    }
    if (mpIters->pElemIter)
    {
        delete mpIters->pElemIter;
    }
    if (mpIters->pBoundaryElemIter)
    {
        delete mpIters->pBoundaryElemIter;
    }

    delete mpIters;

    if (mpNodeMap)
    {
        delete mpNodeMap;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    // if we are writing from a mesh..
    assert(PetscTools::AmMaster());
    if (mpMesh)
    {
        std::vector<double> coords(SPACE_DIM);

        //Iterate over the locally-owned nodes
        if ((*(mpIters->pNodeIter)) != mpMesh->GetNodeIteratorEnd())
        {
            // Either this is a sequential mesh (and we own it all)
            // or it's parallel (and the master owns the first chunk)
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = (*(mpIters->pNodeIter))->GetPoint()[j];
            }

            mNodeCounterForParallelMesh=(*(mpIters->pNodeIter))->GetIndex() + 1;//Ready for when we run out of local nodes

            ++(*(mpIters->pNodeIter));
            return coords;
        }

        // If we didn't return then the iterator has reached the end of the local nodes.
        // It must be a parallel mesh and we are expecting messages...

        assert( mpDistributedMesh != nullptr );

        MPI_Status status;
        status.MPI_ERROR = MPI_SUCCESS; //For MPICH2
        // do receive, convert to std::vector on master
        boost::scoped_array<double> raw_coords(new double[SPACE_DIM]);
        MPI_Recv(raw_coords.get(), SPACE_DIM, MPI_DOUBLE, MPI_ANY_SOURCE, mNodeCounterForParallelMesh, PETSC_COMM_WORLD, &status);
        assert(status.MPI_ERROR == MPI_SUCCESS);
        for (unsigned j=0; j<coords.size(); j++)
        {
            coords[j] = raw_coords[j];
        }

        mNodeCounterForParallelMesh++;
        return coords;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextNode();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    assert(PetscTools::AmMaster());
    // if we are writing from a mesh..
    if (mpMesh)
    {
        ElementData elem_data;
        elem_data.NodeIndices.resize(mNodesPerElement);

        if (mpDistributedMesh == nullptr) // not using parallel mesh
        {
            // Use the iterator
            assert(this->mNumElements == mpMesh->GetNumElements());

            for (unsigned j=0; j<elem_data.NodeIndices.size(); j++)
            {
                unsigned old_index = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
                elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
            }
            // Set attribute
            elem_data.AttributeValue = (*(mpIters->pElemIter))->GetAttribute();

            ++(*(mpIters->pElemIter));

            return elem_data;
        }
        else // Parallel mesh
        {
            // Use the mElementCounterForParallelMesh variable to identify next element
            if (mpDistributedMesh->CalculateDesignatedOwnershipOfElement(mElementCounterForParallelMesh) == true)
            {
                // Master owns this element
                Element<ELEMENT_DIM, SPACE_DIM>* p_element = mpDistributedMesh->GetElement(mElementCounterForParallelMesh);
                assert(elem_data.NodeIndices.size() == mNodesPerElement);
                assert(!p_element->IsDeleted());
                // Master can use the local data to recover node indices & attribute
                for (unsigned j=0; j<mNodesPerElement; j++)
                {
                    elem_data.NodeIndices[j] = p_element->GetNodeGlobalIndex(j);
                }
                elem_data.AttributeValue = p_element->GetAttribute();
            }
            else
            {
                //Master doesn't own this element.
                UnpackElement(elem_data, mElementCounterForParallelMesh, mNodesPerElement, this->mNumNodes +  mElementCounterForParallelMesh);
            }
            // increment element counter
            mElementCounterForParallelMesh++;

            return elem_data; // non-master processors will return garbage here - but they should never write to file
        }
    }
    else // not writing from a mesh
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextElement();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextBoundaryElement()
{
    assert(PetscTools::AmMaster());
    // if we are writing from a mesh..
    if (mpMesh)
    {
        ElementData boundary_elem_data;
        boundary_elem_data.NodeIndices.resize(mNodesPerBoundaryElement);

        if (mpDistributedMesh == nullptr) // not using parallel mesh
        {
            // Use the iterator
            assert(this->mNumBoundaryElements==mpMesh->GetNumBoundaryElements());

            for (unsigned j=0; j<boundary_elem_data.NodeIndices.size(); j++)
            {
                unsigned old_index = (*(*(mpIters->pBoundaryElemIter)))->GetNodeGlobalIndex(j);
                boundary_elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
            }
            boundary_elem_data.AttributeValue = (*(*(mpIters->pBoundaryElemIter)))->GetAttribute();

            ++(*(mpIters->pBoundaryElemIter));
            return boundary_elem_data;
        }
        else // Parallel mesh
        {
            // Use the mElementCounterForParallelMesh variable to identify next element
            if (mpDistributedMesh->CalculateDesignatedOwnershipOfBoundaryElement(mBoundaryElementCounterForParallelMesh) == true)
            {
                // Master owns this boundary element
                BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = mpDistributedMesh->GetBoundaryElement(mBoundaryElementCounterForParallelMesh);
                assert(boundary_elem_data.NodeIndices.size() == ELEMENT_DIM);
                assert(!p_boundary_element->IsDeleted());

                // Master can use the local data to recover node indices & attribute
                for (unsigned j=0; j<ELEMENT_DIM; j++)
                {
                    boundary_elem_data.NodeIndices[j] = p_boundary_element->GetNodeGlobalIndex(j);
                }
                boundary_elem_data.AttributeValue = p_boundary_element->GetAttribute();
            }
            else
            {
                //Master doesn't own this boundary element.
                UnpackElement(boundary_elem_data, mBoundaryElementCounterForParallelMesh, ELEMENT_DIM, this->mNumNodes + this->mNumElements + mBoundaryElementCounterForParallelMesh);
            }
            // increment element counter
            mBoundaryElementCounterForParallelMesh++;

            return boundary_elem_data; // non-master processors will return garbage here - but they should never write to file
        }
    }
    else // not writing from a mesh
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextBoundaryElement();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextCableElement()
{
    assert(PetscTools::AmMaster());

    // if we are writing from a mesh..
    if (mpMesh)
    {
        // Need to be using a MixedDimensionMesh or there will be no cable data
        assert(mpMixedMesh);

        ElementData elem_data;
        elem_data.NodeIndices.resize(2);

        // Use the mCableElementCounterForParallelMesh variable to identify next element
        if (mpMixedMesh->CalculateDesignatedOwnershipOfCableElement(mCableElementCounterForParallelMesh) == true)
        {
            // Master owns this element
            Element<1, SPACE_DIM>* p_element = mpMixedMesh->GetCableElement(mCableElementCounterForParallelMesh);
            assert(!p_element->IsDeleted());

            // Master can use the local data to recover node indices & attribute
            for (unsigned j=0; j<2; j++)
            {
                elem_data.NodeIndices[j] = p_element->GetNodeGlobalIndex(j);
            }
            elem_data.AttributeValue = p_element->GetAttribute();
        }
        else
        {
            // Master doesn't own this element
            UnpackElement(elem_data, mCableElementCounterForParallelMesh, 2, this->mNumNodes + this->mNumElements + this->mNumBoundaryElements + mCableElementCounterForParallelMesh);
        }
        // Increment element counter
        mCableElementCounterForParallelMesh++;

        return elem_data; // non-master processors will return garbage here - but they should never write to file
    }
    else // not writing from a mesh
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextCableElement();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteNclFile(
        AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        bool invertMeshPermutation)
{
    MeshEventHandler::BeginEvent(MeshEventHandler::NCL);
    unsigned max_elements_all;
    if (PetscTools::IsSequential())
    {
        max_elements_all = rMesh.CalculateMaximumContainingElementsPerProcess();
    }
    else
    {
        unsigned max_elements_per_process = rMesh.CalculateMaximumContainingElementsPerProcess();
        MPI_Allreduce(&max_elements_per_process, &max_elements_all, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD);
    }

    std::string node_connect_list_file_name = this->mBaseName + ".ncl";
    if (invertMeshPermutation && !rMesh.rGetNodePermutation().empty())
    {
        node_connect_list_file_name += "-tmp";
    }

    PetscTools::BeginRoundRobin();
    {
        out_stream p_ncl_file = out_stream(nullptr);

        if (PetscTools::AmMaster())
        {
            //Open the file for the first time
            p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(node_connect_list_file_name, std::ios::binary | std::ios::trunc);

            // Write the ncl header
            *p_ncl_file << rMesh.GetNumNodes() << "\t";
            *p_ncl_file << max_elements_all << "\t";
            *p_ncl_file << "\tBIN\n";
        }
        else
        {
            // Append to the existing file
            p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(node_connect_list_file_name, std::ios::binary | std::ios::app);
        }

        // Write each node's data
        unsigned default_marker = UINT_MAX;

        typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
        for (NodeIterType iter = rMesh.GetNodeIteratorBegin();
             iter != rMesh.GetNodeIteratorEnd();
             ++iter)
        {
            // Get the containing element indices from the node's set and sort them
            std::set<unsigned>& r_elem_set = iter->rGetContainingElementIndices();
            std::vector<unsigned> elem_vector(r_elem_set.begin(), r_elem_set.end());
            std::sort(elem_vector.begin(), elem_vector.end());
            // Pad the vector with unsigned markers
            for (unsigned elem_index=elem_vector.size();  elem_index<max_elements_all; elem_index++)
            {
                elem_vector.push_back(default_marker);
            }
            assert(elem_vector.size() == max_elements_all);
            // Write raw data out of std::vector into the file
            if (max_elements_all > 0u)
            {
                p_ncl_file->write((char*)&elem_vector[0], elem_vector.size()*sizeof(unsigned));
            }
        }

        if (PetscTools::AmTopMost())
        {
            *p_ncl_file << "#\n# " + ChasteBuildInfo::GetProvenanceString();
        }

        p_ncl_file->close();
    }
    PetscTools::EndRoundRobin();

    if (invertMeshPermutation && PetscTools::AmMaster() && !rMesh.rGetNodePermutation().empty() && max_elements_all > 0u)
    {
        // Open files
        const std::string real_node_connect_list_file_name = this->mBaseName + ".ncl";
        out_stream p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(real_node_connect_list_file_name, std::ios::binary | std::ios::trunc);
        FileFinder temp_ncl_path = this->mpOutputFileHandler->FindFile(node_connect_list_file_name);
        std::ifstream temp_ncl_file(temp_ncl_path.GetAbsolutePath().c_str(), std::ios::binary);
        // Copy the header
        std::string header_line;
        getline(temp_ncl_file, header_line, '\n');
        (*p_ncl_file) << header_line << "\n";
        const std::streampos data_start = temp_ncl_file.tellg();
        const std::streamoff item_width = max_elements_all * sizeof(unsigned);
        // Copy the binary data, permuted
        std::vector<unsigned> elem_vector(max_elements_all);
        for (unsigned node_index=0; node_index<rMesh.GetNumAllNodes(); node_index++)
        {
            unsigned permuted_index = rMesh.rGetNodePermutation()[node_index];
            temp_ncl_file.seekg(data_start + item_width * permuted_index, std::ios_base::beg);
            temp_ncl_file.read((char*)&elem_vector[0], max_elements_all*sizeof(unsigned));
            p_ncl_file->write((char*)&elem_vector[0], max_elements_all*sizeof(unsigned));
        }
        // Footer
        *p_ncl_file << "#\n# " + ChasteBuildInfo::GetProvenanceString();
        p_ncl_file->close();
        // Remove temp file
        remove(temp_ncl_path.GetAbsolutePath().c_str());
    }
    PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteNclFile");
    MeshEventHandler::EndEvent(MeshEventHandler::NCL);
}

///\todo #1322 Mesh should be const
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
      AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
      bool keepOriginalElementIndexing)
{
    this->mpMeshReader = nullptr;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumElements = mpMesh->GetNumElements();
    this->mNumBoundaryElements = mpMesh->GetNumBoundaryElements();
    this->mNumCableElements = mpMesh->GetNumCableElements();

    typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator ElemIterType;
    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator BoundaryElemIterType;
    mpIters->pBoundaryElemIter = new BoundaryElemIterType(mpMesh->GetBoundaryElementIteratorBegin());

    // Use this process's first element to gauge the size of all the elements
    if ((*(mpIters->pElemIter)) != mpMesh->GetElementIteratorEnd())
    {
        mNodesPerElement = (*(mpIters->pElemIter))->GetNumNodes();
    }

    // Use this process's first boundary element to gauge the size of all the boundary elements
    if ((*(mpIters->pBoundaryElemIter)) != mpMesh->GetBoundaryElementIteratorEnd())
    {
        mNodesPerBoundaryElement = (*(*(mpIters->pBoundaryElemIter)))->GetNumNodes();
    }

    //Connectivity file is written when we write to a binary file (only available for TrianglesMeshWriter) and if we are preserving the element order
    if (this->mFilesAreBinary && keepOriginalElementIndexing)
    {
        WriteNclFile(*mpMesh);
    }

    // Have we got a parallel mesh?
    ///\todo #1322 This should be const too
    mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);

    // Have we got a MixedDimensionMesh?
    ///\todo #1322,  This should be const too
    mpMixedMesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* >(this->mpMesh);

    if (mpDistributedMesh != nullptr)
    {
        // It's a parallel mesh
        WriteFilesUsingParallelMesh(keepOriginalElementIndexing);
        return;
    }

    if (!PetscTools::AmMaster())
    {
        PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingMesh"); //Paired with Master process writing files
        return;
    }

    // Set up node map if we might have deleted nodes
    unsigned node_map_index = 0;
    if (mpMesh->IsMeshChanging())
    {
        mpNodeMap = new NodeMap(1 + mpMesh->GetMaximumNodeIndex());
        for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
        {
            mpNodeMap->SetNewIndex(it->GetIndex(), node_map_index++);
        }
    }

    this->WriteFiles();
    PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingMesh"); // Paired with waiting Slave processes
    delete mpIters->pNodeIter;
    mpIters->pNodeIter = nullptr;
    delete mpIters->pElemIter;
    mpIters->pElemIter = nullptr;
    delete mpIters->pBoundaryElemIter;
    mpIters->pBoundaryElemIter = nullptr;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReaderAndMesh(
        AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
        AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    WriteNclFile(rMesh, true);
    this->WriteFilesUsingMeshReader(rMeshReader);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingParallelMesh(bool keepOriginalElementIndexing)
{
    if (keepOriginalElementIndexing)
    {
        // Master goes on to write as usual
        if (PetscTools::AmMaster())
        {
            this->WriteFiles();
        }
        else
        {
//            PetscTools::Barrier("DodgyBarrierBeforeNODE");
            MeshEventHandler::BeginEvent(MeshEventHandler::NODE);
            boost::scoped_array<double> raw_coords(new double[SPACE_DIM]);
            // Slaves concentrate the Nodes
            typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
            for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
            {
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    raw_coords[j] = it->GetPoint()[j];
                }
                MPI_Ssend(raw_coords.get(), SPACE_DIM, MPI_DOUBLE, 0, it->GetIndex(), PETSC_COMM_WORLD);//Nodes sent with positive tags
            }
//            PetscTools::Barrier("DodgyBarrierAfterNODE");
            MeshEventHandler::EndEvent(MeshEventHandler::NODE);

            MeshEventHandler::BeginEvent(MeshEventHandler::ELE);
            // Slaves concentrate the Elements for which they are owners
            boost::scoped_array<unsigned> raw_indices(new unsigned[mNodesPerElement]); // Assuming that we don't have parallel quadratic elements
            typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator ElementIterType;
            for (ElementIterType it = mpMesh->GetElementIteratorBegin(); it != mpMesh->GetElementIteratorEnd(); ++it)
            {
                unsigned index = it->GetIndex();
                if (mpDistributedMesh->CalculateDesignatedOwnershipOfElement(index) == true)
                {
                    for (unsigned j=0; j<mNodesPerElement; j++)
                    {
                        raw_indices[j] = it->GetNodeGlobalIndex(j);
                    }
                    PostElement(index, raw_indices.get(), mNodesPerElement, this->mNumNodes + index, it->GetAttribute());
                }
            }
//            PetscTools::Barrier("DodgyBarrierAfterELE");
            MeshEventHandler::EndEvent(MeshEventHandler::ELE);
            MeshEventHandler::BeginEvent(MeshEventHandler::FACE);
            // Slaves concentrate the Faces for which they are owners (not in 1-d)
            if (ELEMENT_DIM != 1)  /// \todo #2351 Also exclude VTK writer
            {
                boost::scoped_array<unsigned> raw_face_indices(new unsigned[ELEMENT_DIM]); // Assuming that we don't have parallel quadratic meshes
                typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator BoundaryElementIterType;
                for (BoundaryElementIterType it = mpMesh->GetBoundaryElementIteratorBegin(); it != mpMesh->GetBoundaryElementIteratorEnd(); ++it)
                {
                    unsigned index = (*it)->GetIndex();
                    if (mpDistributedMesh->CalculateDesignatedOwnershipOfBoundaryElement(index) == true)
                    {
                        for (unsigned j=0; j<ELEMENT_DIM; j++)
                        {
                            raw_face_indices[j] = (*it)->GetNodeGlobalIndex(j);
                        }
                        PostElement(index, raw_face_indices.get(), ELEMENT_DIM, this->mNumNodes + this->mNumElements + index, (*it)->GetAttribute());
                    }
                }
            }
//            PetscTools::Barrier("DodgyBarrierAfterFACE");
            MeshEventHandler::EndEvent(MeshEventHandler::FACE);

            // Slaves concentrate the cable elements for which they are owners
            if (mpMixedMesh)
            {
                typedef typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator CableElementIterType;
                for (CableElementIterType it = mpMixedMesh->GetCableElementIteratorBegin(); it != mpMixedMesh->GetCableElementIteratorEnd(); ++it)
                {
                    unsigned index =(*it)->GetIndex();
                    if (mpMixedMesh->CalculateDesignatedOwnershipOfCableElement(index) == true)
                    {
                        unsigned raw_cable_element_indices[2];
                        for (unsigned j=0; j<2; j++)
                        {
                            raw_cable_element_indices[j] = (*it)->GetNodeGlobalIndex(j);
                        }
                        PostElement(index, raw_cable_element_indices, 2, this->mNumNodes + this->mNumElements + this->mNumBoundaryElements + index, (*it)->GetAttribute());
                    }
                }
            }
        }
        PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingParallelMesh");
    }
    else
    {
        PetscTools::BeginRoundRobin();

        if (PetscTools::AmMaster())
        {
            // Make sure headers are written first
            assert(PetscTools::GetMyRank() == 0);
            CreateFilesWithHeaders();
        }

        AppendLocalDataToFiles();

        if (PetscTools::AmTopMost())
        {
            // Make sure footers are written last
            assert(PetscTools::GetMyRank() == PetscTools::GetNumProcs()-1);
            WriteFilesFooter();
        }

        PetscTools::EndRoundRobin();
    }
}

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::CreateFilesWithHeaders()
{
    // If control reaches this line you haven't implemented the optimised
    // parallel write for whichever visualiser you are writing for.
    NEVER_REACHED;
}
// LCOV_EXCL_STOP


// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::AppendLocalDataToFiles()
{
    // If control reaches this line you haven't implemented the optimised
    // parallel write for whichever visualiser you are writing for.
    NEVER_REACHED;
}
// LCOV_EXCL_STOP


// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesFooter()
{
    // If control reaches this line you haven't implemented the optimised
    // parallel write for whichever visualiser you are writing for.
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

// Explicit instantiation
template class AbstractTetrahedralMeshWriter<1,1>;
template class AbstractTetrahedralMeshWriter<1,2>;
template class AbstractTetrahedralMeshWriter<1,3>;
template class AbstractTetrahedralMeshWriter<2,2>;
template class AbstractTetrahedralMeshWriter<2,3>;
template class AbstractTetrahedralMeshWriter<3,3>;
