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

#include <boost/scoped_array.hpp>
#include "VtkMeshWriter.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "MixedDimensionMesh.hpp"
#include "NodesOnlyMesh.hpp"

#ifdef CHASTE_VTK
#include "vtkQuadraticTetra.h"
#include "vtkQuadraticTriangle.h"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshWriter<ELEMENT_DIM, SPACE_DIM>::VtkMeshWriter(const std::string& rDirectory,
                     const std::string& rBaseName,
                     const bool& rCleanDirectory)
    : AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, rCleanDirectory),
      mWriteParallelFiles(false)
{
    this->mIndexFromZero = true;

    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::~VtkMeshWriter()
{
    mpVtkUnstructedMesh->Delete(); // Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::MakeVtkMesh()
{
    //Construct nodes aka as Points
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->GetNextNode(); //this->mNodeData[item_num];
        // Add zeroes if the dimension is below 3
        for (unsigned dim=SPACE_DIM; dim<3; dim++)
        {
            current_item.push_back(0.0);//For y and z-coordinates if necessary
        }
        assert(current_item.size() == 3);
        p_pts->InsertPoint(item_num, current_item[0], current_item[1], current_item[2]);
    }
    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); //Reference counted

    //Construct elements aka Cells
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        std::vector<unsigned> current_element = this->GetNextElement().NodeIndices; // this->mElementData[item_num];

        assert((current_element.size() == ELEMENT_DIM + 1) || (current_element.size() == (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2));

        vtkCell* p_cell=nullptr;
        if (ELEMENT_DIM == 3 && current_element.size() == 4)
        {
            p_cell = vtkTetra::New();
        }
        else if (ELEMENT_DIM == 3 && current_element.size() == 10)
        {
            p_cell = vtkQuadraticTetra::New();
        }
        else if (ELEMENT_DIM == 2 && current_element.size() == 3)
        {
            p_cell = vtkTriangle::New();
        }
        else if (ELEMENT_DIM == 2 && current_element.size() == 6)
        {
            p_cell = vtkQuadraticTriangle::New();
        }
        else if (ELEMENT_DIM == 1)
        {
            p_cell = vtkLine::New();
        }

        //Set the linear nodes
        vtkIdList* p_cell_id_list = p_cell->GetPointIds();
        for (unsigned j = 0; j < current_element.size(); ++j)
        {
            p_cell_id_list->SetId(j, current_element[j]);
        }

        //VTK defines the node ordering in quadratic triangles differently to Chaste, so they must be treated as a special case
        if (SPACE_DIM == 2 && current_element.size() == 6)
        {
            p_cell_id_list->SetId(3, current_element[5]);
            p_cell_id_list->SetId(4, current_element[3]);
            p_cell_id_list->SetId(5, current_element[4]);
        }

        mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
        p_cell->Delete(); //Reference counted
    }

    if (SPACE_DIM > 1)
    {
        /// \todo #2351 Temporary workaround for parallel writer
        for (unsigned item_num=0; item_num<this->GetNumBoundaryFaces(); item_num++)
        {
            this->GetNextBoundaryElement();
        }
    }

    //If necessary, construct cables
    if (this->GetNumCableElements() > 0)
    {
        AugmentCellData();
        //Make a blank cell radius data for the regular elements
        std::vector<double> radii(this->GetNumElements(), 0.0);
        for (unsigned item_num=0; item_num<this->GetNumCableElements(); item_num++)
        {
            ElementData cable_element_data = this->GetNextCableElement();
            std::vector<unsigned> current_element = cable_element_data.NodeIndices;
            radii.push_back(cable_element_data.AttributeValue);
            assert(current_element.size() == 2);
            vtkCell* p_cell=vtkLine::New();
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            for (unsigned j = 0; j < 2; ++j)
            {
                p_cell_id_list->SetId(j, current_element[j]);
            }
            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); //Reference counted
        }
        AddCellData("Cable radius", radii);

    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddProvenance(std::string fileName)
{
    std::string comment = "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->";

    out_stream p_vtu_file = this->mpOutputFileHandler->OpenOutputFile(fileName, std::ios::out | std::ios::app);

    *p_vtu_file << "\n" << comment << "\n";
    p_vtu_file->close();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteFiles()
{
    // Using separate scope here to make sure file is properly closed before re-opening it to add provenance info.
    {
        MakeVtkMesh();
        assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
        vtkXMLUnstructuredGridWriter* p_writer = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION >= 6
        p_writer->SetInputData(mpVtkUnstructedMesh);
#else
        p_writer->SetInput(mpVtkUnstructedMesh);
#endif
        std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName+".vtu";
        p_writer->SetFileName(vtk_file_name.c_str());
        //p_writer->PrintSelf(std::cout, vtkIndent());
        p_writer->Write();
        p_writer->Delete(); //Reference counted
    }

    AddProvenance(this->mBaseName + ".vtu");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AugmentCellData()
{
    unsigned num_cell_arrays = mpVtkUnstructedMesh->GetCellData()->GetNumberOfArrays();
    for (unsigned i = 0; i < num_cell_arrays; i++)
    {
        vtkDataArray* array = mpVtkUnstructedMesh->GetCellData()->GetArray(i);

        //Check data was the correct size before the cables were added
        unsigned num_cable_pads = this->GetNumCableElements();
        if (mWriteParallelFiles)
        {
            assert((unsigned)array->GetNumberOfTuples() == this->mpDistributedMesh->GetNumLocalElements());
            num_cable_pads =  this->mpMixedMesh->GetNumLocalCableElements();
        }
        else
        {
            assert((unsigned)array->GetNumberOfTuples() == this->GetNumElements());
        }

        //Check that tuples of size 3 will be big enough for padding the rest of the data
        assert(array->GetNumberOfComponents() <= 3);
        double null_data[3] = {0.0, 0.0, 0.0};

        //Pad data
        for (unsigned new_index = 0; new_index <  num_cable_pads; new_index++)
        {
            array->InsertNextTuple(null_data);
        }
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddCellData(std::string dataName, std::vector<c_vector<double, SPACE_DIM> > dataPayload)
{
    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(3);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            p_vectors->InsertNextValue(dataPayload[i][j]);
        }
        //When SPACE_DIM<3, then pad
        for (unsigned j=SPACE_DIM; j<3; j++)
        {
            p_vectors->InsertNextValue(0.0);
        }
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddTensorCellData(std::string dataName, std::vector<c_vector<double,SPACE_DIM*(SPACE_DIM+1)/2> > dataPayload)
{
    assert(SPACE_DIM != 1);    // LCOV_EXCL_LINE

    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(SPACE_DIM*SPACE_DIM);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        if (SPACE_DIM == 2)
        {
            p_vectors->InsertNextValue(dataPayload[i](0)); //a11
            p_vectors->InsertNextValue(dataPayload[i](1)); //a12
            p_vectors->InsertNextValue(dataPayload[i](1)); //a21
            p_vectors->InsertNextValue(dataPayload[i](2)); //a22
        }
        else if (SPACE_DIM == 3)
        {
            p_vectors->InsertNextValue(dataPayload[i](0)); //a11
            p_vectors->InsertNextValue(dataPayload[i](1)); //a12
            p_vectors->InsertNextValue(dataPayload[i](2)); //a13
            p_vectors->InsertNextValue(dataPayload[i](1)); //a21
            p_vectors->InsertNextValue(dataPayload[i](3)); //a22
            p_vectors->InsertNextValue(dataPayload[i](4)); //a23
            p_vectors->InsertNextValue(dataPayload[i](2)); //a31
            p_vectors->InsertNextValue(dataPayload[i](4)); //a32
            p_vectors->InsertNextValue(dataPayload[i](5)); //a33
        }
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddTensorCellData(std::string dataName, std::vector<c_matrix<double,SPACE_DIM,SPACE_DIM> > dataPayload)
{
    assert(SPACE_DIM != 1);    // LCOV_EXCL_LINE

    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(SPACE_DIM*SPACE_DIM);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        if (SPACE_DIM == 2)
        {
            p_vectors->InsertNextValue(dataPayload[i](0,0)); //a11
            p_vectors->InsertNextValue(dataPayload[i](0,1)); //a12
            p_vectors->InsertNextValue(dataPayload[i](1,0)); //a21
            p_vectors->InsertNextValue(dataPayload[i](1,1)); //a22
        }
        else if (SPACE_DIM == 3)
        {
            p_vectors->InsertNextValue(dataPayload[i](0,0)); //a11
            p_vectors->InsertNextValue(dataPayload[i](0,1)); //a12
            p_vectors->InsertNextValue(dataPayload[i](0,2)); //a13
            p_vectors->InsertNextValue(dataPayload[i](1,0)); //a21
            p_vectors->InsertNextValue(dataPayload[i](1,1)); //a22
            p_vectors->InsertNextValue(dataPayload[i](1,2)); //a23
            p_vectors->InsertNextValue(dataPayload[i](2,0)); //a31
            p_vectors->InsertNextValue(dataPayload[i](2,1)); //a32
            p_vectors->InsertNextValue(dataPayload[i](2,2)); //a33
        }
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());

    if (mWriteParallelFiles && this->mpDistributedMesh != nullptr)
    {
        // In parallel, the vector we pass will only contain the values from the privately owned nodes.
        // To get the values from the halo nodes (which will be inserted at the end of the vector we need to
        // communicate with the equivalent vectors on other processes.

        // resize the payload data to include halos
        assert( dataPayload.size() == this->mpDistributedMesh->GetNumLocalNodes() );
        dataPayload.resize( this->mpDistributedMesh->GetNumLocalNodes() + this->mpDistributedMesh->GetNumHaloNodes() );


        // then do the communication
        for ( unsigned rank_offset = 1; rank_offset < PetscTools::GetNumProcs(); rank_offset++ )
        {
            unsigned send_to      = (PetscTools::GetMyRank() + rank_offset) % (PetscTools::GetNumProcs());
            unsigned receive_from = (PetscTools::GetMyRank() + PetscTools::GetNumProcs()- rank_offset ) % (PetscTools::GetNumProcs());

            unsigned number_of_nodes_to_send    = mNodesToSendPerProcess[send_to].size();
            unsigned number_of_nodes_to_receive = mNodesToReceivePerProcess[receive_from].size();

            boost::scoped_array<double> send_data(new double[number_of_nodes_to_send]);
            boost::scoped_array<double> receive_data(new double[number_of_nodes_to_receive]);
            // Pack
            for (unsigned node = 0; node < number_of_nodes_to_send; node++)
            {
                unsigned global_node_index = mNodesToSendPerProcess[send_to][node];
                unsigned local_node_index = global_node_index
                            - this->mpDistributedMesh->GetDistributedVectorFactory()->GetLow();
                send_data[node] = dataPayload[local_node_index];
            }
            {
                // Send
                int ret;
                MPI_Status status;
                ret = MPI_Sendrecv(send_data.get(), number_of_nodes_to_send,
                                   MPI_DOUBLE,
                                   send_to, 0,
                                   receive_data.get(),  number_of_nodes_to_receive,
                                   MPI_DOUBLE,
                                   receive_from, 0,
                                   PETSC_COMM_WORLD, &status);
                UNUSED_OPT(ret);
                assert ( ret == MPI_SUCCESS );
            }

            // Unpack
            for ( unsigned node = 0; node < number_of_nodes_to_receive; node++ )
            {
                unsigned global_node_index = mNodesToReceivePerProcess[receive_from][node];
                unsigned halo_index = mGlobalToNodeIndexMap[global_node_index];
                assert( halo_index >= this->mpDistributedMesh->GetNumLocalNodes() );
                dataPayload[halo_index] = receive_data[node];
            }

        }
    }

    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); //Reference counted
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddPointData(std::string dataName, std::vector<c_vector<double, SPACE_DIM> > dataPayload)
{
    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());

    if (mWriteParallelFiles)
    {
        // In parallel, the vector we pass will only contain the values from the privately owned nodes.
        // To get the values from the halo nodes (which will be inserted at the end of the vector we need to
        // communicate with the equivalent vectors on other processes.

        // resize the payload data to include halos
        assert( dataPayload.size() == this->mpDistributedMesh->GetNumLocalNodes() );
        dataPayload.resize( this->mpDistributedMesh->GetNumLocalNodes() + this->mpDistributedMesh->GetNumHaloNodes() );

        // then do the communication
        for ( unsigned rank_offset = 1; rank_offset < PetscTools::GetNumProcs(); rank_offset++ )
        {
            unsigned send_to      = (PetscTools::GetMyRank() + rank_offset) % (PetscTools::GetNumProcs());
            unsigned receive_from = (PetscTools::GetMyRank() + PetscTools::GetNumProcs()- rank_offset ) % (PetscTools::GetNumProcs());

            unsigned number_of_nodes_to_send    = mNodesToSendPerProcess[send_to].size();
            unsigned number_of_nodes_to_receive = mNodesToReceivePerProcess[receive_from].size();

            boost::scoped_array<double> send_data(new double[number_of_nodes_to_send * SPACE_DIM]);
            boost::scoped_array<double> receive_data(new double[number_of_nodes_to_receive * SPACE_DIM]);

            for (unsigned node = 0; node < number_of_nodes_to_send; node++)
            {
                unsigned global_node_index = mNodesToSendPerProcess[send_to][node];
                unsigned local_node_index = global_node_index
                            - this->mpDistributedMesh->GetDistributedVectorFactory()->GetLow();
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    send_data[ node*SPACE_DIM + j ] = dataPayload[local_node_index][j];
                }
            }

                int ret;
                MPI_Status status;
                ret = MPI_Sendrecv(send_data.get(), number_of_nodes_to_send * SPACE_DIM,
                                   MPI_DOUBLE,
                                   send_to, 0,
                                   receive_data.get(),  number_of_nodes_to_receive * SPACE_DIM,
                                   MPI_DOUBLE,
                                   receive_from, 0,
                                   PETSC_COMM_WORLD, &status);
                UNUSED_OPT(ret);
                assert ( ret == MPI_SUCCESS );

            // Unpack
            for ( unsigned node = 0; node < number_of_nodes_to_receive; node++ )
            {
                unsigned global_node_index = mNodesToReceivePerProcess[receive_from][node];
                unsigned halo_index = mGlobalToNodeIndexMap[global_node_index];
                assert( halo_index >= this->mpDistributedMesh->GetNumLocalNodes() );
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    dataPayload[halo_index][j] = receive_data[ node*SPACE_DIM + j ];
                }
            }
        }
    }

    p_vectors->SetNumberOfComponents(3);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            p_vectors->InsertNextValue(dataPayload[i][j]);
        }
        //When SPACE_DIM<3, then pad
        for (unsigned j=SPACE_DIM; j<3; j++)
        {
            p_vectors->InsertNextValue(0.0);
        }
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddTensorPointData(std::string dataName, std::vector<c_matrix<double,SPACE_DIM,SPACE_DIM> > dataPayload)
{
    assert(SPACE_DIM != 1);    // LCOV_EXCL_LINE

    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(SPACE_DIM*SPACE_DIM);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        if (SPACE_DIM == 2)
        {
            p_vectors->InsertNextValue(dataPayload[i](0,0)); //a11
            p_vectors->InsertNextValue(dataPayload[i](0,1)); //a12
            p_vectors->InsertNextValue(dataPayload[i](1,0)); //a21
            p_vectors->InsertNextValue(dataPayload[i](1,1)); //a22
        }
        else if (SPACE_DIM == 3)
        {
            p_vectors->InsertNextValue(dataPayload[i](0,0)); //a11
            p_vectors->InsertNextValue(dataPayload[i](0,1)); //a12
            p_vectors->InsertNextValue(dataPayload[i](0,2)); //a13
            p_vectors->InsertNextValue(dataPayload[i](1,0)); //a21
            p_vectors->InsertNextValue(dataPayload[i](1,1)); //a22
            p_vectors->InsertNextValue(dataPayload[i](1,2)); //a23
            p_vectors->InsertNextValue(dataPayload[i](2,0)); //a31
            p_vectors->InsertNextValue(dataPayload[i](2,1)); //a32
            p_vectors->InsertNextValue(dataPayload[i](2,2)); //a33
        }
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::SetParallelFiles( AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh )
{
    //Have we got a distributed mesh?
    this->mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);
    mpNodesOnlyMesh = dynamic_cast<NodesOnlyMesh<SPACE_DIM>* >(&rMesh);

    if (this->mpDistributedMesh == nullptr && mpNodesOnlyMesh == nullptr)
    {
        EXCEPTION("Cannot write parallel files using a sequential mesh");
    }

    if (PetscTools::IsSequential())
    {
        return;     // mWriteParallelFiles is not set sequentially (so we don't set up data exchange machinery)
    }

    mWriteParallelFiles = true;

    // Populate the global to node index map (as this will be required to add point data)

    //Node index that we are writing to VTK (index into mNodes and mHaloNodes as if they were concatenated)
    unsigned index = 0;

    // Owned nodes
    for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
         node_iter != rMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        mGlobalToNodeIndexMap[node_iter->GetIndex()] = index;
        index++;
    }

    // Halo nodes
    if (this->mpDistributedMesh)
    {
        for (typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::HaloNodeIterator halo_iter=this->mpDistributedMesh->GetHaloNodeIteratorBegin();
                halo_iter != this->mpDistributedMesh->GetHaloNodeIteratorEnd();
                ++halo_iter)
        {
            mGlobalToNodeIndexMap[(*halo_iter)->GetIndex()] = index;
            index++;
        }

        //Calculate the halo exchange so that node-wise payloads can be communicated
        this->mpDistributedMesh->CalculateNodeExchange( mNodesToSendPerProcess, mNodesToReceivePerProcess );
    }
}

///\todo #1322 Mesh should be const
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
      AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
      bool keepOriginalElementIndexing)
{
    // Have we got a parallel mesh?
    this->mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);
    this->mpMixedMesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);

    if (PetscTools::IsSequential() || !mWriteParallelFiles || (this->mpDistributedMesh == nullptr && mpNodesOnlyMesh == nullptr))
    {
        AbstractTetrahedralMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteFilesUsingMesh( rMesh,keepOriginalElementIndexing );
    }
    else
    {
        //Make the local mesh into a VtkMesh
        vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
        p_pts->GetData()->SetName("Vertex positions");

        // Owned nodes
        for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
             node_iter != rMesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double, SPACE_DIM> current_item = node_iter->rGetLocation();
            if (SPACE_DIM == 3)
            {
                p_pts->InsertNextPoint(current_item[0], current_item[1], current_item[2]);
            }
            else if (SPACE_DIM == 2)
            {
                p_pts->InsertNextPoint(current_item[0], current_item[1], 0.0);
            }
            else // (SPACE_DIM == 1)
            {
                p_pts->InsertNextPoint(current_item[0], 0.0, 0.0);
            }
        }

        // Halo nodes
        if (this->mpDistributedMesh)
        {
            for (typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::HaloNodeIterator halo_iter=this->mpDistributedMesh->GetHaloNodeIteratorBegin();
                    halo_iter != this->mpDistributedMesh->GetHaloNodeIteratorEnd();
                    ++halo_iter)
            {
                c_vector<double, SPACE_DIM> current_item = (*halo_iter)->rGetLocation();
                if (SPACE_DIM == 3)
                {
                    p_pts->InsertNextPoint(current_item[0], current_item[1], current_item[2]);
                }
                else if (SPACE_DIM == 2)
                {
                    p_pts->InsertNextPoint(current_item[0], current_item[1], 0.0);
                }
                else // (SPACE_DIM == 1)
                {
                    p_pts->InsertNextPoint(current_item[0], 0.0, 0.0);
                }
            }
        }

        mpVtkUnstructedMesh->SetPoints(p_pts);
        p_pts->Delete(); //Reference counted

        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
             elem_iter != rMesh.GetElementIteratorEnd();
             ++elem_iter)
        {

            vtkCell* p_cell=nullptr;
            ///\todo This ought to look exactly like the other MakeVtkMesh
            if (ELEMENT_DIM == 3)
            {
                p_cell = vtkTetra::New();
            }
            else if (ELEMENT_DIM == 2)
            {
                p_cell = vtkTriangle::New();
            }
            else //(ELEMENT_DIM == 1)
            {
                p_cell = vtkLine::New();
            }
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            for (unsigned j = 0; j < ELEMENT_DIM+1; ++j)
            {
                unsigned global_node_index = elem_iter->GetNodeGlobalIndex(j);
                p_cell_id_list->SetId(j, mGlobalToNodeIndexMap[global_node_index]);
            }
            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); //Reference counted
        }
        //If necessary, construct cables
        if (this->mpMixedMesh )
        {
            AugmentCellData();
            //Make a blank cell radius data for the regular elements
            std::vector<double> radii(this->mpMixedMesh->GetNumLocalElements(), 0.0);
            for (typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator elem_iter = this->mpMixedMesh->GetCableElementIteratorBegin();
                 elem_iter != this->mpMixedMesh->GetCableElementIteratorEnd();
                 ++elem_iter)
            {
                radii.push_back((*elem_iter)->GetAttribute());
                vtkCell* p_cell=vtkLine::New();
                vtkIdList* p_cell_id_list = p_cell->GetPointIds();
                for (unsigned j = 0; j < 2; ++j)
                {
                    unsigned global_node_index = (*elem_iter)->GetNodeGlobalIndex(j);
                    p_cell_id_list->SetId(j, mGlobalToNodeIndexMap[global_node_index]);
                }
                mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
                p_cell->Delete(); //Reference counted
            }
            AddCellData("Cable radius", radii);
        }


        //This block is to guard the mesh writers (vtkXMLPUnstructuredGridWriter) so that they
        //go out of scope, flush buffers and close files
        {
            assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
            vtkXMLPUnstructuredGridWriter* p_writer = vtkXMLPUnstructuredGridWriter::New();

            p_writer->SetDataModeToBinary();

            p_writer->SetNumberOfPieces(PetscTools::GetNumProcs());
            //p_writer->SetGhostLevel(-1);
            p_writer->SetStartPiece(PetscTools::GetMyRank());
            p_writer->SetEndPiece(PetscTools::GetMyRank());


#if VTK_MAJOR_VERSION >= 6
            p_writer->SetInputData(mpVtkUnstructedMesh);
#else
            p_writer->SetInput(mpVtkUnstructedMesh);
#endif
            std::string pvtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName+ ".pvtu";
            p_writer->SetFileName(pvtk_file_name.c_str());
            //p_writer->PrintSelf(std::cout, vtkIndent());
            p_writer->Write();
            p_writer->Delete(); //Reference counted
        }

        // Add provenance to the individual files
        std::stringstream filepath;
        filepath << this->mBaseName << "_" << PetscTools::GetMyRank() << ".vtu";
        AddProvenance(filepath.str());
        /// Add to the main file \todo #1494 Do we need a barrier?
        if (PetscTools::AmMaster())
        {
            AddProvenance(this->mBaseName+ ".pvtu");
        }
    }
}

// Explicit instantiation
template class VtkMeshWriter<1,1>;
template class VtkMeshWriter<1,2>;
template class VtkMeshWriter<1,3>;
template class VtkMeshWriter<2,2>; // Actually used
template class VtkMeshWriter<2,3>;
template class VtkMeshWriter<3,3>; // Actually used

#endif //CHASTE_VTK
