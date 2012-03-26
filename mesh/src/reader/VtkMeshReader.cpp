/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/

/*

Copyright (c) 2005-2012, University of Oxford.
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



#ifdef CHASTE_VTK

#include "VtkMeshReader.hpp"
#include "Exception.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::VtkMeshReader(std::string pathBaseName) :
    mIndexFromZero(true),
    mNumNodes(0),
    mNumElements(0),
    mNumFaces(0),
    mNumCableElements(0),
    mNodesRead(0),
    mElementsRead(0),
    mFacesRead(0),
    mBoundaryFacesRead(0),
    mCableElementsRead(0),
    mNumElementAttributes(0),
    mNumFaceAttributes(0),
    mNumCableElementAttributes(0),
    mOrderOfElements(1),
    mNodesPerElement(4),
    mVtkCellType(VTK_TETRA)
{
    vtkXMLUnstructuredGridReader* vtk_xml_unstructured_grid_reader;

    // Check file exists
    mVtuFile.open(pathBaseName.c_str());
    if ( !mVtuFile.is_open() )
    {
        EXCEPTION("Could not open VTU file: " + pathBaseName);
    }
    mVtuFile.close();

    // Load the mesh geometry and data from a file
    vtk_xml_unstructured_grid_reader = vtkXMLUnstructuredGridReader::New();
    vtk_xml_unstructured_grid_reader->SetFileName( pathBaseName.c_str() );
    vtk_xml_unstructured_grid_reader->Update();
    mpVtkUnstructuredGrid = vtk_xml_unstructured_grid_reader->GetOutput();

    mNumNodes = vtk_xml_unstructured_grid_reader->GetNumberOfPoints();

    unsigned num_cells = vtk_xml_unstructured_grid_reader->GetNumberOfCells();


    if (SPACE_DIM == 2)
    {
        mNodesPerElement = 3;
        mVtkCellType = VTK_TRIANGLE;
        mNumFaces = 0; ///\todo Make the above filter work in 2D
    }
    //Determine if we have multiple cell types - such as cable elements in addition to tets/triangles
    vtkCellTypes* cell_types = vtkCellTypes::New();
    mpVtkUnstructuredGrid->GetCellTypes(cell_types);

    if(cell_types->GetNumberOfTypes() > 1)
    {
        for(unsigned cell_id = 0; cell_id < num_cells; ++cell_id)
        {
            if(mpVtkUnstructuredGrid->GetCellType(cell_id) == mVtkCellType)
            {
                ++mNumElements;
                assert(mNumCableElements == 0); //We expect all the simplices first and then the cables at the end of the array
            }
            else if(mpVtkUnstructuredGrid->GetCellType(cell_id) == VTK_LINE)
            {
                ++mNumCableElements;
            }
            else
            {
                NEVER_REACHED;
            }
        }
    }
    else
    {
        //There is only 1 cell type, so all cells are elements
        mNumElements = num_cells;
    }

    cell_types->Delete();


    // Extract the surface faces
    mpVtkGeometryFilter = vtkGeometryFilter::New();
    mpVtkGeometryFilter->SetInput(mpVtkUnstructuredGrid);
    mpVtkGeometryFilter->Update();
    mNumFaces = mpVtkGeometryFilter->GetOutput()->GetNumberOfCells();
    ///\todo #2052 The boundary face filter includes the cable elements
    if (mNumCableElements > 0)
    {
        //Stop reading faces
        mNumFaces =0;
    }


    vtk_xml_unstructured_grid_reader->Delete();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::VtkMeshReader(vtkUnstructuredGrid* p_vtkUnstructuredGrid) :
    mIndexFromZero(true),
    mNumNodes(0),
    mNumElements(0),
    mNumFaces(0),
    mNumCableElements(0),
    mNodesRead(0),
    mElementsRead(0),
    mFacesRead(0),
    mBoundaryFacesRead(0),
    mCableElementsRead(0),
    mNumElementAttributes(0),
    mNumFaceAttributes(0),
    mNumCableElementAttributes(0),
    mOrderOfElements(1),
    mNodesPerElement(4),
    mVtkCellType(VTK_TETRA)
{
    mpVtkUnstructuredGrid = p_vtkUnstructuredGrid;

    mNumNodes = mpVtkUnstructuredGrid->GetNumberOfPoints();
    mNumElements = mpVtkUnstructuredGrid->GetNumberOfCells();

    // Extract the surface faces
    mpVtkGeometryFilter = vtkGeometryFilter::New();
    mpVtkGeometryFilter->SetInput(mpVtkUnstructuredGrid);
    mpVtkGeometryFilter->Update();

    mNumFaces = mpVtkGeometryFilter->GetOutput()->GetNumberOfCells();
    assert (SPACE_DIM == 3);
//    if (SPACE_DIM == 2)
//    {
//        mNodesPerElement = 3;
//        mVtkCellType = VTK_TRIANGLE;
//        mNumFaces = 0; ///\todo Make the above filter work in 2D
//    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::~VtkMeshReader()
{
    mpVtkGeometryFilter->Delete();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumElements() const
{
    return mNumElements;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumCableElements() const
{
    return mNumCableElements;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumFaces() const
{
    return mNumFaces;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumEdges() const
{
    return mNumFaces;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumElementAttributes() const
{
    return mNumElementAttributes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumCableElementAttributes() const
{
    return mNumCableElementAttributes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumFaceAttributes() const
{
    return mNumFaceAttributes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::Reset()
{
    mNodesRead=0;
    mElementsRead=0;
    mFacesRead=0;
    mBoundaryFacesRead=0;
    mCableElementsRead=0;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::Initialize()
{
    mpVtkUnstructuredGrid->Initialize();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNextNode()
{
    if ( mNodesRead >= mNumNodes )
    {
        EXCEPTION( "Trying to read data for a node that doesn't exist" );
    }

    std::vector<double> next_node;

    for (unsigned i = 0; i < 3; i++)
    {
        next_node.push_back( mpVtkUnstructuredGrid->GetPoint(mNodesRead)[i] );
    }

    mNodesRead++;
    return next_node;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNextElementData()
{
    if ( mElementsRead >= mNumElements )
    {
        EXCEPTION( "Trying to read data for an element that doesn't exist" );
    }

    if (mpVtkUnstructuredGrid->GetCellType(mElementsRead)!=  mVtkCellType)
    {
        EXCEPTION("Element is not of expected type (vtkTetra/vtkTriangle)");
    }

    ElementData next_element_data;

    for (unsigned i = 0; i < mNodesPerElement; i++)
    {
        next_element_data.NodeIndices.push_back(mpVtkUnstructuredGrid->GetCell(mElementsRead)->GetPointId(i));
    }

    // \todo implement method to read element data properly (currently returns zero always...)
    next_element_data.AttributeValue = 0;

    mElementsRead++;
    return next_element_data;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNextCableElementData()
{
    if ( mCableElementsRead >=  mNumCableElements )
    {
        EXCEPTION( "Trying to read data for a cable element that doesn't exist" );
    }

    unsigned next_index = mNumElements + mCableElementsRead;
    assert(mpVtkUnstructuredGrid->GetCellType(next_index)==VTK_LINE);

    ElementData next_element_data;

    for (unsigned i = 0; i < 2; i++)
    {
        next_element_data.NodeIndices.push_back(mpVtkUnstructuredGrid->GetCell(next_index)->GetPointId(i));
    }

    /// \todo #2052 Look at the radius data too
    next_element_data.AttributeValue = 0;

    mCableElementsRead++;
    return next_element_data;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNextFaceData()
{
    if ( mBoundaryFacesRead >= mNumFaces)
    {
        EXCEPTION( "Trying to read data for a boundary element that doesn't exist");
    }

    ElementData next_face_data;

    for (unsigned i = 0; i < (mNodesPerElement-1); i++)
    {
        next_face_data.NodeIndices.push_back(mpVtkGeometryFilter->GetOutput()->GetCell(mBoundaryFacesRead)->GetPointId(i));
    }

    mBoundaryFacesRead++;
    return next_face_data;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetCellData(std::string dataName, std::vector<double>& dataPayload)
{
    vtkCellData *p_cell_data = mpVtkUnstructuredGrid->GetCellData();

    if ( !p_cell_data->HasArray(dataName.c_str()) )
    {
        EXCEPTION("No cell data '" + dataName + "'");
    }

    vtkDataArray *p_scalars = p_cell_data->GetArray( dataName.c_str() );
    if (p_scalars->GetNumberOfComponents() != 1)
    {
        EXCEPTION("The cell data '" + dataName + "' is not scalar data.");
    }

    dataPayload.clear();
    for (unsigned i = 0; i < mNumElements; i++)
    {
        dataPayload.push_back( p_scalars->GetTuple(i)[0] );
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetCellData(std::string dataName, std::vector<c_vector<double,SPACE_DIM> >& dataPayload)
{
    vtkCellData *p_cell_data = mpVtkUnstructuredGrid->GetCellData();

    if ( !p_cell_data->HasArray(dataName.c_str()) )
    {
        EXCEPTION("No cell data '" + dataName + "'");
    }

    vtkDataArray *p_scalars = p_cell_data->GetArray( dataName.c_str() );
    if (p_scalars->GetNumberOfComponents() != 3)
    {
        EXCEPTION("The cell data '" + dataName + "' is not 3-vector data.");
    }

    dataPayload.clear();
    for (unsigned i = 0; i < mNumElements; i++)
    {
        c_vector <double, SPACE_DIM> data;
        for (unsigned j = 0; j < SPACE_DIM; j++)
        {
            data[j]=  p_scalars->GetTuple(i)[j];
        }
        dataPayload.push_back(data);
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetPointData(std::string dataName, std::vector<double>& dataPayload)
{
    vtkPointData *p_point_data = mpVtkUnstructuredGrid->GetPointData();

    if ( !p_point_data->HasArray(dataName.c_str()) )
    {
        EXCEPTION("No point data '" + dataName + "'");
    }

    vtkDataArray *p_scalars = p_point_data->GetArray( dataName.c_str() );
    if (p_scalars->GetNumberOfComponents() != 1)
    {
        EXCEPTION("The point data '" + dataName + "' is not scalar data.");
    }

    dataPayload.clear();

    for (unsigned i = 0; i < mNumNodes; i++)
    {
        dataPayload.push_back( p_scalars->GetTuple(i)[0] );
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetPointData(std::string dataName, std::vector<c_vector<double,SPACE_DIM> >& dataPayload)
{
   vtkPointData *p_point_data = mpVtkUnstructuredGrid->GetPointData();

    if ( !p_point_data->HasArray(dataName.c_str()) )
    {
        EXCEPTION("No point data '" + dataName + "'");
    }

    vtkDataArray *p_scalars = p_point_data->GetArray( dataName.c_str() );

    if (p_scalars->GetNumberOfComponents() != 3)
    {
        EXCEPTION("The point data '" + dataName + "' is not 3-vector data.");
    }
    dataPayload.clear();
    for (unsigned i = 0; i < mNumNodes; i++)
    {
        c_vector<double, SPACE_DIM> data;
        for (unsigned j = 0; j< SPACE_DIM; j++)
        {
            data[j] = p_scalars->GetTuple(i)[j];
        }
        dataPayload.push_back( data );
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
vtkUnstructuredGrid* VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::OutputMeshAsVtkUnstructuredGrid()
{
    return mpVtkUnstructuredGrid;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VtkMeshReader<1,1>;
template class VtkMeshReader<1,2>;
template class VtkMeshReader<1,3>;
template class VtkMeshReader<2,2>;
template class VtkMeshReader<2,3>;
template class VtkMeshReader<3,3>;
#endif // CHASTE_VTK
