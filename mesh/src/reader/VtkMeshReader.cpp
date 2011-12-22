/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/

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



#ifdef CHASTE_VTK

#include "VtkMeshReader.hpp"
#include "Exception.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::VtkMeshReader(std::string pathBaseName) :
    mIndexFromZero(true),
    mNumNodes(0),
    mNumElements(0),
    mNumFaces(0),
    mNodesRead(0),
    mElementsRead(0),
    mFacesRead(0),
    mBoundaryFacesRead(0),
    mNumElementAttributes(0),
    mNumFaceAttributes(0),
    mOrderOfElements(1),
    mNodesPerElement(4)
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
    mNumElements = vtk_xml_unstructured_grid_reader->GetNumberOfCells();

    // Extract the surface faces
    mpVtkGeometryFilter = vtkGeometryFilter::New();
    mpVtkGeometryFilter->SetInput(mpVtkUnstructuredGrid);
    mpVtkGeometryFilter->Update();

    mNumFaces = mpVtkGeometryFilter->GetOutput()->GetNumberOfCells();

    vtk_xml_unstructured_grid_reader->Delete();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::VtkMeshReader(vtkUnstructuredGrid* p_vtkUnstructuredGrid) :
    mIndexFromZero(true),
    mNumNodes(0),
    mNumElements(0),
    mNumFaces(0),
    mNodesRead(0),
    mElementsRead(0),
    mFacesRead(0),
    mBoundaryFacesRead(0),
    mNumElementAttributes(0),
    mNumFaceAttributes(0),
    mOrderOfElements(1),
    mNodesPerElement(4)
{
    mpVtkUnstructuredGrid = p_vtkUnstructuredGrid;

    mNumNodes = mpVtkUnstructuredGrid->GetNumberOfPoints();
    mNumElements = mpVtkUnstructuredGrid->GetNumberOfCells();

    // Extract the surface faces
    mpVtkGeometryFilter = vtkGeometryFilter::New();
    mpVtkGeometryFilter->SetInput(mpVtkUnstructuredGrid);
    mpVtkGeometryFilter->Update();

    mNumFaces = mpVtkGeometryFilter->GetOutput()->GetNumberOfCells();
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
unsigned VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNumFaceAttributes() const
{
    return mNumFaceAttributes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::Reset()
{
//    CloseFiles();
//    OpenFiles();
//    ReadHeaders();

    mNodesRead=0;
    mElementsRead=0;
    mFacesRead=0;
    mBoundaryFacesRead=0;
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

    if ( !mpVtkUnstructuredGrid->GetCell(mElementsRead)->IsA("vtkTetra") )
    {
        EXCEPTION("Element is not a vtkTetra");
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
ElementData VtkMeshReader<ELEMENT_DIM,SPACE_DIM>::GetNextFaceData()
{
    if ( mBoundaryFacesRead >= mNumFaces)
    {
        EXCEPTION( "Trying to read data for a boundary element that doesn't exist");
    }

    ElementData next_face_data;

    for (unsigned i = 0; i < 3; i++)
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
