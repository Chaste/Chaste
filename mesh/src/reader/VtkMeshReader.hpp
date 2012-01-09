/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/

/*

Copyright (C) University of Oxford, 2005-2012

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


#ifndef VTKMESHREADER_HPP_
#define VTKMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericGeometryFilter.h>
#include <vtkDataCompressor.h>

#include "UblasVectorInclude.hpp"
#include "AbstractMeshReader.hpp"

/**
 *  VtkMeshReader
 *
 *  Reads a mesh from VTK .vtu format (that's an XML-based, data compressed unstructured mesh)
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VtkMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** vtkUnstructuredGrid object: the rest of the class acts as an interface to this */
    vtkUnstructuredGrid* mpVtkUnstructuredGrid;

    /** vtkGeometry filter object to extract the surface elements */
    vtkGeometryFilter* mpVtkGeometryFilter;

    bool mIndexFromZero;             /**< True if input data is numbered from zero, false otherwise */

    std::ifstream mVtuFile;            /**< Location of the .vtu file */

    unsigned mNumNodes;                /**< Number of nodes in the mesh */
    unsigned mNumElements;            /**< Number of elements in the mesh */
    unsigned mNumFaces;                /**< Number of faces in the mesh */

    unsigned mNodesRead;            /**< Number of nodes read from file so far */
    unsigned mElementsRead;            /**< Number of elements read from file so far */
    unsigned mFacesRead;            /**< Number of faces read from file so far */
    unsigned mBoundaryFacesRead;    /**< Number of boundary faces read from file so far */

    unsigned mNumNodeAttributes;     /**< Is the number of attributes stored at each node */
    unsigned mMaxNodeBdyMarker;     /**< Is the maximum node boundary marker */
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element */
    unsigned mNumFaceAttributes;     /**< Is the number of attributes stored for each face */

    unsigned mOrderOfElements;        /**< Order of the elements (i.e. linear, quadratic, cubic FE basis functions */
    unsigned mNodesPerElement;        /**< Number of nodes per element */

public:

    /**
     * Constructor
     *
     * @param pathBaseName Full file path of the input file
     */
    VtkMeshReader(std::string pathBaseName);

    /**
     * Alternative constructor, takes a vtkUnstructuredGrid that is already in memory as an input
     * parameter rather than a .vtu file
     *
     * @param p_vtkUnstructuredGrid Pointer to a vtkUnstructuredGrid object
     */
    VtkMeshReader(vtkUnstructuredGrid* p_vtkUnstructuredGrid);

    /**
     * Destructor
     */
    virtual ~VtkMeshReader();

    /**
     * Returns the mNumElements
     */
    unsigned GetNumElements() const;

    /**
     * Returns mNumNodes
     */
    unsigned GetNumNodes() const;

    /**
     * Returns mNumFaces (synonym of GetNumEdges() method)
     */
    unsigned GetNumFaces() const;

    /**
     * Returns mNumFaces (synonym of GetNumFaces() method)
     */
    unsigned GetNumEdges() const;

    /**
     * Returns mNumElementAttributes
     */
    unsigned GetNumElementAttributes() const;

    /**
     * Returns mNumFaceAttributes
     */
    unsigned GetNumFaceAttributes() const;

    /**
     * Resets mNodesRead, mElementsRead, mFacesRead and mBoundaryFacesRead to zero (for another pass through
     * the mesh from the beginning
     */
    void Reset();

    /**
     * Deletes the vtkUnstructuredGrid and vtkGeometryFilter in preparation for deletion of the mesh reader (should)
     * not be called if the vtkUnstructuredGrid is still required elsewhere, e.g. in an AdaptiveTetrahedralMesh or an
     * AdaptiveBidomainProblem.
     */
    void Initialize();

    /**
     * Returns a vector of the coordinates of each node in turn
     */
    std::vector<double> GetNextNode();

    /**
     * Returns a vector of the nodes of each element (and any attribute infomation, if there is any) in turn
     */
    ElementData GetNextElementData();

    /**
     * Returns a vector of the nodes of each face in turn (synonym of GetNextEdgeData())
     */
    ElementData GetNextFaceData();


    /**
     * Returns an std::vector containing the vtkCellData with attribute name specified
     * Throws if the attribute name does not exist
     * @param dataName Name of the cell data
     * @param dataPayload in which to store the result
     */
     void GetCellData(std::string dataName, std::vector<double>& dataPayload);

    /**
     * Returns an std::vector containing the vtkPointData with attribute name specified
     * Throws if the attribute name does not exist
     *
     * @param dataName Name of the point data
     * @param dataPayload in which to store the result
     */
    void GetPointData(std::string dataName, std::vector<double>& dataPayload);

    /**
     * Returns an std::vector containing the vector-directed vtkCellData with attribute name specified
     * Throws if the attribute name does not exist
     * @param dataName Name of the cell data
     * @param dataPayload in which to store the result
     */
     void GetCellData(std::string dataName, std::vector<c_vector<double,SPACE_DIM> >& dataPayload);

    /**
     * Returns an std::vector containing the vector-directed vtkPointData with attribute name specified
     * Throws if the attribute name does not exist
     *
     * @param dataName Name of the point data
     * @param dataPayload in which to store the result
     */
    void GetPointData(std::string dataName, std::vector<c_vector<double,SPACE_DIM> >& dataPayload);

    /**
     * Return a pointer to mpVtkUnstructuredGrid
     */
    vtkUnstructuredGrid* OutputMeshAsVtkUnstructuredGrid();

};

#endif/*CHASTE_VTK*/

#endif/*VTKMESHREADER_HPP_*/
