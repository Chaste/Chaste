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

    char* mExpectedElementType;       /**< VTK cell type (vtkTetra)*/

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
