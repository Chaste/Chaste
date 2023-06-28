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

#ifndef SEMMESHWRITER_HPP_
#define SEMMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned DIM>
class SemMesh;

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkConvexPointSet.h>
#include <vtkPolygon.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataCompressor.h>
#endif //CHASTE_VTK

#include "SemMesh.hpp"
#include "AbstractMeshWriter.hpp"
#include "NodeMap.hpp"

// Forward declaration prevents circular include chain
template<unsigned DIM>
class SemMesh;

template<unsigned DIM>
struct MeshWriterIterators;

/**
 * A mesh writer class for the SemMesh class.
 */
template<unsigned DIM>
class SemMeshWriter : public AbstractMeshWriter<DIM, DIM>
{
private:
    /**
     * If writing from a mesh object, the mesh to write to disk.
     * Otherwise NULL.
     */
    SemMesh<DIM>* mpMesh;

    /** Iterators over the mesh. */
    MeshWriterIterators<DIM>* mpIters;

    /** Track deleted nodes so they don't get written. */
    NodeMap* mpNodeMap;

    /** The last index written to mpNodeMap. */
    unsigned mNodeMapCurrentIndex;

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
///\todo Merge into VtkMeshWriter (#1076)
    vtkUnstructuredGrid* mpVtkUnstructedMesh;
#endif //CHASTE_VTK

public:

    /**
     * Constructor.
     *
     * @param rDirectory reference to the output directory, relative to where 
     *                   Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to 
     *                       riting files
     */
    SemMeshWriter(const std::string& rDirectory,
                  const std::string& rBaseName,
                  const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    ~SemMeshWriter();

    ///\todo Mesh should be const
    /**
     * Write files using a mesh.
     *
     * @param rMesh reference to the SemMesh
     */
    void WriteFilesUsingMesh(SemMesh<DIM>& rMesh);

    /**
     * Write VTK file using a mesh.
     *
     * @param rMesh reference to the SemMesh
     * @param stamp optional stamp (e.g. a time-stamp) to put into the filename
     */
    void WriteVtkUsingMesh(SemMesh<DIM>& rMesh, std::string stamp="");

    /**
     * Populate mpVtkUnstructedMesh using a SemMesh.
     * Called by WriteVtkUsingMesh().
     *
     * @param rMesh reference to the SemMesh
     */
    void MakeVtkMesh(SemMesh<DIM>& rMesh);

    /**
     * Add data to a future VTK file.
     *
     * @param dataName a tag to go into the VTK file
     * @param dataPayload a pay-load of length (number of elements)
     */
    void AddCellData(std::string dataName, std::vector<double> dataPayload);

    /**
     * Add data to a future VTK file.
     *
     * @param dataName a tag to go into the VTK file
     * @param dataPayload a pay-load of length (number of nodes)
     */
    void AddPointData(std::string dataName, std::vector<double> dataPayload);

    /**
     * @return the coordinates of the next node to be written to file.
     */
    std::vector<double> GetNextNode();

    /**
     * @return the data (indices/attributes) of the next element to be written 
     *         to file.
     */
    ElementData GetNextElement();

    /**
     * Overriden method to write mesh data to files.
     */
    void WriteFiles();
};

#endif /*SEMMESHWRITER_HPP_*/
