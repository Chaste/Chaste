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

#ifndef VTKMESHWRITER_HPP_
#define VTKMESHWRITER_HPP_

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkDataCompressor.h>
#include "AbstractTetrahedralMeshWriter.hpp"
#include "Version.hpp"


#include <map>

// Forward declaration prevents circular include chain
template<unsigned SPACE_DIM>
class NodesOnlyMesh;

/**
 *  VtkMeshWriter
 *
 *  Writes a mesh in VTK .vtu format (that's an XML-based, data compressed unstructured mesh)
 *
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VtkMeshWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>
{

//Requires  "sudo aptitude install libvtk5-dev" or similar

private:
    bool mWriteParallelFiles; /**< Whether to write parallel (.pvtu + .vtu for each process) files, defaults to false */

    std::map<unsigned, unsigned> mGlobalToNodeIndexMap; /**< Map a global node index into a local index (into mNodes and mHaloNodes as if they were concatenated) */

    std::vector<std::vector<unsigned> > mNodesToSendPerProcess; /**< Used to communicate node-wise halo data */
    std::vector<std::vector<unsigned> > mNodesToReceivePerProcess;  /**< Used to communicate node-wise halo data */

    /** A pointer to a NodesOnlyMesh to write to file, created by dynamic cast */
    NodesOnlyMesh<SPACE_DIM>* mpNodesOnlyMesh;

    /**
     * A VTK mesh data structure.
     * Created at construction, has data associated with it by AddCellData
     * and AddCellPoint, then is filled with mesh geometry by MakeVtkMesh() in
     * WriteFiles().
     */
    vtkUnstructuredGrid* mpVtkUnstructedMesh;

    /**
     * Private helper method which copies the mesh details into the waiting
     * VTK mesh structure.  Called by  WriteFiles().
     */
    void MakeVtkMesh();

    /**
     * At the time of adding VTK cell data, it is assumed that there is one piece of data
     * for each element in the original mesh.  If the mesh is mixed-dimension (elements and
     * cable elements) the VTK mesh makes no distinction between the two types of cells.
     * All data associated with cells must be the same length as the overall number of cells.
     * This method inspects each cell data component and adds dummy data to cover the cable elements.
     */
    void AugmentCellData();

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param rCleanDirectory  whether to clean the directory (defaults to true)
     */
    VtkMeshWriter(const std::string& rDirectory, const std::string& rBaseName, const bool& rCleanDirectory=true);

    /**
     * Write mesh data to files.
     */
    void WriteFiles();


    /**
      * Add a scalar data field to each element (known as "cell" in VTK).
      * @param name is a meaningful name with which to annotate the data
      * @param data is the data which should appear in the same order as the element numbering
      * The length of the data vector is assumed to match the number of elements in the mesh.
      * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
      */
     void AddCellData(std::string name, std::vector<double> data);

     /**
     * Add a vector data field to each element (known as "cell" in VTK).
     * @param name is a meaningful name with which to annotate the data
     * @param data is the data which should appear in the same order as the element numbering
     * The length of the data vector is assumed to match the number of elements in the mesh.
     * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
     */
    void AddCellData(std::string name, std::vector<c_vector<double, SPACE_DIM> > data);

    /**
     * Add a symmetric tensor data field to each element (known as "cell" in VTK).
     * @param name is a meaningful name with which to annotate the data
     * @param data is the data which should appear in the same order as the element numbering
     * The length of the data vector is assumed to match the number of elements in the mesh.
     * The data vector represents the lower half of the tensor
     * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
     */
    void AddTensorCellData(std::string name, std::vector<c_vector<double,SPACE_DIM*(SPACE_DIM+1)/2> > data);

    /**
     * Add a tensor data field to each element (known as "cell" in VTK).
     * @param name is a meaningful name with which to annotate the data
     * @param data is the data which should appear in the same order as the element numbering
     * The length of the data vector is assumed to match the number of elements in the mesh.
     * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
     */
    void AddTensorCellData(std::string name, std::vector<c_matrix<double,SPACE_DIM,SPACE_DIM> > data);


    /**
     * Add a scalar data field to each node (known as "point" in VTK).
     * @param name is a meaningful name with which to annotate the data
     * @param data is the data which should appear in the same order as the node numbering
     * The length of the data vector is assumed to match the number of nodes in the mesh
     * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
     */
    void AddPointData(std::string name, std::vector<double> data);

    /**
     * Add a vector data field to each node (known as "point" in VTK).
     * @param name is a meaningful name with which to annotate the data
     * @param data is the data which should appear in the same order as the node numbering
     * The length of the data vector is assumed to match the number of nodes in the mesh
     * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
     */
    void AddPointData(std::string name, std::vector<c_vector<double, SPACE_DIM> > data);

    /**
     * Add a tensor data field to each point.
     * @param name is a meaningful name with which to annotate the data
     * @param data is the data which should appear in the same order as the node numbering
     * The length of the data vector is assumed to match the number of nodes in the mesh
     * Checking cannot be done at this stage since the data is associated with an empty VTK mesh structure.
     */
    void AddTensorPointData(std::string name, std::vector<c_matrix<double,SPACE_DIM,SPACE_DIM> > data);

    /**
     * Should be called to enable files to be written in parallel (i.e. a .pvtu file and .vtu files for each
     * process's sub-mesh).
     *
     * @param rMesh the mesh (must be a DistributedTetrahedralMesh)
     */
     void SetParallelFiles(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh);

    /**
     * Write files. Overrides the method implemented in AbstractTetrahedralMeshWriter, which concentrates mesh
     * data onto a single file in order to output a monolithic file. For VTK, a DistributedTetrahedralMesh in
     * parallel is instead written out as a set of .vtu files (one for each sub-mesh) and a .pvtu file that
     * provides the visualizer with information about them.
     *
     * @param rMesh the mesh
     * @param keepOriginalElementIndexing  Whether to write the mesh with the same element ordering.
     *                                     Optimisations can be applied if this is not needed.
     *
     * \todo #1322 Mesh should really be const!
     */
    void WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                     bool keepOriginalElementIndexing=true);

    /**
     * Add Chaste provenance data to a VTK file as an XML comment string
     * @param fileName is the file name relative to mpOutputFileHandler
     * The file is assumed have been written to and to be closed - so that it can safely be appended to.
     */
    void AddProvenance(std::string fileName);

    /**
     * Destructor.
     */
    virtual ~VtkMeshWriter();
};

#endif //CHASTE_VTK

#endif /*VTKMESHWRITER_HPP_*/
