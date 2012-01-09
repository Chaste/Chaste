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
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkDataCompressor.h>
#include "AbstractTetrahedralMeshWriter.hpp"
#include "Version.hpp"

#include <map>

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
