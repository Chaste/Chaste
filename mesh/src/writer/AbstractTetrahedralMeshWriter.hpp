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


#ifndef _ABSTRACTTETRAHEDRALMESHWRITER_HPP_
#define _ABSTRACTTETRAHEDRALMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistributedTetrahedralMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MixedDimensionMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators;

#include <fstream>
#include <sstream>
#include <boost/scoped_array.hpp>

#include "AbstractMeshWriter.hpp"
#include "AbstractMesh.hpp"
#include "NodeMap.hpp"

#include "GenericEventHandler.hpp"
/** \todo Temporary for #2351 */
class MeshEventHandler : public GenericEventHandler<11, MeshEventHandler>
{
public:
    static const char* EventName[11]; /**< See #2351*/

    /** See #2351*/
    typedef enum
    {
        TRIANGLES=0,
        BINTRI,
        VTK,
        PVTK,
        NODE,
        ELE,
        FACE,
        NCL,
        COMM1,
        COMM2
    } EventType;
};

/**
 * An abstract tetrahedral mesh writer class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * Post an element from a slave process to the master for concentration
     * @param globalIndex index of this element
     * @param indices the node index data for this element
     * @param numIndices number of indices
     * @param tag temporary
     * @param attribute
     */
    void PostElement(unsigned globalIndex, unsigned indices[], unsigned numIndices, unsigned tag, double attribute)
    {
        MeshEventHandler::BeginEvent(MeshEventHandler::COMM1);
        MPI_Ssend(indices, numIndices, MPI_UNSIGNED, 0,
                  tag, //Elements sent with tags offset
                  PETSC_COMM_WORLD);
        MeshEventHandler::EndEvent(MeshEventHandler::COMM1);
        // Attribute value has the same tag (assume that it doesn't overtake the previous message)
        MeshEventHandler::BeginEvent(MeshEventHandler::COMM2);
        MPI_Ssend(&attribute, 1, MPI_DOUBLE, 0,
                  tag, //Elements sent with tags offset
                  PETSC_COMM_WORLD);
        MeshEventHandler::EndEvent(MeshEventHandler::COMM2);
    }
    /**
     * Unpack an element from a slave process on the master for concentration
     * @param rElementData the output structure to fill (should have the NodeIndices structure of the correct size
     * @param globalIndex index of this element
     * @param numIndices number of indices
     * @param tag temporary
     */
    void UnpackElement(ElementData& rElementData, unsigned globalIndex, unsigned numIndices, unsigned tag)
    {
        assert( numIndices == rElementData.NodeIndices.size());
        boost::scoped_array<unsigned> raw_indices(new unsigned[numIndices]);
        MPI_Status status;
        // Get it from elsewhere
        MeshEventHandler::BeginEvent(MeshEventHandler::COMM1);
        MPI_Recv(raw_indices.get(), numIndices, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 tag,
                 PETSC_COMM_WORLD, &status);
        MeshEventHandler::EndEvent(MeshEventHandler::COMM1);
        // Convert to std::vector
        for (unsigned j=0; j< rElementData.NodeIndices.size(); j++)
        {
            rElementData.NodeIndices[j] = raw_indices[j];
        }

        // Attribute value has the same tag (assume that it doesn't overtake the previous message)
        double attribute;
        MeshEventHandler::BeginEvent(MeshEventHandler::COMM2);
        MPI_Recv(&attribute, 1U, MPI_DOUBLE, MPI_ANY_SOURCE,
                 tag,
                 PETSC_COMM_WORLD, &status);
        MeshEventHandler::EndEvent(MeshEventHandler::COMM2);

        // Attribute value
        rElementData.AttributeValue = attribute;
    }

    /**
     * Write a parallel mesh to file. Used by the serialization methods.
     *
     * @param keepOriginalElementIndexing  Whether to write the mesh with the same element ordering as in memory.
     *                                     Optimisations can be applied if this is not needed.
     */
    virtual void WriteFilesUsingParallelMesh(bool keepOriginalElementIndexing=true);

    /**
     * Write out a node connectivity information file (collectively called,
     * involves some communication).
     *
     * @param rMesh the mesh object in memory
     * @param invertMeshPermutation  whether to permute the file using the inverse of the mesh's permutation
     */
    void WriteNclFile(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                      bool invertMeshPermutation=false);

    /**
     * Create output files and add headers.
     */
    virtual void CreateFilesWithHeaders();

    /**
     * Append local mesh data to output files.
     */
    virtual void AppendLocalDataToFiles();

    /**
     * Append footers to output files.
     */
    virtual void WriteFilesFooter();

    NodeMap* mpNodeMap; /**<Node map to be used when writing a mesh that has deleted nodes*/

protected:

    unsigned mNodesPerElement; /**< Same as (ELEMENT_DIM+1), except when writing a quadratic mesh!*/
    unsigned mNodesPerBoundaryElement; /**< Same as (ELEMENT_DIM), except when writing a quadratic mesh!*/

    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh; /**<Pointer to the mesh (if we are writing from a mesh)*/
    DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpDistributedMesh; /**< Another pointer to the mesh, produced by dynamic cast*/
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* mpMixedMesh; /**< Another pointer to the mesh, produced by dynamic cast*/
    MeshWriterIterators<ELEMENT_DIM,SPACE_DIM>* mpIters; /**< Handy iterators so that we know the next node/element to be written */

    bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
    bool mWriteMetaFile; /**< Whether to write a metafile (only used by MeshylazerMeshWriter) */
    unsigned mNodeCounterForParallelMesh; /**< Used by master process for polling processes for the next node */
    unsigned mElementCounterForParallelMesh;/**< Used by master process for polling processes for the next element */
    unsigned mBoundaryElementCounterForParallelMesh;/**< Used by master process for polling processes for the next boundary element */
    unsigned mCableElementCounterForParallelMesh;/**< Used by master process for polling processes for the next cable element */
    bool mFilesAreBinary;  /**< Whether all data is to be written as binary - used in derived class TrianglesMeshWriter*/

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    AbstractTetrahedralMeshWriter(const std::string& rDirectory,
                       const std::string& rBaseName,
                       const bool clearOutputDir=true);

    /**
     *  Destructor just deletes the node map if memory has been allocated for it
     */
    virtual ~AbstractTetrahedralMeshWriter();

    /**
     * Write a const mesh to file. Used by the serialization methods and avoids iterators...
     *
     * @param rMesh the mesh
     * @param keepOriginalElementIndexing  Whether to write the mesh with the same element ordering as in memory.
     *                                     Optimisations can be applied if this is not needed.
     *
     * \todo #1322 Mesh should really be const!
     */
    virtual void WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                     bool keepOriginalElementIndexing=true);

    /**
     * Write a const mesh to file. Used by the serialization methods and avoids iterators.
     * The master process will use the mesh reader to copy over most of the mesh files,
     * converting them to binary format. However, the mesh is used to write a .ncl file
     * which contains node connectivity information.
     *
     * @param rMeshReader a reader of the original mesh files on disk
     * @param rMesh the mesh object in memory
     */
    void WriteFilesUsingMeshReaderAndMesh(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                          AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

   /**
     * @return the coordinates of the next node to be written to file
     */
    std::vector<double> GetNextNode();


    /**
     * @return the data (indices/attributes) of the next element to be written to file
     */
    ElementData GetNextElement();

    /**
     * @return the data (indices) of the next boundary element to be written to file
     */
    ElementData GetNextBoundaryElement();

    /**
     * @return the data (indices/attributes) of the next cable element to be written to file
     */
    ElementData GetNextCableElement();
};

#endif //_ABSTRACTTETRAHEDRALMESHWRITER_HPP_
