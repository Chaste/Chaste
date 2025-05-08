/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef IMMERSEDBOUNDARYMESHWRITER_HPP_
#define IMMERSEDBOUNDARYMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMesh;

#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
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

#include "AbstractMeshWriter.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>


// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators;

/**
 * A mesh writer class for immersed boundary meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * If writing from a mesh object, the mesh to write to disk.
     * Otherwise NULL.
     */
    ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** Iterators over the mesh */
    MeshWriterIterators<ELEMENT_DIM,SPACE_DIM>* mpIters;

    /** The number of laminas in the mesh */
    unsigned mNumLaminas;

    /** Vector containing, for each element, the node indices at which it must be split for visualisation */
    std::vector<std::vector<unsigned>> mElementParts;

    // Node that this assumes SPACE_DIM == 2
    /** Shorthand for a 2d cartesian point */
    using geom_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;

    /** Shorthand for a segment of geom_points */
    using geom_segment = boost::geometry::model::segment<geom_point>;

    /** Array of geom_segments; a helper array, filled by the constructor, for GetIntersectionOfEdgeWithBoundary */
    std::array<geom_segment, 4> mBoundaryEdges;

    /**
     * When elements cross the periodic boundary they are broken up into pieces for visualisation. Some additional
     * points are needed where the element edge crosses the boundary to produce good output polygons.
     *
     * This function returns the location on the boundary where an edge (given by rStart and rEnd) crosses d[0,1]x[0,1].
     *
     * @param rStart the start point of the edge that crosses the boundary
     * @param rEnd the end point of the edge that crosses the boundary
     * @return the point at which the edge rStart->rEnd crosses d[0,1]x[0,1]
     */
    c_vector<double, SPACE_DIM> GetIntersectionOfEdgeWithBoundary(const c_vector<double, SPACE_DIM>& rStart,
                                                                  const c_vector<double, SPACE_DIM>& rEnd);

    /**
     * Helper function to get the nearest corner to the average of two points on the boundary.
     * @param rA one of the boundary points
     * @param rB the other boundary point
     * @return the coordinates of the nearest corner of [0,1]x[0,1]
     */
    c_vector<double, SPACE_DIM> GetNearestCorner(const c_vector<double, SPACE_DIM>& rA,
                                                 const c_vector<double, SPACE_DIM>& rB) const;

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
///\todo Merge into VtkMeshWriter (#1076)
    vtkUnstructuredGrid* mpVtkUnstructedMesh;
#endif //CHASTE_VTK

public:

    /**
     * Constructor.
     *
     * @param rDirectory reference to the output directory, relative to where Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to writing files (defaults to true)
     */
    ImmersedBoundaryMeshWriter(const std::string& rDirectory,
                               const std::string& rBaseName,
                               bool clearOutputDir=true);

    /**
     * Destructor.
     */
    ~ImmersedBoundaryMeshWriter();

    /**
     * Write files using a mesh.
     *
     * @param rMesh reference to the vertex-based mesh
     */
    void WriteFilesUsingMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Write VTK file using a mesh.
     *
     * @param rMesh reference to the immersed boundary mesh
     * @param stamp is an optional stamp (like a time-stamp) to put into the name of the file
     */
    void WriteVtkUsingMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp="");

    /**
     * Populate mpVtkUnstructedMesh using a vertex-based mesh.
     * Called by WriteVtkUsingMesh().
     *
     * @param rMesh reference to the vertex-based mesh
     */
    void MakeVtkMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

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
     * @return the coordinates of the next node to be written to file
     */
    std::vector<double> GetNextNode();

    /**
     * @return the data (indices/attributes) of the next element to be written to file
     */
    ImmersedBoundaryElementData GetNextImmersedBoundaryElement();

    /**
     * @return the data (indices/attributes) of the next lamina to be written to file
     */
    ImmersedBoundaryElementData GetNextImmersedBoundaryLamina();

    /**
     * Write mesh data to files.
     * This method must be overridden in concrete classes.
     */
    void WriteFiles();

    /**
     * Analyses the mesh to determine where (if at all) each element overlaps due to the periodic boundaries.
     *
     * @param rMesh the immersed boundary mesh
     */
    void FindElementOverlaps(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /** @return reference to mElementParts */
    const std::vector<std::vector<unsigned>>& rGetElementParts() const;
};

#endif /*IMMERSEDBOUNDARYMESHWRITER_HPP_*/
