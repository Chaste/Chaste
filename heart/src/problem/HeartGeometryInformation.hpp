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
#ifndef HEARTGEOMETRYINFORMATION_HPP_
#define HEARTGEOMETRYINFORMATION_HPP_

#include <vector>
#include <string>
#include <set>
#include "DistanceMapCalculator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "ChasteCuboid.hpp"

/** Names for layers in the heart wall */
typedef enum HeartLayerType_
{
    ENDO = 0,
    MID,
    EPI
} HeartLayerType;

/** Type for region codes (LEFT/RIGHT VENTRICLE or SEPTUM) walls and surfaces */
typedef unsigned HeartRegionType;


/**
 * This class provides a method to calculate the relative position of a node with respect to two (or three)
 * given surfaces.
 */
template<unsigned SPACE_DIM>
class HeartGeometryInformation
{
private:

    /** Area of the septum considered to belong to the left septum (relative to 1)*/
    static const double LEFT_SEPTUM_SIZE;
    /** Area of the septum considered to belong to the right septum (relative to 1)*/
    static const double RIGHT_SEPTUM_SIZE;

    /** The nodes on the epicardial surface */
    std::vector<unsigned> mEpiSurface;

    /** The nodes on the endocardial surface (only used in the 2 surface case) */
    std::vector<unsigned> mEndoSurface;

    /** The nodes on the endocardial left ventricular surface (only used in the 3 surface case) */
    std::vector<unsigned> mLVSurface;

    /** The nodes on the endocardial right ventricular surface (only used in the 3 surface case) */
    std::vector<unsigned> mRVSurface;

    /**
     *  Takes in a file of all the nodes on ONE PARTICULAR surface of the
     *  mesh (eg the right ventricular endo-cardial surface) and collects all the nodes
     *  on that surface in one vector.
     *
     *  @param rSurfaceFileName  The surface file, lists global node indices on this surface.
     *         The number of lines in this file and entries per line doesn't matter. Just have to
     *         be separated by whitespace and returns. But these files would tend to be:
     *         EITHER: multiple entries on one line -- all nodes in each boundary element,
     *         OR: simply a list of the nodes on the surface.
     *  @param rSurfaceNodes  The returned vector of nodes indices on this surface
     *  @param indexFromZero  True for native triangles files. False for Memfem files which are indexed from 1.
     */
    void GetNodesAtSurface(const std::string& rSurfaceFileName,
                           std::vector<unsigned>& rSurfaceNodes,
                           bool indexFromZero=true) const;

    /**
     *  Helper function for GetNodesAtSurface
     *  @param rLineFromFile  A line in a surface file.
     *  @param rSurfaceNodeIndexSet  The nodes in the element corresponding to this line.
     *  @param offset  is the lowest index of a node in the original mesh (0 for native triangles or 1 for MEMFEM).
     */
    void ProcessLine(const std::string& rLineFromFile,
                     std::set<unsigned>& rSurfaceNodeIndexSet,
                     unsigned offset) const;

    /**
     * Helper method to calculate the distance between the node and the Endocardial surface
     * as defined to be the closest surface to the node out of left ventricle and right ventricle.
     *
     * @param nodeIndex is the index of the node in the mesh
     * @return the distance
     */
    double GetDistanceToEndo(unsigned nodeIndex);

     /**
     * Helper method to calculate the distance between the node and the Epicardial surface
     *
     * @param nodeIndex is the index of the node in the mesh
     * @return the distance
     */
    double GetDistanceToEpi(unsigned nodeIndex);

    /** The mesh of the problem*/
    AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

    /** Vector to store the distance map to epicardium*/
    std::vector<double> mDistMapEpicardium;

    /** Vector to store the distance map to endocardium*/
    std::vector<double> mDistMapEndocardium;

    /** Vector to store the distance map to the right ventricle surface*/
    std::vector<double> mDistMapRightVentricle;

    /** Vector to store the distance map to the left ventricle surface*/
    std::vector<double> mDistMapLeftVentricle;

    /** Flag used to tell the methods whether two or three surfaces have been supplied*/
    unsigned mNumberOfSurfacesProvided;

    /** Vector to store the layer for each node*/
    std::vector<HeartLayerType> mLayerForEachNode;

    /**
     * @return a bounding box for a group of node indices (such as the epi-surface)
     *
     * @param rSurfaceNodes The indices of the nodes which represent this surface
     */
    ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfSurface(const std::vector<unsigned>& rSurfaceNodes);

public:
/// \todo #1703 Perhaps add these constants to HeartConfig...
    /** Left ventricular wall */
    static const HeartRegionType LEFT_VENTRICLE_WALL=1001;
    /** Right ventricular wall */
    static const HeartRegionType RIGHT_VENTRICLE_WALL=1002;
    /** Left portion of the septum */
    static const HeartRegionType LEFT_SEPTUM=1003;
    /** Right portion of the septum */
    static const HeartRegionType RIGHT_SEPTUM=1004;
    /** Endocardial surface of the left ventricle */
    static const HeartRegionType LEFT_VENTRICLE_SURFACE=1005;
    /** Endocardial surface of the right ventricle */
    static const HeartRegionType RIGHT_VENTRICLE_SURFACE=1006;
    /** Unknown node type (should never occur...) */
    static const HeartRegionType UNKNOWN=1007;

    /**
     * Constructor for a two surface mesh.
     * File formats: list of nodes, either one per line or multiple (e.g. nodes in each boundary element on surface).
     *
     * @param rMesh: reference to the mesh
     * @param rEpiFile: a file containing a list of global node indices on the epicardial surface
     * @param rEndoFile: a file containing a list of global node indices on the endocardial surface
     * @param indexFromZero  true for native triangles files. false for Memfem files which are indexed from 1.
     */
    HeartGeometryInformation (AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              const std::string& rEpiFile,
                              const std::string& rEndoFile,
                              bool indexFromZero);


    /**
     * Constructor for a three surface mesh.
     * File formats: list of nodes, either one per line or multiple (e.g. nodes in each boundary element on surface).
     *
     * @param rMesh: reference to the mesh
     * @param rEpiFile: a file containing a list of global node indices on the epicardial surface
     * @param rRVFile: a file containing a list of global node indices on the endocardial right ventricular surface (can be empty string)
     * @param rLVFile: a file containing a list of global node indices on the endocardial left ventricular surface (can be empty string)
     * @param indexFromZero  true for native triangles files. false for Memfem files which are indexed from 1.
     *
     * If either rRVFile or rLVfile are the empty string, then it is assumed that this is a
     * wedge preparation for left or right ventricle, respectively.  That is, the ventricle with a non-empty string.
     * If both are empty strings then throws exception.
     */
    HeartGeometryInformation (AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              const std::string& rEpiFile,
                              const std::string& rLVFile,
                              const std::string& rRVFile,
                              bool indexFromZero);

    /**
     * Alternative constructor that takes in the file containing a list of numbers (as many as the number of nodes).
     * Each number specifies the layer for the corresponding node.
     *
     * This constructor should be called if the heterogeneities have /already/ been computed
     * by an instance of this class and written to file by the WriteLayerForEachNode() method.
     *
     * @param nodeHeterogeneityFileName the file name.
     */
    HeartGeometryInformation (std::string nodeHeterogeneityFileName);

    /**
     * @param nodeIndex index is the index of the node in the mesh
     * @return the region type based on the relative distances to epi and endocardial surfaces
     */
    HeartRegionType GetHeartRegion (unsigned nodeIndex) const;

    /**
     *
     * @return the distance map to the epicardium
     */
    std::vector<double>& rGetDistanceMapEpicardium()
    {
        return mDistMapEpicardium;
    }

    /**
     *
     * @return the distance map to the endocardium
     */
    std::vector<double>& rGetDistanceMapEndocardium()
    {
        assert(mNumberOfSurfacesProvided==2);
        return mDistMapEndocardium;
    }

    /**
     *
     * @return the distance map to the right ventricle
     */
    std::vector<double>& rGetDistanceMapRightVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapRightVentricle;
    }

    /**
     *
     * @return the distance map to the left ventricle
     */
    std::vector<double>& rGetDistanceMapLeftVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapLeftVentricle;
    }

    /** @return the nodes on the epicardial surface */
    const std::vector<unsigned>& rGetNodesOnEpiSurface()
    {
        return mEpiSurface;
    }


    /** @return the nodes on the endocardial surface */
    const std::vector<unsigned>& rGetNodesOnEndoSurface()
    {
        assert(mNumberOfSurfacesProvided==2);
        return mEndoSurface;
    }

    /** @return the nodes on the endocardial left ventricular surface */
    const std::vector<unsigned>& rGetNodesOnLVSurface()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mLVSurface;
    }

    /** @return the nodes on the endocardial right ventricular surface */
    const std::vector<unsigned>& rGetNodesOnRVSurface()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mRVSurface;
    }

    /**
     * @return the layer for every node in the mesh.
     */
    const std::vector<HeartLayerType>& rGetLayerForEachNode()
    {
        assert(mLayerForEachNode.size()>0);
        return mLayerForEachNode;
    }

    /**
     * @return the relative position within the wall thickness (normalised to [0,1])
     * @param nodeIndex index is the index of the node in the mesh
     * @return the relative position
     */
    double CalculateRelativeWallPosition(unsigned nodeIndex);

    /**
     *  Compute which layer (endocardial, midmyocardial or epicardial) each node is in
     *  @param epiFraction is the fraction of wall designed to be epicardial layer
     *  @param endoFraction is the fraction of wall designed to be endocardial layer
     */
    void DetermineLayerForEachNode(double epiFraction, double endoFraction);

    /**
     *  Write the layer for each node. DetermineLayerForEachNode() must have been
     *  called first.
     *
     *  @param outputDir Output directory - note not cleaned
     *  @param file Output file
     */
    void WriteLayerForEachNode(std::string outputDir, std::string file);

    /**
     * @return bounding box
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * epicardial surface
     *
     */
    inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfEpi()
    {
        return CalculateBoundingBoxOfSurface(mEpiSurface);
    }
    /**
     * @return bounding box
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * endocardial surface
     *
     */
    inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfEndo()
    {
        return CalculateBoundingBoxOfSurface(mEndoSurface);
    }
    /**
     * @return bounding box
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * endocardial left ventricular surface
     *
     */
    inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfLV()
    {
        return CalculateBoundingBoxOfSurface(mLVSurface);
    }
    /**
     * @return bounding box
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * endocardial left ventricular surface
     *
     */
     inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfRV()
    {
        return CalculateBoundingBoxOfSurface(mRVSurface);
    }
};
#endif //HEARTGEOMETRYINFORMATION_HPP_

