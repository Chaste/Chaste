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

#ifndef FINECOARSEMESHPAIR_HPP_
#define FINECOARSEMESHPAIR_HPP_

#include "AbstractTetrahedralMesh.hpp"
#include "DistributedBoxCollection.hpp"
#include "QuadraturePointsGroup.hpp"
#include "GaussianQuadratureRule.hpp"
#include "Warnings.hpp"
#include "CommandLineArguments.hpp"

/**
 * At the beginning of a two mesh simulation we need to figure out and store
 * which fine-mesh element each (coarse-mesh) quadrature point is in, and
 * what the weight of that Gauss point for that particular element is. This struct
 * just contains this two pieces of data
 */
template<unsigned DIM>
struct ElementAndWeights
{
    unsigned ElementNum; /**< Which element*/
    c_vector<double, DIM+1> Weights; /**<Gauss weights for this element*/
};

/**
 * Class for a pair of meshes, one fine, one coarse, which should cover the same domain (or very nearly match).
 * This class is used to set up interpolation information from one mesh to the other.
 *
 * The functionality is based on the four information-transfers required in
 * cardiac electro-mechanics problems
 *
 * -# Calcium (or voltage) to induce deformation:
 *          FINE(electrics) MESH NODEs  --->  COARSE(mechanics) MESH QUADRATURE POINTS
 * -#  Deformation gradient (assume constant in any coarse element) for altering conductivities:
 *          COARSE ELEMENTS  --->  FINE ELEMENTS
 * -# Deformation gradient/fibre-stretch (assume constant in any coarse element) for cell-model
 *       stretch activated channels
 *           COARSE ELEMENTS  --->  FINE NODES
 * -#  Voltage visualisation on coarse mesh
 *          FINE NODES ---> COARSE NODES
 *
 * The usage of this class for each of these tasks is:
 *
 * -# FINE NODEs  --->  COARSE QUADRATURE POINTS
 *          FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
 *          mesh_pair.SetUpBoxesOnFineMesh();
 *          mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, false);
 *          mesh_pair.rGetElementsAndWeights();
 * -# COARSE ELEMENTS  --->  FINE ELEMENTS
 *          FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
 *          mesh_pair.SetUpBoxesOnCoarseMesh();
 *          mesh_pair.ComputeCoarseElementsForFineElementCentroids();
 *          mesh_pair.rGetCoarseElementsForFineElementCentroids();
 * -# COARSE ELEMENTS  --->  FINE NODES
 *          FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
 *          mesh_pair.SetUpBoxesOnCoarseMesh();
 *          mesh_pair.ComputeCoarseElementsForFineNodes();
 *          mesh_pair.rGetCoarseElementsForFineNodes();
 * -#  FINE NODES ---> COARSE NODES
 * Note the this should not be done at the same time as (1), because the results are stored in the same place
 *          FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
 *          mesh_pair.SetUpBoxesOnFineMesh();
 *          mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(false);
 *          mesh_pair.rGetElementsAndWeights();
 *
 *
 * To see progression for any of these methods, run test from the command line with '-mesh_pair_verbose' as
 * a command line parameter
 *
 */
template <unsigned DIM>
class FineCoarseMeshPair
{
friend class TestFineCoarseMeshPair;

private:

    /** Fine mesh. */
    AbstractTetrahedralMesh<DIM,DIM>& mrFineMesh;

    /** Coarse mesh (often this will be a quadratic mesh). */
    AbstractTetrahedralMesh<DIM,DIM>& mrCoarseMesh;

    /**
     * Boxes on the fine mesh domain, for easier determination of
     * containing element for a given point. */
    DistributedBoxCollection<DIM>* mpFineMeshBoxCollection;

    /**
     * Boxes on the coarse mesh domain, for easier determination of
     * containing element for a given point.
     */
    DistributedBoxCollection<DIM>* mpCoarseMeshBoxCollection;

    /**
     * The containing elements and corresponding weights in the fine
     * mesh for the set of points given. The points may have been
     * quadrature points in the coarse mesh, or nodes in coarse mesh, etc.
     */
    std::vector<ElementAndWeights<DIM> > mFineMeshElementsAndWeights;

    /** Indices of the points which were found to be outside the fine mesh. */
    std::vector<unsigned> mNotInMesh;

    /**
     * The corresponding weights, for the nearest elements, of the
     * points which were found to be outside the fine mesh.
     */
    std::vector<c_vector<double,DIM+1> > mNotInMeshNearestElementWeights;

    /**
     * 2 values,
     *   [0] number of points for which the containing element was found
     *   [1] num points not in the searched mesh at all
     *
     *  Note, after ComputeFineElementsAndWeightsForCoarseQuadPoints() or
     *  ComputeFineElementsAndWeightsForCoarseNodes(), then
     *   mStatisticsCounters[1] = mNotInMesh.size() = mNotInMeshNearestElementWeights.size();
     */
    std::vector<unsigned> mStatisticsCounters;

    /**
     * The element in the coarse mesh that each fine mesh node is contained
     * in (or nearest to). ComputeCoarseElementsForFineNodes() needs to be
     * called for this to be set up.
     */
    std::vector<unsigned> mCoarseElementsForFineNodes;

    /**
     * The element in the coarse mesh that each fine element centroid is contained
     * in (or nearest to). ComputeCoarseElementsForFineElementCentroids() needs to
     * be called for this to be set up.
     */
    std::vector<unsigned> mCoarseElementsForFineElementCentroids;

    /**
     * For a given point, compute the containing element and corresponding weight
     * in the fine mesh.
     *
     * @param rPoint The point
     * @param safeMode See documentation for ComputeFineElementsAndWeightsForCoarseQuadPoints()
     * @param boxForThisPoint The box in the fine box collection containing this point
     * @param index The index into the mFineMeshElementsAndWeights std::vector
     */
    void ComputeFineElementAndWeightForGivenPoint(ChastePoint<DIM>& rPoint,
                                                  bool safeMode,
                                                  unsigned boxForThisPoint,
                                                  unsigned index);

    /**
     * @return for a given point the containing element in the coarse mesh
     *
     * @param rPoint The point
     * @param safeMode See documentation in ComputeCoarseElementsForFineNodes
     * @param boxForThisPoint The box in coarse box collection containing this point
     */
    unsigned ComputeCoarseElementForGivenPoint(ChastePoint<DIM>& rPoint,
                                               bool safeMode,
                                               unsigned boxForThisPoint);

    /**
     * Set up a box collection on the given mesh. Should only be called using either
     *   SetUpBoxes(*mpFineMesh, boxWidth, mpFineBoxCollection)  (from SetUpBoxesOnFineMesh)
     * or
     *   SetUpBoxes(*mpCoarseMesh, boxWidth, mpCoarseBoxCollection)  (from SetUpBoxesOnCoarseMesh)
     *
     * @param rMesh The mesh, either *mpFineMesh or *mpCoarseMesh)
     * @param boxWidth box width (see SetUpBoxesOnCoarseMesh() dox)
     * @param rpBoxCollection reference to either mpFineBoxCollection or mpCoarseBoxCollection
     */
    void SetUpBoxes(AbstractTetrahedralMesh<DIM,DIM>& rMesh,
                    double boxWidth,
                    DistributedBoxCollection<DIM>*& rpBoxCollection);

    /**
     * Helper method. Gets all the elements in the given box, in the given box collection
     * and puts them in the returned std::vector.
     *
     * @param rpBoxCollection Reference to the box collection to use (either mpFineBoxCollection
     *  or mpCoarseBoxCollection)
     * @param boxIndex box index
     * @param rElementIndices The returned vector of element indices in that box of the box collection. Not
     *  cleared before use.
     */
    void CollectElementsInContainingBox(DistributedBoxCollection<DIM>*& rpBoxCollection,
                                        unsigned boxIndex,
                                        std::set<unsigned>& rElementIndices);
    /**
     * Helper method. Gets all the elements in the given box, or in a box local to the given box,
     * in the given box collection, and puts them in the returned std::vector.
     *
     * @param rpBoxCollection Reference to the box collection to use (either mpFineBoxCollection
     *  or mpCoarseBoxCollection)
     * @param boxIndex box index
     * @param rElementIndices The returned vector of element indices in that box or a local box. Not
     *  cleared before use.
     */
    void CollectElementsInLocalBoxes(DistributedBoxCollection<DIM>*& rpBoxCollection,
                                     unsigned boxIndex,
                                     std::set<unsigned>& rElementIndices);

    /**
     * Resets mNotInMesh, mNotInMeshNearestElementWeights and
     * mStatisticsCounters.
     */
    void ResetStatisticsVariables();
    /**
     * In parallel: share mStatisticsCounters and all the "fine element weights for..." information
     */
    void ShareFineElementData();
    /**
     * In parallel: share mStatisticsCounters and all the "this coarse element index for..." information
     */
    void ShareCoarseElementData();

public:

    /**
     * Constructor sets up domain size.
     *
     * @param rFineMesh Fine mesh (reference)
     * @param rCoarseMesh Coarse mesh (reference)
     */
    FineCoarseMeshPair(AbstractTetrahedralMesh<DIM,DIM>& rFineMesh, AbstractTetrahedralMesh<DIM,DIM>& rCoarseMesh);

    /**
     * Destructor just deletes the box collection.
     */
    ~FineCoarseMeshPair();

    /**
     * Set up boxes on fine mesh. The elements contained in each box is stored, which makes
     * finding the containing element for a given point much faster.
     * This should be called before ComputeFineElementsAndWeightsForCoarseQuadPoints() or
     * ComputeFineElementsAndWeightsForCoarseNodes().
     *
     * @param boxWidth width to use for the boxes (which will be cubes). Note that a domain
     *    which is a touch larger than the smallest containing cuboid of the fine mesh is used.
     *    boxWidth defaults to a negative value, in which case a box width such that there are
     *    approximately 20 boxes in the x-direction, unless this width is less than maximum (fine
     *    mesh edge length), in which case it is chosen accordingly.
     */
    void SetUpBoxesOnFineMesh(double boxWidth = -1);

    /**
     * Set up boxes on coarse mesh. The elements contained in each box is stored, which makes
     * finding the containing element for a given point much faster.
     * This should be called before ComputeCoarseElementsForFineNodes() or
     * ComputeCoarseElementsForFineElementCentroids()
     *
     * @param boxWidth width to use for the boxes (which will be cubes). Note that a domain
     *    which is a touch larger than the smallest containing cuboid of the fine mesh is used.
     *    boxWidth defaults to a negative value, in which case a box width such that there are
     *    approximately 20 boxes in the x-direction, unless this width is less than maximum (fine
     *    mesh edge length), in which case it is chosen accordingly.
     */
    void SetUpBoxesOnCoarseMesh(double boxWidth = -1);

    /**
     * Set up the containing (fine) elements and corresponding weights for all the
     * quadrature points in the coarse mesh. Call GetElementsAndWeights() after calling this
     * with the index of the quad point (=the index of the quad point in a QuadraturePointsGroup=
     * the index if the quad points were listed by looping over all the element and then
     * looping over all the quad points).
     *
     * If calling this DO NOT call ComputeFineElementsAndWeightsForCoarseNodes
     * until you do done with this data
     *
     * @param rQuadRule The quadrature rule, used to determine the number of quadrature points per element.
     * @param safeMode This method uses the elements in the boxes to guess which element a quad point is in. If a
     *   quad point is in none of these elements, then if safeMode==true, it will then search the whole mesh.
     *   If safeMode==false it will assume immediately the quad point isn't in the mesh at all. safeMode=false is
     *   will far more efficient with big meshes. It should be fine to use safeMode=false if SetUpBoxesOnFineMesh() is
     *   called with default values.
     */
    void ComputeFineElementsAndWeightsForCoarseQuadPoints(GaussianQuadratureRule<DIM>& rQuadRule,
                                                          bool safeMode);

    /**
     * Set up the containing (fine) elements and corresponding weights for all the
     * nodes in the coarse mesh. Call GetElementsAndWeights() after calling this
     * with the index of the nodes.
     *
     * If calling this DO NOT call ComputeFineElementsAndWeightsForCoarseQuadPoints
     * until you do done with this data.
     *
     * @param safeMode This method uses the elements in the boxes to guess which element a point is in. If a
     *   point is in none of these elements, then if safeMode==true, it will then search the whole mesh.
     *   If safeMode==false it will assume immediately the point isn't in the coarse mesh at all. safeMode=false is
     *   will far more efficient with big meshes. It should be fine to use safeMode=false if SetUpBoxesOnFineMesh() is
     *   called with default values.
     */
    void ComputeFineElementsAndWeightsForCoarseNodes(bool safeMode);

    /**
     * Print information about the number of points found in the
     * searched mesh. What the points are and which mesh was searched
     * depends on whichever of the main Compute methods was last called.
     *
     * If ComputeFineElementsAndWeightsForCoarseQuadPoints() or
     * ComputeFineElementsAndWeightsForCoarseNodes() were last called,
     * the indices of the points that were not found to be contained in
     * the searched mesh are also printed, along with the weights for
     * that point in the nearest element (which indicates how far from
     * the mesh the points are).
     */
    void PrintStatistics();

    /**
     * Compute the element in the coarse mesh that each fine mesh node is contained in (or nearest to).
     * Call SetUpBoxesOnCoarseMesh() before, and rGetCoarseElementsForFineNodes() afterwards.
     *
     * @param safeMode This method uses the elements in the boxes to guess which element a point is in. If a
     *   point is in none of these elements, then if safeMode==true, it will then search the whole mesh.
     *   If safeMode==false it will assume immediately the point isn't in the coarse mesh at all. safeMode=false is
     *   will far more efficient with big meshes. It should be fine to use safeMode=false if SetUpBoxesOnFineMesh() is
     *   called with default values.
     */
    void ComputeCoarseElementsForFineNodes(bool safeMode);

    /**
     * Compute the element in the coarse mesh that each fine element centroid is contained in
     * (or nearest to). Call SetUpBoxesOnCoarseMesh() before, and
     * rGetCoarseElementsForFineElementCentroids() afterwards.
     *
     * @param safeMode This method uses the elements in the boxes to guess which element a point is in. If a
     *   point is in none of these elements, then if safeMode==true, it will then search the whole mesh.
     *   If safeMode==false it will assume immediately the point isn't in the coarse mesh at all. safeMode=false is
     *   will far more efficient with big meshes. It should be fine to use safeMode=false if SetUpBoxesOnFineMesh() is
     *   called with default values.
     */
    void ComputeCoarseElementsForFineElementCentroids(bool safeMode);

    /**
     * @return  A reference to the elements/weights information
     */
    std::vector<ElementAndWeights<DIM> >& rGetElementsAndWeights()
    {
        return mFineMeshElementsAndWeights;
    }

    /**
     * @return the elements in the coarse mesh that each fine mesh node is contained in (or nearest to).
     * ComputeCoarseElementsForFineNodes() needs to be called before calling this.
     */
    std::vector<unsigned>& rGetCoarseElementsForFineNodes()
    {
        assert(mCoarseElementsForFineNodes.size() > 0);
        return mCoarseElementsForFineNodes;
    }

    /**
     * @return the elements in the coarse mesh that each fine mesh element centroid is contained in (or nearest to).
     * ComputeCoarseElementsForFineElementCentroids() needs to be called before calling this.
     */
    std::vector<unsigned>& rGetCoarseElementsForFineElementCentroids()
    {
        assert(mCoarseElementsForFineElementCentroids.size() > 0);
        return mCoarseElementsForFineElementCentroids;
    }

    /**
     * Destroy the box collection for the fine mesh - can be used to free memory once
     * ComputeFineElementsAndWeightsForCoarseQuadPoints (etc) has been called.
     */
    void DeleteFineBoxCollection();

    /**
     * Destroy the box collection for the coarse mesh - can be used to free memory once
     * ComputeCoarseElementsForFineNodes (etc) has been called.
     */
    void DeleteCoarseBoxCollection();

    /**
     * Access the fine mesh of this mesh pair
     * @return the fine mesh
     */
    const AbstractTetrahedralMesh<DIM, DIM>& GetFineMesh() const;

    /**
     * Access the coarse mesh of this mesh pair
     * @return the coarse mesh
     */
    const AbstractTetrahedralMesh<DIM, DIM>& GetCoarseMesh() const;
};

#endif /*FINECOARSEMESHPAIR_HPP_*/
