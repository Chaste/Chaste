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

#ifndef IMMERSEDBOUNDARYMESH_HPP_
#define IMMERSEDBOUNDARYMESH_HPP_

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMeshWriter;

#include <set>
#include <vector>

#include <boost/polygon/voronoi.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ChasteSerialization.hpp"
#include "FluidSource.hpp"
#include "ImmersedBoundaryArray.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "ImmersedBoundaryMeshReader.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "Node.hpp"


/**
 * An immersed boundary mesh class, in which elements may contain different numbers of nodes.
 * This is facilitated by the ImmersedBoundaryElement class.
 *
 * This class allows immersed boundary simulations.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestImmersedBoundaryMesh;

protected:
    /** Number of grid points in x direction */
    unsigned mNumGridPtsX;

    /** Number of grid points in y direction */
    unsigned mNumGridPtsY;

    /** Characteristic node spacing */
    double mCharacteristicNodeSpacing;

    /** The required spacing when an element divides */
    double mElementDivisionSpacing;

    /** The distance at which two elements are neighbouring */
    double mNeighbourDist;

    /**
     * A summary statistic used by UpdateNodeLocationsVoronoiDiagramIfOutOfDate()
     * to determine if mNodeLocationsVoronoiDiagram is out of date
     */
    double mSummaryOfNodeLocations;

    /** The distance above which a cell vertex moving will trigger a step size exception */
    double mCellRearrangementThreshold;

    /**
     * A halo distance around the unit square, used when calculating the node voronoi diagram
     */
    static constexpr double mVoronoiHalo = 0.1;

    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;

    /** 2D grid for fluid x velocity */
    multi_array<double, 3> m2dVelocityGrids;

    /** 3D grid for fluid x velocity */
    multi_array<double, 4> m3dVelocityGrids;

    /** Vector of pointers to ImmersedBoundaryElements. */
    std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** Vector of pointers to lamina ImmersedBoundaryElements. */
    std::vector<ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>*> mLaminas;

    /** Vector of fluid sources related to elements. */
    std::vector<std::shared_ptr<FluidSource<SPACE_DIM>>> mElementFluidSources;

    /** Vector of fluid sources used to balance those of the elements. */
    std::vector<std::shared_ptr<FluidSource<SPACE_DIM>>> mBalancingFluidSources;

    /** A voronoi diagram of the node locations, used to determine the element neighbours */
    boost::polygon::voronoi_diagram<double> mNodeLocationsVoronoiDiagram;

    /**
     * A vector keeping track voronoi cell IDs indexed by node index.  This vector is the length of the largest node
     * index, and is updated by UpdateNodeLocationsVoronoiDiagramIfOutOfDate().
     *
     * The value mVoronoiCellIdsIndexedByNodeIndex[node_idx] gives the ID of the voronoi cell corresponding to the node
     * with index node_idx.
     */
    std::vector<unsigned> mVoronoiCellIdsIndexedByNodeIndex;

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     * @return local index
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     * @return local index
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     * @return local index
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /**
     * Determine whether each element is on the boundary or not, and call the element's SetIsBoundaryElement() method.
     *
     * It is not possible to define in precise terms whether an element is on the boundary (as there is no commonality
     * of nodes as in a vertex population).  Instead we use information from the node locations voronoi diagram.  Any
     * nodes with infinite voronoi edges are certainly in boundary elements, as well as nodes in elements whose voronoi
     * cells are "too big".
     */
    void TagBoundaryElements();

    /**
     * Divide an element along the axis passing through two of its nodes.
     *
     * @param pElement the element to divide
     * @param nodeAIndex the local index of one node within this element
     * @param nodeBIndex the local index of another node within this element
     * @param centroid the centroid of the element being divided
     * @param axisOfDivision the specified division axis
     *
     * @return the index of the new element
     */
    unsigned DivideElement(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                           unsigned nodeAIndex,
                           unsigned nodeBIndex,
                           c_vector<double, SPACE_DIM> centroid,
                           c_vector<double, SPACE_DIM> axisOfDivision);

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the ImmersedBoundaryMesh and its member variables. Note that this will
     * write out an ImmersedBoundaryMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void save(Archive& archive, const unsigned int version) const
    {
        archive& boost::serialization::base_object<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                                       ArchiveLocationInfo::GetMeshFilename(),
                                                                       false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Load a mesh by using ImmersedBoundaryMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void load(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    /** Forward declaration of element iterator. */
    class ImmersedBoundaryElementIterator;

    /** Forward declaration of lamina iterator. */
    class ImmersedBoundaryLaminaIterator;

    /**
     * @return a reference to the const collection of nodes
    */
    const std::vector<Node<SPACE_DIM>*>& rGetNodes() const;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted elements
     */
    inline ImmersedBoundaryElementIterator GetElementIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline ImmersedBoundaryElementIterator GetElementIteratorEnd();

    /**
     * @return an iterator to the first lamina in the mesh.
     *
     * @param skipDeletedLaminas whether to include deleted laminas
     */
    inline ImmersedBoundaryLaminaIterator GetLaminaIteratorBegin(bool skipDeletedLaminas = true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline ImmersedBoundaryLaminaIterator GetLaminaIteratorEnd();

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param elements vector of pointers to ImmersedBoundaryElements
     * @param laminas vector of pointers to lamina ImmersedBoundaryElements
     * @param numGridPtsX the number of grid points in the x direction
     * @param numGridPtsY the number of grid points in the y direction
     */
    ImmersedBoundaryMesh(std::vector<Node<SPACE_DIM>*> nodes,
                         std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*> elements,
                         std::vector<ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>*> laminas = {},
                         unsigned numGridPtsX = 128u,
                         unsigned numGridPtsY = 128u);

    /**
     * Default constructor for use by serializer.
     */
    ImmersedBoundaryMesh();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of ImmersedBoundaryElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of ImmersedBoundaryElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @return the number of laminas in the mesh.
     */
    unsigned GetNumLaminas() const;

    /**
     * @return the number of fluid mesh points in the x direction.
     */
    unsigned GetNumGridPtsX() const;

    /**
     * @return the number of fluid mesh points in the y direction.
     */
    unsigned GetNumGridPtsY() const;

    /**
     * @return the characteristic node spacing.
     */
    double GetCharacteristicNodeSpacing() const;

    /**
     * @return the spacing ratio: Characteristic Node Spacing / Fluid Grid Spacing
     */
    double GetSpacingRatio() const;

    /**
     * @return the maximum node index (nodes may not be indexed consecutively), or UINT_MAX if there are none
     */
    unsigned GetMaxNodeIndex() const;

    /**
     * @return the maximum element index (elements may not be indexed consecutively), or UINT_MAX if there are none
     */
    unsigned GetMaxElementIndex() const;

    /**
     * @return the maximum lamina index (laminas may not be indexed consecutively), or UINT_MAX if there are none
     */
    unsigned GetMaxLaminaIndex() const;

    /**
     * Overridden GetVectorFromAtoB() method.
     *
     * Evaluates the (surface) distance between two points in a periodic
     * geometry.
     *
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     * @return the vector from location1 to location2
     */
    c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2);

    /**
     * Move the node with a particular index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param point the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * Transform a location so that it lies within the domain [0,1)x[0,1).  This method assumes the location is already
     * 'close' to the domain, for instance because a location has previously been 'straightened out' for geometric
     * calculations that require no looping due to periodicity.
     *
     * @param rLocation the location to transform
     */
    void ConformToGeometry(c_vector<double, SPACE_DIM>& rLocation);

    /**
     * @return reference to non-modifiable 2d fluid velocity grids.
     */
    const multi_array<double, 3>& rGet2dVelocityGrids() const;

    /**
     * @return reference to non-modifiable 3d fluid velocity grids.
     */
    //const multi_array<double, 4>& rGet3dVelocityGrids() const;

    /**
     * @return reference to modifiable 2d fluid velocity grids.
     */
    multi_array<double, 3>& rGetModifiable2dVelocityGrids();

    /**
     * @param meshPointsX the new number of fluid mesh points in the x direction.
     */
    void SetNumGridPtsX(unsigned meshPointsX);

    /**
     * @param meshPointsY the new number of fluid mesh points in the x direction.
     */
    void SetNumGridPtsY(unsigned meshPointsY);

    /**
     * @param numGridPts the new number of fluid mesh points in both directions.
     */
    void SetNumGridPtsXAndY(unsigned numGridPts);

    /**
     * @param nodeSpacing the new characteristic node spacing.
     */
    void SetCharacteristicNodeSpacing(double nodeSpacing);

    /**
     * Add a node to the mesh
     * @param pNewNode a pointer to the node to be added to the mesh.
     * @return the index of the new node
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * @return reference to vector of element-associated fluid sources
     */
    std::vector<std::shared_ptr<FluidSource<SPACE_DIM>>>& rGetElementFluidSources();

    /**
     * @return reference to vector of balancing fluid sources
     */
    std::vector<std::shared_ptr<FluidSource<SPACE_DIM>>>& rGetBalancingFluidSources();

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * @param index  the global index of a specified immersed boundary element.
     *
     * @return a pointer to the immersed boundary element
     */
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * @param index  the global index of a specified immersed boundary lamina.
     *
     * @return a pointer to the immersed boundary lamina
     */
    ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* GetLamina(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * A formula for the centroid of a plane polygon may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return (centroid_x, centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();

    /**
     * Get the volume (or area in 2D, or length in 1D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return the surface area of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of the Voronoi cell of an element.
     * The voronoi cell of an element is the union of the voronoi cells of all nodes
     * in the element.
     *
     * @param elemIdx the global index of a specified immersed boundary element
     * @return the voronoi surface area of the element
     */
    double GetVoronoiSurfaceAreaOfElement(unsigned elemIdx);

    /**
     * Compute the average node spacing of an element.
     *
     * @param index  the global index of a specified immersed boundary element
     * @param recalculate whether or not to recalculate the value
     * @return the average node spacing of the element
     */
    double GetAverageNodeSpacingOfElement(unsigned index, bool recalculate = true);

    /**
     * Compute the average node spacing of a lamina.
     *
     * @param index  the global index of a specified immersed boundary lamina
     * @param recalculate whether or not to recalculate the value
     * @return the average node spacing of the lamina
     */
    double GetAverageNodeSpacingOfLamina(unsigned index, bool recalculate = true);

    /**
     * Compute the second moments and product moment of area for a given 2D element
     * about its centroid. These are:
     *
     * I_xx, the second moment of area about an axis through the centroid of the
     * element parallel to the x-axis;
     *
     * I_yy, the second moment of area about an axis through the centroid of the
     * element parallel to the y-axis;
     *
     * and I_xy, product moment of area through the centroid of the element.
     *
     * Formulae for these quantities may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This method is used within GetShortAxisOfElement() to compute the direction
     * of the shortest principal axis passing through the centroid, or 'short axis',
     * of the element.
     *
     * Note that by definition, the second moments of area must be non-negative,
     * while the product moment of area may not be.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

    /**
     * Get the elongation shape factor of a given element.
     * This is defined as the square root of the ratio of
     * the two second moments of the element around its
     * principal axes.
     *
     * @param elementIndex index of an element in the mesh
     *
     * @return the elongation shape factor of the element.
     */
    double GetElongationShapeFactorOfElement(unsigned elementIndex);

    /**
     * Compute tortuosity, defined as the ratio of total length to straight-line length, of piecewise linear curve
     * through centroids of successive elements.
     *
     * @return the tortuosity of the mesh
     */
    double GetTortuosityOfMesh();

    /**
     * Calculate the skewness of the mass-distribution of a given element along a line perpendicular to a given axis,
     * as a measure of element asymmetry.
     *
     * This method is only valid in 2D.
     *
     * @param elemIndex the index of the element in this mesh
     * @param axis the axis perpendicular to which the distribution is found
     * @return the skewness of the mass distribution for the given element about the given axis
     */
    double GetSkewnessOfElementMassDistributionAboutAxis(unsigned elemIndex, c_vector<double, SPACE_DIM> axis);

    /**
     * Calculate the bounding box of an element specified by its index.
     *
     * @param index the index of the element
     * @return a Chaste cuboid representing the element bounding box
     */
    ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfElement(unsigned index);

    /**
     * Compute the direction of the shortest principal axis passing through the centroid,
     * or 'short axis', of a given element. This is the eigenvector associated with the
     * eigenvalue of largest magnitude of the inertia matrix
     *
     * J = (  I_xx  -I_xy )
     *     ( -I_xy   I_yy )
     *
     * whose entries are computed by calling the method CalculateMomentsOfElement().
     *
     * Note that if the nodes owned by the element are supplied in clockwise rather than
     * anticlockwise manner, or if this arises when any periodicity is enforced, then the
     * sign of each moment may be incorrect change. This means that we need to consider the eigenvalue
     * of largest magnitude rather than largest value when computing the short axis of the
     * element.
     *
     * If the element is a regular polygon then the eigenvalues of the inertia tensor are
     * equal: in this case we return a random unit vector.
     *
     * This method is only implemented in 2D at present.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return a unit vector giving the direction of the short axis
     */
    c_vector<double, SPACE_DIM> GetShortAxisOfElement(unsigned index);

    /**
     * Divide an element along a specified axis.
     *
     * If the new nodes (intersections of axis with element) are within
     * mCellRearrangementThreshold of existing nodes then they are
     * moved 2*mCellRearrangementThreshold away.
     *
     * @param pElement the element to divide
     * @param axisOfDivision axis to divide the element by
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongGivenAxis(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                         c_vector<double, SPACE_DIM> axisOfDivision,
                                         bool placeOriginalElementBelow = false);

    /**
     * Divide an element along its short axis.
     *
     * @param pElement the element to divide
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongShortAxis(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                         bool placeOriginalElementBelow = false);

    /**
     * @return mElementDivisionSpacing
     */
    double GetElementDivisionSpacing();

    /**
     * @param elementDivisionSpacing the new value of mElementDivisionSpacing
     */
    void SetElementDivisionSpacing(double elementDivisionSpacing);

    /**
     * @return mNeighbourDist
     */
    double GetNeighbourDist() const;

    /**
     * @param cellRearrangementThreshold the new value of mCellRearrangementThreshold
     */
    void SetCellRearrangementThreshold(double cellRearrangementThreshold);

    /**
     * @return the maximum distance a cell vertex can move without triggering a step size exception
     */
    double GetCellRearrangementThreshold();

    /**
     * @param neighbourDist the new value of mNeighbourDist
     */
    void SetNeighbourDist(double neighbourDist);

    /**
     * Update mNodeLocationsVoronoiDiagram if it is out of date, which is the case when mSummaryOfNodeLocations is
     * out-of-date, i.e. nodes have moved location.  This ensures the update is performed at most once per time step.
     */
    void UpdateNodeLocationsVoronoiDiagramIfOutOfDate();

    /**
     * ReMesh method that evenly redistributes nodes around each element and lamina.
     *
     * @param randomOrder whether to remesh elements and laminas starting from a random location. Defaults to false.
     */
    void ReMesh(bool randomOrder=false);

    /**
     * ReMesh method that evenly redistributes nodes around a specific element.
     *
     * @param pElement the element to remesh
     * @param randomOrder whether to remesh elements and laminas starting from a random location
     */
    void ReMeshElement(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* pElement, bool randomOrder);

    /**
     * ReMesh method that evenly redistributes nodes along a specific lamina.
     *
     * @param pLamina the lamina to remesh
     * @param randomOrder whether to remesh elements and laminas starting from a random location
     */
    void ReMeshLamina(ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* pLamina, bool randomOrder);

    /**
     * Determine whether two nodes belong to the same element, or the same lamina
     *
     * @param pNodeA a pointer to a node
     * @param pNodeB a pointer to a different node
     * @return whether the two nodes belong to the same element or lamina
     */
    bool NodesInDifferentElementOrLamina(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Given an element index, find a set containing the indices of its neighbouring elements.
     *
     * @param elemIdx global index of the element
     * @return a set  of its neighbouring element indices
     */
    std::set<unsigned> GetNeighbouringElementIndices(unsigned elemIdx);

    /**
     * Get the length of an edge in the voronoi diagram mNodeLocationsVoronoiDiagram. This gives the "real" length after
     * undoing the scaling required when calculating the diagram.
     * @param rEdge the edge to calculate the length of
     * @return the length of rEdge
     */
    double CalculateLengthOfVoronoiEdge(const boost::polygon::voronoi_diagram<double>::edge_type& rEdge);

    /**
     * Calculate the polygon distribution for the mesh: number of {0, 1, 2, 3, 4, 5,..., 12+}-gons.
     * Note that the vector will always begin {0, 0, 0, ...} as there can be no 0, 1, or 2-gons, but this choice means
     * that accessing the nth element of the vector gives you the number of n-gons which seems to be most natural.
     * All 12-sided and higher order polygons are accumulated in the array[12] position.
     *
     * @return an array of length 13 representing the polygon distribution.
     */
    std::array<unsigned, 13> GetPolygonDistribution();

    /**
     * Get the voronoi diagram of node locations. This may be needed for population writers and others.
     * @param update whether to update the diagram before returning the reference (default true)
     * @return mNodeLocationsVoronoiDiagram
     */
    const boost::polygon::voronoi_diagram<double>& rGetNodeLocationsVoronoiDiagram(bool update=true);

    /** @return mVoronoiCellIdsIndexedByNodeIndex */
    const std::vector<unsigned int>& GetVoronoiCellIdsIndexedByNodeIndex() const;

    /**
     * Helper method for voronoi functions.  Scale a location up to the integer grid needed by boost voronoi.
     *
     * @param location a location in [-mVoronoiHalo, 1.0 + mVoronoiHalo]
     * @return a corresponding integer location scaled to [INT_MIN, INT_MAX]
     */
    int ScaleUpToVoronoiCoordinate(double location) const;

    /**
     * Helper method for voronoi functions.  Scale a distance down from the voronoi coordinates.
     *
     * @param distance a distance in the voronoi diagram
     * @return a distance in the immersed boundary domain
     */
    double ScaleDistanceDownFromVoronoi(const double distance) const;


    /**
     * A smart iterator over the elements in the mesh.
     */
    class ImmersedBoundaryElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline ImmersedBoundaryElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * ImmersedBoundaryMesh::GetElementIteratorBegin and ImmersedBoundaryMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        ImmersedBoundaryElementIterator(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                        typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*>::iterator elementIter,
                                        bool skipDeletedElements = true);

    private:
        /** The mesh we're iterating over. */
        ImmersedBoundaryMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedElement();
    };

    /**
     * A smart iterator over the laminas in the mesh.
     */
    class ImmersedBoundaryLaminaIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current lamina.
         * @return reference
         * Make sure to use a reference for the result to avoid copying laminas unnecessarily.
         */
        inline ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline ImmersedBoundaryLaminaIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * ImmersedBoundaryMesh::GetlaminaIteratorBegin and ImmersedBoundaryMesh::GetlaminaIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param laminaIter where to start iterating
         * @param skipDeletedLaminas whether to include deleted laminas
         */
        ImmersedBoundaryLaminaIterator(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                       typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>*>::iterator laminaIter,
                                       bool skipDeletedLaminas = true);

    private:
        /** The mesh we're iterating over. */
        ImmersedBoundaryMesh& mrMesh;

        /** The actual lamina iterator. */
        typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>*>::iterator mLaminaIter;

        /** Whether to skip deleted laminas. */
        bool mSkipDeletedLaminas;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this lamina.
         * @return true if allowed
         */
        inline bool IsAllowedLamina();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ImmersedBoundaryMesh)

////////////////////////////////////////////////////////////////////////////////////////
// ImmersedBoundaryElementIterator class implementation - most methods are inlined    //
////////////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
    bool skipDeletedElements)
{
    return ImmersedBoundaryElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return ImmersedBoundaryElementIterator(*this, mElements.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator!=(const typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    } while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::ImmersedBoundaryElementIterator(
    ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
    typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*>::iterator elementIter,
    bool skipDeletedElements)
        : mrMesh(rMesh),
          mElementIter(elementIter),
          mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.empty())
    {
        // Cope with empty meshes
        mElementIter = mrMesh.mElements.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mElementIter == mrMesh.mElements.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

///////////////////////////////////////////////////////////////////////////////////////
// ImmersedBoundaryLaminaIterator class implementation - most methods are inlined    //
///////////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetLaminaIteratorBegin(
    bool skipDeletedLaminas)
{
    return ImmersedBoundaryLaminaIterator(*this, mLaminas.begin(), skipDeletedLaminas);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetLaminaIteratorEnd()
{
    return ImmersedBoundaryLaminaIterator(*this, mLaminas.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::operator*()
{
    assert(!IsAtEnd());
    return **mLaminaIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::operator->()
{
    assert(!IsAtEnd());
    return *mLaminaIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::operator!=(const typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator& rOther)
{
    return mLaminaIter != rOther.mLaminaIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::operator++()
{
    do
    {
        ++mLaminaIter;
    } while (!IsAtEnd() && !IsAllowedLamina());

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::ImmersedBoundaryLaminaIterator(
    ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
    typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>*>::iterator laminaIter,
    bool skipDeletedLaminas)
        : mrMesh(rMesh),
          mLaminaIter(laminaIter),
          mSkipDeletedLaminas(skipDeletedLaminas)
{
    if (mrMesh.mLaminas.empty())
    {
        // Cope with empty meshes
        mLaminaIter = mrMesh.mLaminas.end();
    }
    else
    {
        // Make sure we start at an allowed lamina
        if (mLaminaIter == mrMesh.mLaminas.begin() && !IsAllowedLamina())
        {
            ++(*this);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::IsAtEnd()
{
    return mLaminaIter == mrMesh.mLaminas.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryLaminaIterator::IsAllowedLamina()
{
    return !(mSkipDeletedLaminas && (*this)->IsDeleted());
}

#endif /*IMMERSEDBOUNDARYMESH_HPP_*/
