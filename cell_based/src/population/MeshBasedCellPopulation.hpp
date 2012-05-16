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

#ifndef MESHBASEDCELLPOPULATION_HPP_
#define MESHBASEDCELLPOPULATION_HPP_

#include <map>
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MutableMesh.hpp"
#include "VertexMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>


/**
 * A facade class encapsulating a mesh-based 'cell population'
 *
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh.
 *
 */
template<unsigned DIM>
class MeshBasedCellPopulation : public AbstractCentreBasedCellPopulation<DIM>
{
    friend class TestMeshBasedCellPopulation;
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedCellPopulation<DIM> >(*this);

        /*
         * In its current form the code does not allow the direct serialization
         * of the VertexMesh class, so instead we delete mpVoronoiTessellation.
         */
        delete mpVoronoiTessellation;
        mpVoronoiTessellation = NULL;

        archive & mMarkedSprings;
        archive & mSpringRestLengths;
        archive & mUseAreaBasedDampingConstant;
        archive & mAreaBasedDampingConstantParameter;
        archive & mOutputVoronoiData;
        archive & mOutputCellPopulationVolumes;
        archive & mWriteVtkAsPoints;
        archive & mHasVariableRestLength;


        this->Validate();
    }

protected:
#define COVERAGE_IGNORE // Avoid prototypes being treated as code by gcov

    /**
     * Pointer to a VertexMesh object that stores the Voronoi tessellation that is dual to
     * mrMesh. The tessellation is created by calling CreateVoronoiTessellation() and can
     * be accessed by calling GetVoronoiTessellation().
     *
     * The tessellation can be used to compute the area and perimeter (in 2D) or volume and
     * surface area (in 3D) of the Voronoi element corresponding to each node in the Delaunay
     * mesh (including ghost nodes) by calling the methods GetVolumeOfVoronoiElement() and
     * GetSurfaceAreaOfVoronoiElement() respectively. Each of these methods should be called
     * rather than the relevant method on the VertexMesh. This is because the index of a given
     * Node in mrMesh may not equal the index of the corresponding VertexElement in
     * mpVoronoiTessellation; a map between these indices may be accessed by calling the methods
     * GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex()
     * and GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex() on mpVoronoiTessellation.
     */
    VertexMesh<DIM, DIM>* mpVoronoiTessellation;

    /** Static cast of the mesh from AbstractCellPopulation */
    MutableMesh<DIM, DIM>* mpMutableMesh;

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this cell population has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * Special springs that we want to keep track of for some reason.
     * Currently used to track cells in the process of dividing
     * (which are represented as two cells joined by a shorter spring).
     */
    std::set<std::pair<CellPtr,CellPtr> > mMarkedSprings;

    /**
     * Keeps track of the rest lengths of springs if these are being used in the simulation.
     */
    std::map<std::pair<unsigned,unsigned>, double> mSpringRestLengths;

    /** Results file for elements. */
    out_stream mpVizElementsFile;

    /** Results file for Voronoi data. */
    out_stream mpVoronoiFile;

    /** Results file for cell population volume (in 3D) or area (in 2D) data. */
    out_stream mpCellPopulationVolumesFile;

    /** Whether to use a viscosity that is linear in the cell area, rather than constant. */
    bool mUseAreaBasedDampingConstant;

    /** Non-dimensional parameter d0 for use in area-based damping constant calculations. */
    double mAreaBasedDampingConstantParameter;

    /** Whether to write cell volume and surface area (in 3D) or area and perimeter (in 2D) information to file. */
    bool mOutputVoronoiData;

    /** Whether to write the cell population volumes (in 3D) or areas (in 2D)  to file. */
    bool mOutputCellPopulationVolumes;

    /** Whether to write cells as points in VTK. */
    bool mWriteVtkAsPoints;

    /** Whether springs have variable rest lengths. */
    bool mHasVariableRestLength;

#undef COVERAGE_IGNORE // Avoid prototypes being treated as code by gcov

    /**
     * Update mIsGhostNode if required by a remesh.
     *
     * @param rMap A map between node indices before and after remesh
     */
    virtual void UpdateGhostNodesAfterReMesh(NodeMap& rMap);

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it.
     */
    virtual void Validate();

public:
#define COVERAGE_IGNORE // Avoid prototypes being treated as code by gcov

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely 1 cell for each node of the mesh.
     *
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells cells corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     * @param validate whether to validate the cell population
     */
    MeshBasedCellPopulation(MutableMesh<DIM, DIM>& rMesh,
                    std::vector<CellPtr>& rCells,
                    const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                    bool deleteMesh=false,
                    bool validate=true);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     */
    MeshBasedCellPopulation(MutableMesh<DIM, DIM>& rMesh);

    /**
     * Destructor.
     */
    ~MeshBasedCellPopulation();

    /**
     * @return reference to  mrMesh.
     */
    MutableMesh<DIM, DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const MutableMesh<DIM, DIM>& rGetMesh() const;

    /** @return mUseAreaBasedDampingConstant. */
    bool UseAreaBasedDampingConstant();

    /**
     * Set method for mOutputVoronoiData.
     *
     * @param outputVoronoiData whether to output cell area and perimeter information
     */
    void SetOutputVoronoiData(bool outputVoronoiData);

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * Overridden SetNode() method.
     *
     * Move the node with a given index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation);

    /**
     * Overridden GetDampingConstant() method that includes the
     * case of a cell-area-based damping constant.
     *
     * @param nodeIndex the global index of this node
     * @return the damping constant for the given Cell.
     */
    double GetDampingConstant(unsigned nodeIndex);

    /**
     * Set method for mUseAreaBasedDampingConstant.
     *
     * @param useAreaBasedDampingConstant  whether to use a viscosity that is linear in the cell area, rather than constant
     */
    void SetAreaBasedDampingConstant(bool useAreaBasedDampingConstant);

    /**
     * Remove all cells that are labelled as dead.
     *
     * Note that this now calls MutableMesh::DeleteNodePriorToReMesh()
     * and therefore a ReMesh(map) must be called before any element
     * information is used.
     *
     * Note also that after calling this method the cell population will be in an inconsistent state until
     * Update() is called! So don't try iterating over cells or anything like that.
     *
     * @return number of cells removed.
     */
    virtual unsigned RemoveDeadCells();

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population and update mIsGhostNode.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell - this is required for
     *  mesh-based cell populations
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell);

    /**
     * Overridden CreateOutputFiles() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param cleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     */
    void CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory);

    /**
     * Overridden CloseOutputFiles() method.
     */
    void CloseOutputFiles();

    /**
     * Overridden WriteResultsToFiles() method.
     */
    void WriteResultsToFiles();

    /**
     * Overridden Update(bool hasHadBirthsOrDeaths) method.
     * Fixes up the mappings between cells and nodes.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had Births Or Deaths
     * not needed in this cell population class
     */
    virtual void Update(bool hasHadBirthsOrDeaths=true);

       /**
     *  Tessellates when required
     *  If areas or volumes are needed for mUseAreaBasedDampingConstant or mOutputCellPopulationVolumes or mOutputCellVolumes
     *  If Voronoi data is to be output
     */
    void TessellateIfNeeded();

    /**
     * Overridden GetNode() method.
     *
     * @param index  global index of the specified Node
     *
     * @return pointer to the Node with given index.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    virtual void WriteVtkResultsToFile();

    /**
     * Write the current index and location of each node in mrMesh (including ghost nodes),
     * as well as the area and perimeter (in 2D) or volume and surface area (in 3D)
     * of its corresponding element in mpVoronoiTessellation, to mpVoronoiFile.
     */
    void WriteVoronoiResultsToFile();

    /**
     * Write the current total area (in 2D) or volume (in 3D) of mrMesh, and of the set of
     * apoptotic cells in the cell population (using mpVoronoiTessellation), to mpCellPopulationVolumesFile.
     */
    void WriteCellPopulationVolumeResultsToFile();

    /**
     * Overridden WriteCellVolumeResultsToFile() method.
     */
    void WriteCellVolumeResultsToFile();

    /**
     * Overridden GetVolumeOfCell() method.
     * 
     * @param pCell boost shared pointer to a cell
     */
    double GetVolumeOfCell(CellPtr pCell);

    /**
     * Create a Voronoi tessellation of the mesh.
     */
    void CreateVoronoiTessellation();

    /**
     * Get a reference to mpVoronoiTessellation.
     */
    VertexMesh<DIM, DIM>* GetVoronoiTessellation();

    /**
     * Get the volume (or area in 2D, or length in 1D) of the element of mpVoronoiTessellation associated with
     * the node with this global index in the Delaunay mesh.
     *
     * This method should be called instead of calling GetVoronoiTessellation()->GetVolumeOfElement()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * \todo This method is somewhat redundant following the introduction of the method GetVolumeOfCell() (see #1985).
     * 
     * @param index a node global index
     */
    double GetVolumeOfVoronoiElement(unsigned index);

    /**
     * Get the surface area of the element of mpVoronoiTessellation associated with
     * the node with this global index in the Delaunay mesh.
     *
     * This method should be called instead of calling GetVoronoiTessellation()->GetSurfaceAreaOfElement()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index a node global index
     */
    double GetSurfaceAreaOfVoronoiElement(unsigned index);

    /**
     * Get the length of the edge of mpVoronoiTessellation associated with
     * the two nodes with these global indices in the Delaunay mesh.
     *
     * This method should be called instead of calling GetVoronoiTessellation()->GetEdgeLength()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index1 a node global index
     * @param index2 a node global index
     */
    double GetVoronoiEdgeLength(unsigned index1, unsigned index2);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the cell population by calling
     * GetWidth() on the mesh.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension);

    /**
     * Iterator over edges in the mesh, which correspond to springs between cells.
     *
     * This class takes care of the logic to make sure that you consider each edge exactly once.
     */
    class SpringIterator
    {
    public:

        /**
         * Get a pointer to the node in the mesh at end A of the spring.
         */
        Node<DIM>* GetNodeA();

        /**
         * Get a pointer to the node in the mesh at end B of the spring.
         */
        Node<DIM>* GetNodeB();

        /**
         * Get the cell at end A of the spring.
         */
        CellPtr GetCellA();

        /**
         * Get the cell at end B of the spring.
         */
        CellPtr GetCellB();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther SpringIterator with which comparison is made
         */
        bool operator!=(const MeshBasedCellPopulation<DIM>::SpringIterator& rOther);

        /**
         * Prefix increment operator.
         */
        SpringIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * @param rCellPopulation the cell population
         * @param edgeIter iterator over edges in the mesh
         */
        SpringIterator(MeshBasedCellPopulation<DIM>& rCellPopulation, typename MutableMesh<DIM,DIM>::EdgeIterator edgeIter);

    private:

        /** Keep track of what edges have been visited. */
        std::set<std::set<unsigned> > mSpringsVisited;

        /** The cell population member. */
        MeshBasedCellPopulation<DIM>& mrCellPopulation;

        /** The edge iterator member. */
        typename MutableMesh<DIM, DIM>::EdgeIterator mEdgeIter;
    };

    /**
     * @return iterator pointing to the first spring in the cell population
     */
    SpringIterator SpringsBegin();

    /**
     * @return iterator pointing to one past the last spring in the cell population
     */
    SpringIterator SpringsEnd();

    /**
     * Helper method for use in debugging.
     */
    void CheckCellPointers();

    /**
     * Helper method that returns a set of pointers to two given Cells.
     * Used by the spring marking routines.
     * Elements in the returned pair are ordered by cell ID number - the
     * cell in the pair will have a smaller ID.
     *
     * @param pCell1 a Cell
     * @param pCell2 a Cell
     */
    std::pair<CellPtr,CellPtr> CreateCellPair(CellPtr pCell1, CellPtr pCell2);

    /**
     * @param rCellPair a set of pointers to Cells
     *
     * @return whether the spring between two given cells is marked.
     */
    bool IsMarkedSpring(const std::pair<CellPtr,CellPtr>& rCellPair);

    /**
     * Mark the spring between the given cells.
     *
     * @param rCellPair a set of pointers to Cells
     */
    void MarkSpring(std::pair<CellPtr,CellPtr>& rCellPair);

    /**
     * Stop marking the spring between the given cells.
     *
     * @param rCellPair a set of pointers to Cells
     */
    void UnmarkSpring(std::pair<CellPtr,CellPtr>& rCellPair);

    /**
     * @return mAreaBasedDampingConstantParameter
     */
    double GetAreaBasedDampingConstantParameter();

    /**
     * Set mAreaBasedDampingConstantParameter.
     *
     * @param areaBasedDampingConstantParameter the new value of mAreaBasedDampingConstantParameter
     */
    void SetAreaBasedDampingConstantParameter(double areaBasedDampingConstantParameter);

    /**
     * @return mOutputVoronoiData
     */
    bool GetOutputVoronoiData();

    /**
     * @return mOutputCellPopulationVolumes
     */
    bool GetOutputCellPopulationVolumes();

    /**
     * Set mOutputCellPopulationVolumes.
     *
     * @param outputCellPopulationVolumes the new value of mOutputCellPopulationVolumes
     */
    void SetOutputCellPopulationVolumes(bool outputCellPopulationVolumes);

    /**
     * Outputs CellPopulation parameters to file
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Set mWriteVtkAsPoints.
     *
     * @param writeVtkAsPoints whether to write cells as points in VTK
     */
    void SetWriteVtkAsPoints(bool writeVtkAsPoints);

    /**
     * @return mWriteVtkAsPoints.
     */
    bool GetWriteVtkAsPoints();

    /**
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);

    /**
     * Populate mSpringRestLengths by looping over all springs and calculating the current length
     */
    void CalculateRestLengths();

    /**
     *  Return the rest length for a given spring
     *
     *  @param indexA index of first node in pair
     *  @param indexB index of second node in pair
     */
    double GetRestLength(unsigned indexA, unsigned indexB);
};
#undef COVERAGE_IGNORE // Avoid prototypes being treated as code by gcov

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MeshBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const MutableMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a MeshBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MeshBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MutableMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MeshBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*MESHBASEDCELLPOPULATION_HPP_*/
