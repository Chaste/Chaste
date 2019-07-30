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

#ifndef _TETRAHEDRALMESH_HPP_
#define _TETRAHEDRALMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "UblasVectorInclude.hpp"
#include "UblasMatrixInclude.hpp"

#include <vector>
#include <string>
#include <set>

#include "AbstractTetrahedralMesh.hpp"
#include "AbstractMeshReader.hpp"
#include "ChastePoint.hpp"


struct triangulateio; /**< Forward declaration for triangle helper methods (used in MutableMesh QuadraticMesh)*/

//////////////////////////////////////////////////////////////////////////
//   DECLARATION
//////////////////////////////////////////////////////////////////////////

/**
 * A concrete tetrahedral mesh class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TetrahedralMesh : public AbstractTetrahedralMesh< ELEMENT_DIM, SPACE_DIM>
{
    friend class TestTetrahedralMesh; // to give access to private methods (not variables)
    friend class TestCryptSimulation2d; // to give access to private methods (not variables)
private:
    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

protected:

    /**
     * Overridden solve node mapping method.
     *
     * @param index the global index of the node
     * @return local index
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Overridden solve element mapping method.
     *
     * @param index the global index of the element
     * @return local index
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Overridden solve boundary element mapping method.
     *
     * @param index the global index of the boundary element
     * @return local index
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /** Vector storing the weighted direction for each element in the mesh. */
    std::vector< c_vector<double, SPACE_DIM> > mElementWeightedDirections;

    /** Vector storing the Jacobian matrix for each element in the mesh. */
    std::vector< c_matrix<double, SPACE_DIM, ELEMENT_DIM> > mElementJacobians;

    /** Vector storing the inverse Jacobian matrix for each element in the mesh. */
    std::vector< c_matrix<double, ELEMENT_DIM, SPACE_DIM> > mElementInverseJacobians;

    /** Vector storing the Jacobian determinant for each element in the mesh. */
    std::vector<double> mElementJacobianDeterminants;

    /** Vector storing the weighted direction for each boundary element in the mesh. */
    std::vector< c_vector<double, SPACE_DIM> > mBoundaryElementWeightedDirections;

    /** Vector storing the determinant of the Jacobian matrix for each boundary element in the mesh. */
    std::vector<double> mBoundaryElementJacobianDeterminants;

    /**
     * Export the mesh (currently only the nodes) to an external mesher
     * This is determined at compile time when the MESHER_IO template is
     * instantiated to either
     *   - triangulateio (for triangle remesher in 2D)
     *   - tetgenio (for tetgen remesher in 3D)
     * Since conditional compilation is not allowed, care must be taken to only use
     * common data members in this method
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     * @param mesherInput is a triangulateio or tetgenio class (decided at compile time)
     *  Note that only nodes are exported and thus any late mesh is based on the convex hull
     * @param elementList is a pointer to either mesherInput.trianglelist or mesherInput.tetrahedronlist (used when we are not remeshing, but converting the form of an existing mesh)
     * Note that this should have been re-malloced and mesherInput.numberoftriangles or mesherInput.numberoftetrahedra should be set prior to the call when elementList is non-NULL
     *
     */
    template <class MESHER_IO>
    void ExportToMesher(NodeMap& map, MESHER_IO& mesherInput, int *elementList=nullptr);

    /**
     * Import the mesh from an external mesher
     * This is determined at compile time when the MESHER_IO template is
     * instantiated to either
     *   - triangulateio (for triangle remesher in 2D)
     *   - tetgenio (for tetgen remesher in 3D)
     * Since conditional compilation is not allowed, care must be taken to only use
     * common data members in this method
     * @param mesherOutput is a triangulateio or tetgenio class (decided at compile time)
     * @param numberOfElements is a copy of either mesherOutput.numberoftriangles or mesherOutput.numberoftetrahedra
     * @param elementList is a pointer to either mesherOutput.trianglelist or mesherOutput.tetrahedronlist
     * @param numberOfFaces is a copy of either mesherOutput.edges or mesherOutput.numberoftrifaces
     * @param faceList is a pointer to either mesherOutput.edgelist or mesherOutput.trifacelist
     * @param edgeMarkerList is a pointer to either mesherOutput.edgemarkerlist or NULL
     * \todo #1545: (or add arguments ...)
     *
     */
    template <class MESHER_IO>
    void ImportFromMesher(MESHER_IO& mesherOutput, unsigned numberOfElements, int *elementList, unsigned numberOfFaces, int *faceList, int *edgeMarkerList);


    /**
     * Convenience method to tidy up a triangleio data structure before use
     * @param mesherIo is a triangulateio class
     */
    void InitialiseTriangulateIo(triangulateio& mesherIo);

    /**
     * Convenience method to tidy up a triangleio data structure after use
     * @param mesherIo is a triangulateio class
     */
    void FreeTriangulateIo(triangulateio& mesherIo);

public:

    /**
     * Constructor.
     */
    TetrahedralMesh();

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

    /**
     * Read in the number of nodes per processor from file.
     *
     * @param rNodesPerProcessorFile the name of the file
     */
    void ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile);

    /**
     * Check whether mesh is conforming
     * Conforming (defn.): the intersection of two elements should be
     * either the empty set, a vertex, an edge or a face.
     *
     * It may be possible to construct non-conforming meshes which contain
     * internal faces owned by only one element: two coplanar triangular faces
     * of two elements form a square, but the same square on the adjacent pair of
     * elements is formed by splitting the diagonal the other way.
     *
     * @return false if there are any orphaned internal faces
     */
    bool CheckIsConforming();

    /**
     * @return the volume of the mesh, calculated by adding the determinant of each element
     * and dividing by n!, where n is the element dimension.
     */
    double GetVolume();

    /**
     * @return the surface area of the mesh.
     */
    double GetSurfaceArea();

    /**
     * Overridden RefreshMesh method. This method calls RefreshJacobianCachedData.
     */
    void RefreshMesh();

    /**
     * Permute the nodes randomly so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     */
    void PermuteNodes();

    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     * @param perm is a vector containing the new indices
     */
    void PermuteNodes(const std::vector<unsigned>& perm);

    /**
     * @return the element index for the first element that contains a test point. Like GetContainingElementIndex
     * but uses the user given element (M say) as the first element checked, and then checks M+1,M+2,..,Ne,0,1..
     *
     * @param rTestPoint reference to the point
     * @param startingElementGuess Which element to try first.
     * @param strict  Should the element returned contain the point in the interior and
     *      not on an edge/face/vertex (default = not strict)
     */
    unsigned GetContainingElementIndexWithInitialGuess(const ChastePoint<SPACE_DIM>& rTestPoint, unsigned startingElementGuess, bool strict=false);

    /**
     * @return the element index for an element is closest to the testPoint.
     *
     * "Closest" means that the minimum interpolation weights for the testPoint are
     * maximised for this element.
     *
     * @param rTestPoint reference to the point
     */
    unsigned GetNearestElementIndex(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * @return all element indices for elements that are known to contain a test point.
     *
     * @param rTestPoint reference to the point
     */
    std::vector<unsigned> GetContainingElementIndices(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * Clear all the data in the mesh.
     */
    virtual void Clear();

    /**
     * @return calculated the angle between the node at indexB and the x axis about
     * the node at indexA. The angle returned is in the range (-pi,pi].
     *
     * @param indexA a node index
     * @param indexB a node index
     */
    double GetAngleBetweenNodes(unsigned indexA, unsigned indexB);

    /** Update mElementJacobians, mElementWeightedDirections and mBoundaryElementWeightedDirections. */
    virtual void RefreshJacobianCachedData();

    /**
     * @return the Jacobian matrix and its determinant for a given element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rJacobian the Jacobian matrix
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     */
    virtual void GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double& rJacobianDeterminant) const;

    /**
     * @return the Jacobian matrix, its inverse and its determinant for a given element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rJacobian the Jacobian matrix
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     * @param rInverseJacobian the inverse Jacobian matrix
     */
    virtual void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const;

    /**
     * @return the weighted direction and the determinant of the Jacobian for a given element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rWeightedDirection the weighted direction
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     */
    virtual void GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const;

    /**
     * @return the weighted direction and the determinant of the Jacobian for a given boundary element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rWeightedDirection the weighted direction
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     */
    virtual void GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const;

    /**
     * Iterator over edges in the mesh.
     *
     * This class takes care of the logic to make sure that you consider each edge exactly once.
     *
     */
    class EdgeIterator
    {
    public:
        /**
         * @return a pointer to the node in the mesh at end A of the edge.
         */
        Node<SPACE_DIM>* GetNodeA();
        /**
         * @return a pointer to the node in the mesh at end B of the edge.
         */
        Node<SPACE_DIM>* GetNodeB();

        /**
         * @return comparison not-equal-to.
         *
         * @param rOther edge iterator with which comparison is made
         */
        bool operator!=(const typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        EdgeIterator& operator++();

        /**
         * Constructor for a new edge iterator.
         *
         * @param rMesh  The mesh
         * @param elemIndex  An element index
         */
        EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex);

    private:
        /**
         * Keep track of what edges have been visited.
         * Each edge is stored as a pair of ordered indices.
         */
        std::set< std::pair<unsigned, unsigned> > mEdgesVisited;

        TetrahedralMesh& mrMesh;   /**< The mesh. */

        unsigned mElemIndex;       /**< Element index. */
        unsigned mNodeALocalIndex; /**< Index of one node on the edge. */
        unsigned mNodeBLocalIndex; /**< Index of the other node on the edge. */
    };

    /**
     * @return iterator pointing to the first edge (i.e. connection between 2 nodes) of the mesh
     */
    EdgeIterator EdgesBegin();

    /**
     * @return iterator pointing to one past the last edge (i.e. connection between 2 nodes)
     * of the mesh
     */
    EdgeIterator EdgesEnd();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TetrahedralMesh)

#endif //_TETRAHEDRALMESH_HPP_
