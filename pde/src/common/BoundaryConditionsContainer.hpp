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

#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/map.hpp>

#include "AbstractBoundaryConditionsContainer.hpp"
#include "AbstractBoundaryCondition.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "BoundaryElement.hpp"
#include "Node.hpp"
#include "LinearSystem.hpp"
#include "PetscException.hpp"
#include "ChastePoint.hpp"
#include "ConstBoundaryCondition.hpp"
#include "DistributedVectorFactory.hpp"

/**
 * Boundary Conditions Container.
 *
 * This class contains a list of nodes on the Dirichlet boundary and associated Dirichlet
 * boundary conditions, and a list of surface elements on the Neumann boundary and associated
 * Neumann boundary conditions.
 *
 * \todo #1321
 * Various operations are currently very inefficient - there is certainly scope for
 * optimisation here!
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class BoundaryConditionsContainer : public AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
public:

    /** Type of a read-only iterator over Neumann boundary conditions. */
    typedef typename std::map< const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>*, const AbstractBoundaryCondition<SPACE_DIM>* >::const_iterator
        NeumannMapIterator;

    /** Base class type. */
    typedef AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> BaseClassType;

private:

    std::map< const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>* >*
        mpNeumannMap[PROBLEM_DIM]; /**< List (map) of Neumann boundary conditions. */

    /** Nodes to be identified when periodic boundary conditions are applied */
    std::map< const Node<SPACE_DIM> *, const Node<SPACE_DIM> * >*  mpPeriodicBcMap[PROBLEM_DIM];

    /**
     * Neumann boundary condition iterator.
     */
    NeumannMapIterator mLastNeumannCondition[PROBLEM_DIM];

    /**
     * Array storing whether there are any Neumann boundary conditions for each unknown.
     */
    bool mAnyNonZeroNeumannConditionsForUnknown[PROBLEM_DIM];

    /** A zero boundary condition, used for other unknowns in ApplyNeumannBoundaryCondition */
    ConstBoundaryCondition<SPACE_DIM>* mpZeroBoundaryCondition;

    /** Whether the contents of this container were originally loaded from an archive. */
    bool mLoadedFromArchive;

public:

    /**
     * Constructor calls base constuctor and allocates memory for the Neumann boundary
     * conditions lists.
     *
     * @param deleteConditions whether to delete BCs in destructor (defaults to true)
     */
    BoundaryConditionsContainer(bool deleteConditions=true);

    /**
     * Note that the destructor will delete memory for each boundary condition object, as
     * well as for the internal bookkeeping of this class.
     */
    ~BoundaryConditionsContainer();

    /**
     * Add a Dirichlet boundary condition specifying two parameters, a pointer to a node,
     * and a pointer to a boundary condition object associated with that node.
     *
     * The destructor for the BoundaryConditionsContainer will destroy the boundary
     * conditions objects.
     *
     * @param pBoundaryNode Pointer to a node on the boundary.
     * @param pBoundaryCondition Pointer to the Dirichlet boundary condition at that node.
     * @param indexOfUnknown defaults to 0
     * @param checkIfBoundaryNode defaults to true
     */
    void AddDirichletBoundaryCondition(const Node<SPACE_DIM>* pBoundaryNode,
                                       const AbstractBoundaryCondition<SPACE_DIM>* pBoundaryCondition,
                                       unsigned indexOfUnknown = 0,
                                       bool checkIfBoundaryNode = true);

    /**
     * Add a Neumann boundary condition specifying two parameters, a pointer to a
     * surface element, and a pointer to a boundary condition object associated with
     * that element.
     *
     * The destructor for the BoundaryConditionsContainer will destroy the boundary
     * conditions objects.
     *
     * Note that the value of a Neumann boundary condition should specify
     * D * grad(u).n, not just grad(u).n.
     *
     * Take care if using non-zero Neumann boundary conditions in 1d. If applied at
     * the left hand end you need to multiply the value by -1 to get the right answer.
     *
     * @param pBoundaryElement Pointer to an element on the boundary
     * @param pBoundaryCondition Pointer to the Neumann boundary condition on that element
     * @param indexOfUnknown defaults to 0
     */
    void AddNeumannBoundaryCondition(const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* pBoundaryElement,
                                     const AbstractBoundaryCondition<SPACE_DIM>* pBoundaryCondition,
                                     unsigned indexOfUnknown = 0);


    /**
     *  Add a periodic boundary condition: provide two nodes to be identified when solving
     *  @param pNode1 node 1
     *  @param pNode2 node 2
     *
     *  This method identifies the nodes for all unknowns, so doesn't have to be called for each unknown.
     */
    void AddPeriodicBoundaryCondition(const Node<SPACE_DIM>* pNode1,
                                      const Node<SPACE_DIM>* pNode2);


    /**
     * This function defines zero Dirichlet boundary conditions on every boundary node
     * of the mesh.
     *
     * @param pMesh Pointer to a mesh object, from which we extract the boundary
     * @param indexOfUnknown defaults to 0
     */
    void DefineZeroDirichletOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                           unsigned indexOfUnknown = 0);

    /**
     * This function defines constant Dirichlet boundary conditions on every boundary node
     * of the mesh.
     *
     * @param pMesh Pointer to a mesh object, from which we extract the boundary
     * @param value the value of the constant Dirichlet boundary condition
     * @param indexOfUnknown defaults to 0
     */
    void DefineConstantDirichletOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                               double value,
                                               unsigned indexOfUnknown = 0);

    /**
     * This function defines zero Neumann boundary conditions on every boundary element
     * of the mesh.
     *
     * @param pMesh Pointer to a mesh object, from which we extract the boundary
     * @param indexOfUnknown defaults to 0
     */
    void DefineZeroNeumannOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                         unsigned indexOfUnknown = 0);

    /**
     *  Alter the given linear system to satisfy Dirichlet boundary conditions.
     *
     *  If the number of unknowns is greater than one, it is assumed the solution vector is
     *  of the form (in the case of two unknowns u and v, and N nodes):
     *  solnvec = (U_1, V_1, U_2, V_2, ...., U_N, V_N)
     *
     *  @param rLinearSystem Linear system on which to apply boundary conditions
     *
     *  @param applyToMatrix This optional parameter can be set as false to
     *  ensure that the matrix of the linear system is not updated. To
     *  be used when the matrix does not change between time steps.
     *  @param applyToRhsVector Similarly, whether to apply the changes to the RHS vector (b in Ax=b).
     */
    void ApplyDirichletToLinearProblem(LinearSystem& rLinearSystem,
                                       bool applyToMatrix = true,
                                       bool applyToRhsVector = true);

    /**
     *  Alter the given linear system to satisfy periodic boundary conditions.
     *
     *  For one of the two nodes that have been identified, the row corresponding to the
     *  unknown which has periodic BCs, for one of the nodes, ie replaced with
     *  [0 0 0 ... -1 0 0 .. 0 1 0 .. 0]
     *  where the 1 is the diagonal entry and the -1 on the column corresponding to that
     *  unknown and the other node.
     *
     *  The entry in the RHS vector is zeroed.
     *
     *  @param rLinearSystem Linear system on which to apply boundary conditions
     *
     *  @param applyToMatrix This optional parameter can be set as false to
     *  ensure that the matrix of the linear system is not updated. To
     *  be used when the matrix does not change between time steps.
     *  @param applyToRhsVector Similarly, whether to apply the changes to the RHS vector (b in Ax=b).
     *
     */
    void ApplyPeriodicBcsToLinearProblem(LinearSystem& rLinearSystem,
                                         bool applyToMatrix = true,
                                         bool applyToRhsVector = true);

    /**
     * Alter the residual vector for a nonlinear system to satisfy
     * Dirichlet boundary conditions.
     *
     * If the number of unknowns is greater than one, it is assumed the solution vector is
     * of the form (in the case of two unknowns u and v, and N nodes):
     * solnvec = (U_1, V_1, U_2, V_2, ...., U_N, V_N)
     *
     * @param currentSolution
     * @param residual
     * @param rFactory  the factory to use to create DistributedVector objects
     */
    void ApplyDirichletToNonlinearResidual(const Vec currentSolution, Vec residual,
                                           DistributedVectorFactory& rFactory);

    /**
     * Alter the Jacobian matrix vector for a nonlinear system to satisfy
     * Dirichlet boundary conditions.
     *
     * If the number of unknowns is greater than one, it is assumed the solution vector is
     * of the form (in the case of two unknowns u and v, and N nodes):
     * solnvec = (U_1, V_1, U_2, V_2, ...., U_N, V_N)
     *
     * @param jacobian
     */
    void ApplyDirichletToNonlinearJacobian(Mat jacobian);

    /**
     * Check that we have boundary conditions defined everywhere on mesh boundary.
     *
     * We iterate over all surface elements, and check either that they have an
     * associated Neumann condition, or that each node in the element has an
     * associated Dirichlet condition.
     *
     * @param pMesh Pointer to the mesh to check for validity.
     * @return true iff all boundaries have boundary conditions defined.
     */
    bool Validate(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * @return value of Neumann boundary condition at a specified point in a given surface element
     *
     * It is up to the user to ensure that the point x is contained in the surface element.
     *
     * @param pSurfaceElement pointer to a boundary element
     * @param rX a point
     * @param indexOfUnknown defaults to 0
     */
    double GetNeumannBCValue(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* pSurfaceElement,
                             const ChastePoint<SPACE_DIM>& rX,
                             unsigned indexOfUnknown = 0);

    /**
     * @return true if there is a Neumann boundary condition defined on the given element.
     *
     * \todo #1321
     * This is a horrendously inefficient fix. Perhaps have flag in element object?
     *
     * @param pSurfaceElement pointer to a boundary element
     * @param indexOfUnknown defaults to 0
     */
    bool HasNeumannBoundaryCondition(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* pSurfaceElement,
                                     unsigned indexOfUnknown = 0);

    /**
     * @return whether there are any non-zero Neumann boundary conditions
     */
    bool AnyNonZeroNeumannConditions();

    /**
     * @return iterator pointing to the first Neumann boundary condition
     */
    NeumannMapIterator BeginNeumann();

    /**
     * @return iterator pointing to one past the last Neumann boundary condition
     */
    NeumannMapIterator EndNeumann();

    /**
     * Load a collection of boundary conditions from an archive.
     *
     * @note We assume this collection is empty prior to being called.  If it is not, any boundary
     * conditions already present may get replaced by conditions loaded from the archive, which may
     * lead to a memory leak.
     *
     * This method only loads data if #mLoadedFromArchive is false, to allow for multiple pointers
     * to the same container to be handled correctly.  It sets #mLoadedFromArchive when done.
     *
     * @param archive  the archive to load from
     * @param pMesh  the mesh to use to resolve Node and BoundaryElement indices
     */
    template <class Archive>
    void LoadFromArchive(Archive & archive, AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    {
        if (mLoadedFromArchive)
        {
            return;
        }

        MergeFromArchive(archive, pMesh);
    }

    /**
     * Load extra boundary conditions from an archive to add to this collection.
     *
     * Multiple pointers to the same container need to be handled by the caller - we assume there
     * will be conditions to load.  Sets #mLoadedFromArchive when done.
     *
     * @param archive  the archive to load from
     * @param pMesh  the mesh to use to resolve Node and BoundaryElement indices
     */
    template <class Archive>
    void MergeFromArchive(Archive & archive, AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save this container and its contents.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const;

    /**
     * Load this container, but not its content.
     *
     * Objects loading a boundary conditions container should call LoadFromArchive
     * on the new object immediately after loading it from the archive.
     *
     * Note that boundary conditions should be saved to the ProcessSpecificArchive,
     * since if a DistributedTetrahedralMesh is used each process will only know a
     * portion of the mesh, and hence a portion of the boundary conditions.
     *
     * Extra care needs to be taken when migrating to ensure that boundary conditions
     * are loaded appropriately.  See BidomainProblem::LoadExtraArchiveForBidomain
     * and AbstractCardiacProblem::LoadExtraArchive for examples.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
template<class Archive>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::save(
        Archive & archive, const unsigned int version) const
{
    typedef typename std::map<unsigned, const AbstractBoundaryCondition<SPACE_DIM> *> archive_map_type;

    // Save Dirichlet conditions
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        archive_map_type bc_map;
        typename BaseClassType::DirichletIteratorType it = this->mpDirichletMap[index_of_unknown]->begin();
        while (it != this->mpDirichletMap[index_of_unknown]->end() )
        {
            unsigned node_index = it->first->GetIndex();
            const AbstractBoundaryCondition<SPACE_DIM> * p_cond = it->second;
            bc_map[node_index] = p_cond;

            it++;
        }
        archive & bc_map;
    }

    // Save Neumann conditions
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        archive_map_type bc_map;
        for (NeumannMapIterator it = mpNeumannMap[index_of_unknown]->begin();
             it != mpNeumannMap[index_of_unknown]->end();
             ++it)
        {
            unsigned elem_index = it->first->GetIndex();
            const AbstractBoundaryCondition<SPACE_DIM>* p_cond = it->second;
            bc_map[elem_index] = p_cond;
        }
        archive & bc_map;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
template<class Archive>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::MergeFromArchive(
        Archive & archive, AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mLoadedFromArchive = true;

    typedef typename std::map<unsigned, AbstractBoundaryCondition<SPACE_DIM>*> archive_map_type;

    // Keep track of conditions that might need deleting
    std::set<const AbstractBoundaryCondition<SPACE_DIM>*> maybe_unused_bcs;
    std::set<const AbstractBoundaryCondition<SPACE_DIM>*> used_bcs;

    // Load Dirichlet conditions
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        archive_map_type bc_map;
        archive & bc_map;
        for (typename archive_map_type::iterator it = bc_map.begin();
             it != bc_map.end();
             ++it)
        {
            unsigned node_index = it->first;
            this->mHasDirichletBCs=true;  //We know that a Dirichlet is being added, even if not by this process
            Node<SPACE_DIM>* p_node;
            try
            {
                p_node = pMesh->GetNodeFromPrePermutationIndex(node_index);
            }
            catch (Exception&)
            {
                // It's a distributed mesh and we don't own this node - skip to the next BC
                maybe_unused_bcs.insert(it->second);
                continue;
            }
            AddDirichletBoundaryCondition(p_node, it->second, index_of_unknown, false);
            used_bcs.insert(it->second);
        }
    }
    this->mCheckedAndCommunicatedIfDirichletBcs=true; // Whether the Dirichlet BCC was empty or not, all processes know the status.

    // Load Neumann conditions
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        archive_map_type bc_map;
        archive & bc_map;
        for (typename archive_map_type::iterator it = bc_map.begin();
             it != bc_map.end();
             ++it)
        {
            unsigned boundary_element_index = it->first;
            BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element;
            try
            {
                p_boundary_element = pMesh->GetBoundaryElement(boundary_element_index);
            }
            catch (Exception&)
            {
                // It's a distributed mesh and we don't own this element - skip to the next BC
                maybe_unused_bcs.insert(it->second);
                continue;
            }
            AddNeumannBoundaryCondition(p_boundary_element, it->second, index_of_unknown);
            used_bcs.insert(it->second);
        }
    }

    // Free any unused BCs
    for (typename std::set<const AbstractBoundaryCondition<SPACE_DIM>*>::iterator it=maybe_unused_bcs.begin();
         it != maybe_unused_bcs.end();
         ++it)
    {
        typename std::set<const AbstractBoundaryCondition<SPACE_DIM>*>::iterator used = used_bcs.find(*it);
        if (used == used_bcs.end())
        {
            delete (*it);
        }
    }
}

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
