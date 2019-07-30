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


#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCellInterface.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeartGeometryInformation.hpp"
#include "HeartRegionCodes.hpp"
#include "ZeroStimulus.hpp"

/**
 * A factory to ease creating cardiac cell objects for use in a mono/bidomain simulation.
 *
 * The user should implement their own concrete class, in particular implementing
 * CreateCardiacCellForTissueNode(Node*), which should return the cell corresponding to a
 * given node. The user should also implement GetNumberOfCells() if this isn't equal
 * to the number of nodes. FinaliseCellCreation() can be used to (eg) add stimuli to
 * certain cells after they have been created.
 *
 * This class saves the user having to create cells in parallel, that work is done
 * by the pde instead.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractCardiacCellFactory
{
private:
    /** The mesh is automatically set in MonodomainProblem and BidomainProblem.
     *  This member variable should be accessed through GetMesh(), which will check if it has been set before
     *  and throw an exception otherwise.*/
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /**
     * A pointer to an HeartGeometryInformation information object
     * Can be accessed via get and set methods in this class.
     */
    HeartGeometryInformation<SPACE_DIM>* mpHeartGeometryInformation;

protected:
    /** For use at un-stimulated cells. */
    boost::shared_ptr<ZeroStimulus> mpZeroStimulus;
    /** The solver to give each of the cells */
    boost::shared_ptr<AbstractIvpOdeSolver> mpSolver;

public:
    /**
     * @return a newly created cell object for the given node.
     *
     * The default implementation checks whether the node is in the bath (in which
     * case a pointer to a (unique) fake cell is returned) and if not, calls
     * CreateCardiacCellForTissueNode (which must be defined by subclasses).
     *
     * @param pNode  Pointer to node object.
     */
    virtual AbstractCardiacCellInterface* CreateCardiacCellForNode(Node<SPACE_DIM>* pNode);

    /**
     * Must be overridden by subclasses to return a cell object for the given node.
     *
     * @param pNode  Pointer to node object.
     * @return a newly created cell object for the given tissue node.
     */
    virtual AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<SPACE_DIM>* pNode)=0;

    /**
     * May be overridden by subclasses to perform any necessary work after all cells
     * have been created.
     *
     * @param pCellsDistributed  Pointer to a vector of cardiac cell pointers.
     * @param lo  Lowest index owned by this process.
     * @param hi  Highest index owned by this process.
     */
    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCellInterface* >* pCellsDistributed,
                                      unsigned lo, unsigned hi);

    /**
     * Method that fills in the vector of heterogeneity areas with the NodesLists
     * that correspond to a given layer (implemented in subclasses)
     */
    virtual void FillInCellularTransmuralAreas();

    /**
     * @return  The number of cells
     */
    virtual unsigned GetNumberOfCells();

    /**
     * Default constructor.
     *
     * @param pSolver  the ODE solver to use to simulate this cell.
     */
    AbstractCardiacCellFactory(boost::shared_ptr<AbstractIvpOdeSolver> pSolver = boost::shared_ptr<AbstractIvpOdeSolver>(new EulerIvpOdeSolver));
    /**
     * Destructor: free solver, zero stimulus and fake bath cell.
     */
    virtual ~AbstractCardiacCellFactory();

    /**
     * @param pMesh  the mesh for which to create cardiac cells.
     */
    virtual void SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * @return  the mesh used to create the cells.
     */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* GetMesh();

    /**
     * Set the HeartGeometryInformation object
     *
     * @param pHeartGeometryInformation the HeartGeometryInformation object that is to be set
     */
    void SetHeartGeometryInformation(HeartGeometryInformation<SPACE_DIM>* pHeartGeometryInformation);

    /**
     * @return the HeartGeometryInformation object
     */
    HeartGeometryInformation<SPACE_DIM>* GetHeartGeometryInformation();
};

#endif /*ABSTRACTCARDIACCELLFACTORY_HPP_*/
