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


#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCell.hpp"
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
 * CreateCardiacCellForTissueNode(unsigned), which should return the cell corresponding to a
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
     * Create a cell object for the given node.
     *
     * The default implementation checks whether the node is in the bath (in which
     * case a pointer to a (unique) fake cell is returned) and if not, calls
     * CreateCardiacCellForTissueNode (which must be defined by subclasses).
     *
     * @param nodeIndex  Global node index.
     */
    virtual AbstractCardiacCell* CreateCardiacCellForNode(unsigned nodeIndex);

    /**
     * Must be overridden by subclasses to return a cell object for the given node.
     *
     * @param nodeIndex  Global node index.
     */
    virtual AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)=0;

    /**
     * May be overridden by subclasses to perform any necessary work after all cells
     * have been created.
     *
     * @param pCellsDistributed  Pointer to a vector of cardiac cell pointers.
     * @param lo  Lowest index owned by this process.
     * @param hi  Highest index owned by this process.
     */
    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed,
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

