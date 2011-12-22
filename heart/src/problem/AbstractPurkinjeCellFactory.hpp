/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTPURKINJECELLFACTORY_HPP_
#define ABSTRACTPURKINJECELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "MixedDimensionMesh.hpp"
#include "FakeBathCell.hpp"
/**
 *  Subclass for also creating cell models for Purkinje cells.
 *  The user has to implement the CreatePurkinjeCellForTissueNode() method.
 *
 *  The dimensions should be 2 or 3.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractPurkinjeCellFactory : public AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>
{
protected:
    /** Saved pointer to the mixed dimension mesh */
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* mpMixedDimensionMesh;

    /** A set of local purkinje node indices */
    std::set<unsigned> mLocalPurkinjeNodes;
    /**
     * Must be overridden by subclasses to return a Purkinje cell object for the given node.
     * @param nodeIndex  Global node index.
     */
    virtual AbstractCardiacCell* CreatePurkinjeCellForTissueNode(unsigned nodeIndex)=0;

public:

    /* Constructor does nothing */
    AbstractPurkinjeCellFactory();


    /** Overridden set mesh which must take a MixedDimensionMesh
     *  @param pMesh Pointer to the mesh. */
    void SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * Create a cell object for the given node.
     *
     * The default implementation checks whether the node is in a Purkinje node, in which
     * case it calls CreatePurkinjeCellForTissueNode (which must be defined by subclasses),
     * otherwise it returns a pointer to a (unique) fake cell
     *
     * @param nodeIndex  Global node index.
     */
    AbstractCardiacCell* CreatePurkinjeCellForNode(unsigned nodeIndex);

    /**
     * May be overridden by subclasses to perform any necessary work after all Purkinje cells
     * have been created.
     *
     * @param pPurkinjeCellsDistributed  Pointer to a vector of Purkinje cell pointers.
     * @param lo  Lowest index owned by this process.
     * @param hi  Highest index owned by this process.
     */
    virtual void FinalisePurkinjeCellCreation(std::vector< AbstractCardiacCell* >* pPurkinjeCellsDistributed,
                                              unsigned lo, unsigned hi)
    {
    }

    /**
     *  Get the mixed dimension mesh (for possible use in CreatePurkinjeCellForTissueNode()).
     *  Note: GetMesh() just returns a pointer to an AbstractTetrahedralMesh.
     */
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* GetMixedDimensionMesh();
};



#endif // ABSTRACTPURKINJECELLFACTORY_HPP_
