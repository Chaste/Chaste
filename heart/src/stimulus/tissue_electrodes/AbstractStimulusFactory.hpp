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


#ifndef ABSTRACTSTIMULUSFACTORY_HPP_
#define ABSTRACTSTIMULUSFACTORY_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractStimulusFunction.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"
#include "AbstractChasteRegion.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "Node.hpp"

/**
 * Abstract class that specifies a stimulus factory.
 * The core method is CreateStimulusForNode which creates a stimulus object at each node in the mesh.
 * The default implementation (here) specifies a ZeroStimulus. child classes may override this method.
 *
 * The other crucial method is SetCompatibleExtracellularStimulus which ensures that the stimulus specified obeys compatibility conditions.
 * The implementation here is empty (compatible with ZeroStimulus which implicitly obeys compatibility conditions).
 * child classes may over-ride the method.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractStimulusFactory
{
protected:

    /**
     * The mesh, set by the problem class
     */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /**Vector of regions of nodes to be grounded*/
    std::vector<AbstractChasteRegion<SPACE_DIM>* > mGroundedRegions;

public:

    /**
     * Create a stimulus object for the given node.
     * Default implementation here returns a zero stimulus.
     *
     *  May be overridden by child classes.
     *
     * @param nodeIndex  Global node index.
     */
    virtual boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(unsigned nodeIndex);

    /**
     * @return  The number of cells
     */
    virtual unsigned GetNumberOfCells();

    /**
     * Empty default implementation that is compatible with the
     * default implementation of CreateStimulusForNode (i.e., empty implementation for zero stimulus).
     *
     * If needed, over-ride this
     */
    virtual void SetCompatibleExtracellularStimulus();

    /**
     * Default constructor.
     */
    AbstractStimulusFactory();

    /**
     * Destructor
     */
     virtual ~AbstractStimulusFactory();

    /**
     * @param pMesh  the mesh for which to create stimuli.
     */
    void SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * @return  the mesh used to create the stimuli.
     */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* GetMesh();

//    /*
//     * Allow specification of a chaste region where nodes have extracellular potential grounded
//     *
//     * @param pRegion the region to be grounded.
//     */
//    void SetRegionToBeGrounded(AbstractChasteRegion<SPACE_DIM>* pRegion);

    /**
     * Access to the variable mGroundedRegions which
     * stores the regions to be grounded
     */
    std::vector<AbstractChasteRegion<SPACE_DIM>* > GetRegionsToBeGrounded();
};

#endif /*ABSTRACTSTIMULUSFACTORY_HPP_*/

