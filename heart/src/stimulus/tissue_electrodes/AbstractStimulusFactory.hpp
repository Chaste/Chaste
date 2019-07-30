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
     * @return a stimulus object for the given node.
     * Default implementation here returns a zero stimulus.
     *
     *  May be overridden by child classes.
     *
     * @param pNode  pointer to the node object
     */
    virtual boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<SPACE_DIM>* pNode);

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
     * @return access to the variable mGroundedRegions which
     * stores the regions to be grounded
     */
    std::vector<AbstractChasteRegion<SPACE_DIM>* > GetRegionsToBeGrounded();
};

#endif /*ABSTRACTSTIMULUSFACTORY_HPP_*/

