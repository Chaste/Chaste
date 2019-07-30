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


#ifndef ABSTRACTACINARUNITFACTORY_HPP_
#define ABSTRACTACINARUNITFACTORY_HPP_

#include "AbstractAcinarUnit.hpp"
#include "AbstractTetrahedralMesh.hpp"

/**
 * A factory to ease creating acinar unit objects for use in a ventilation simulation.
 *
 * The user should implement their own concrete class, in particular implementing
 * CreateAcinarUnitForNode(Node*), which should return the acinar unit corresponding to a
 * given node. FinaliseAcinarUnitCreation() can be used to (eg) add global compliance (such
 * as chest wall) to the models
 */
class AbstractAcinarUnitFactory
{
private:
    /** The mesh is automatically set in VentilationProblem
     *  This member variable should be accessed through GetMesh(), which will check if it has been set before
     *  and throw an exception otherwise.*/
    AbstractTetrahedralMesh<1,3>* mpMesh;

public:
    /**
     * @return a newly acinar unit object for the given node.
     *
     * @param pNode  Pointer to node object.
     */
    virtual AbstractAcinarUnit* CreateAcinarUnitForNode(Node<3>* pNode) = 0;

    /**
     * May be overridden by subclasses to perform any necessary work after all acinar
     * have been created.
     *
     * @param pAcinarUnitsDistributed  Pointer to a vector of acinar unit pointers.
     * @param lo  Lowest index owned by this process.
     * @param hi  Highest index owned by this process.
     */
//    virtual void FinaliseAcinarUnitCreation(std::vector< AbstractAcinarUnit* >* pAcinarUnitsDistributed,
//                                            unsigned lo, unsigned hi);

    /**
     * @return  The number of acini
     */
    virtual unsigned GetNumberOfAcini();

    /**
     * As this method is pure virtual, it must be overridden in subclasses.
     *
     * @param time The current time in seconds
     * @param pNode Pointer to node object.
     * @return The pleural pressure at the given node at the given time
     */
    virtual double GetPleuralPressureForNode(double time, Node<3>* pNode)=0;

    /**
     * Default constructor.
     */
    AbstractAcinarUnitFactory();

    /**
     * Destructor
     */
    virtual ~AbstractAcinarUnitFactory();

    /**
     * @param pMesh  the mesh for which to create acinar units.
     */
    virtual void SetMesh(AbstractTetrahedralMesh<1,3>* pMesh);

    /**
     * @return  the mesh used to create the acinar.
     */
    AbstractTetrahedralMesh<1,3>* GetMesh();
};

#endif /*ABSTRACTACINARUNITFACTORY_HPP_*/
