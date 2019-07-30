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

#ifndef ABSTRACTCONTRACTIONCELLFACTORY_HPP_
#define ABSTRACTCONTRACTIONCELLFACTORY_HPP_

#include "AbstractContractionModel.hpp"
#include "QuadraticMesh.hpp"

/**
 * A factory to ease creating contraction cell models for use in a electro-mechanics simulations.
 *
 * The user should implement their own concrete class, in particular implementing
 * CreateContractionCellForElement(Element*), which should return the contraction model corresponding to a
 * given element. The user should also implement GetNumberOfCells() if this isn't equal
 * to the number of nodes. FinaliseCellCreation() can be used to (eg) add stimuli to
 * certain cells after they have been created.
 *
 * DIM is the dimension of the quadratic mesh (element dimension == space dimension for mechanics at present).
 */
template<unsigned DIM>
class AbstractContractionCellFactory
{
private:
    /** This class tests that the mesh is set correctly */
    friend class TestElectroMechanicsProblemDefinition;

protected:
    /** The mechanics mesh */
    QuadraticMesh<DIM>* mpMesh;

public:
    /**
     * @return a newly created contraction model for the given quadrature point.
     *
     * Note that there is no vector of all the quadrature points of the mesh;
     * the quad point index is the index that would be obtained by looping over
     * elements and then looping over quad points.
     *
     * @param pElement  Pointer to element.
     * It is assumed that all the quad points in a given element will be assigned the same kind
     * of contraction cell.
     *
     */
    virtual AbstractContractionModel* CreateContractionCellForElement(Element<DIM, DIM>* pElement) = 0;

    /**
     * Set the mechanics mesh to be used by this cell factory.
     *
     * @param pMesh  A quadratic (mechanics) mesh.
     */
    void SetMechanicsMesh(QuadraticMesh<DIM>* pMesh)
    {
        mpMesh = pMesh;
    }

    //
    // \todo #2370 The methods below are likely to be required for the final version of this class.
    // Thus the definitions have been left in place commented out
    // You'll also need to re-doxygen the comments!
    //

//    /*
//     * May be overridden by subclasses to perform any necessary work after all cells
//     * have been created.
//     *
//     * @param pCellsDistributed  Pointer to a vector of cardiac cell pointers.
//     * @param lo  Lowest index owned by this process.
//     * @param hi  Highest index owned by this process.
//     */
//    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCellInterface* >* pCellsDistributed,
//                                      unsigned lo, unsigned hi);

//    /*
//     * Default constructor.
//     *
//     * @param pSolver  the ODE solver to use to simulate this cell.
//     */
//    AbstractContractionCellFactory(boost::shared_ptr<AbstractIvpOdeSolver> pSolver = boost::shared_ptr<AbstractIvpOdeSolver>(new EulerIvpOdeSolver));

    /**
     * Destructor
     */
    virtual ~AbstractContractionCellFactory(){};
};

#endif /* ABSTRACTCONTRACTIONCELLFACTORY_HPP_ */
