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
#ifndef _ABSTRACTLINEARPARABOLICPDE_HPP_
#define _ABSTRACTLINEARPARABOLICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <petscvec.h>


/**
 * AbstractLinearParabolicPde class.
 *
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractLinearParabolicPde
{
public:

    /**
     * The function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     *
     * @param rX the point in space at which the function c is computed
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * Compute source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     */
    virtual double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * Compute source term at a node.
     *
     * @param rNode the node at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the node
     */
    virtual double ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode, double u);


    /**
     * Compute diffusion term. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX The point in space at which the diffusion term is computed.
     * @param pElement The mesh element that x is contained in (optional).
     * @return A matrix.
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)=0;
    /**
     * Destructor.
     */
    virtual ~AbstractLinearParabolicPde()
    {}
};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearParabolicPde<ELEMENT_DIM, SPACE_DIM>::ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode, double u)
{
    return ComputeSourceTerm(rNode.GetPoint(), u);
}

#endif //_ABSTRACTLINEARPARABOLICPDE_HPP_
